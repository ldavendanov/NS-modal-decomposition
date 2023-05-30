function [Modal,logMarginal,HyperPar,Initial] = MO_DSS_JointEKF_MultiHar_EM(y,Orders,Niter,InitialGuess)
%--------------------------------------------------------------------------
% Joint EKF estimator for Multiple-Output Diagonal State Space
% representation. This function estimates the 'M' modal components of a
% multicomponent non-stationary signal 'y' using the Joint EKF method.
%
% Created by : David Avendano - September 2020
%    Updated : David Avendano - November 2020
%--------------------------------------------------------------------------

if nargin < 3
    Niter = 20;
    InitialGuess = [];
elseif nargin < 4
        InitialGuess = [];
end

%% Pt 1 : Initial set-up

[d,N] = size(y);                                                            % Size of the signal
p = numel(InitialGuess.TargetFrequency);                                    % Number of frequencies to track
M = numel(Orders)*p;                                                        % Total number of modes to track
n = 2*M;                                                                    % Dimension of the modal vector
ind = false(n+p,1);
ind(1:n) = true;          % Indices of the state vector

% Set up initial values
[Initial,HyperPar] = FFT_initialization(y(:,1:512),M,kron(Orders(:),InitialGuess.TargetFrequency));

Initial.x0 = [Initial.x0(ind); InitialGuess.TargetFrequency];
Initial.P0 = Initial.P0(1:2*M+p,1:2*M+p);

HyperPar.Q = blkdiag( InitialGuess.Variances(2)*eye(n), InitialGuess.Variances(3)*eye(p) );
HyperPar.R = InitialGuess.Variances(1)*eye(d);

% Setting up the state space representation
System.ffun = @(z)ffun(z,Orders,p);
System.F = @(z)stm(z,Orders,p);
System.H = [HyperPar.Psi zeros(d,p)];

%% Pt 2 : Expectation-Maximization algorithm for state estimation and hyperparameter estimation

logMarginal = zeros(1,Niter);
T = 2:N;

figure('Position',[100 300 1200 400])

for k=1:Niter
    
    fprintf('EM algorithm - Iteration No. %3d - ',k)
        
    % -- P2.1 : Expectation step - State estimation --
    [State,Covariances,logM,K] = KalmanFilter( y, System, HyperPar, Initial );
    [State,Covariances,J] = KalmanSmoother( State, Covariances, System );
    Covariances = LagOneCovSmoother( State, Covariances, System, K, J );
    logMarginal(k) = mean(logM(T));
    
    fprintf('Log-likelihood %2.4f\n',logMarginal(k))
    
    % -- P2.2 : Maximization step - Updating hyperparameters --
    
    % P2.2.1 : Update initial values
    Initial.x0 = State.xtN(:,1);
    Initial.x0(~ind) = mode(State.xtN(~ind,:),2);
    Initial.P0 = Covariances.PtN(:,:,1);
    
    % P2.2.2 : Update noise covariances
    
    S11th = State.xtN(~ind,T)*State.xtN(~ind,T)' + sum( Covariances.PtN(~ind,~ind,T), 3 );
    S11 = State.xtN(ind,T)*State.xtN(ind,T)' + sum( Covariances.PtN(ind,ind,T), 3 );
    S10th = 0;   S00th = 0;
    S10 = zeros(n);     S00 = zeros(n);
    etN = y - System.H*State.xtN;
    HyperPar.R = etN(:,T)*etN(:,T)'/numel(T);
    
    for i=T
        S10th = S10th + ( State.xtN(~ind,i)*State.xtN(~ind,i-1)' ...
                        + Covariances.PttmN(~ind,~ind,i) );
        S00th = S00th + ( State.xtN(~ind,i-1)*State.xtN(~ind,i-1)' ...
                        + Covariances.PtN(~ind,~ind,i-1) );
                    
        F = System.F(State.xtN(:,i-1));
        F = F(ind,ind);
        S10 = S10 + ( State.xtN(ind,i)*State.xtN(ind,i-1)' ...
                        + Covariances.PttmN(ind,ind,i) )*F';
        S00 = S00 + F*( State.xtN(ind,i-1)*State.xtN(ind,i-1)' ...
                        + Covariances.PtN(ind,ind,i-1) )*F';
                    
        HyperPar.R = HyperPar.R + ( System.H*Covariances.PtN(:,:,i)*System.H' )/numel(T);
        
    end
    
    % Parameter covariance
    Qup_th = ( S11th - S10th - S10th' + S00th )/numel(T);
    HyperPar.Q(~ind,~ind) = Qup_th;
        
    % State covariance
    Qup = ( S11 - S10 - S10' + S00 )/numel(T);
    for i=1:M
        ind1 = (1:2)+2*(i-1);
        HyperPar.Q(ind1,ind1) = mean(diag(Qup(ind1,ind1)))*eye(2);
    end
    
    % P2.3.4 : Update mixing matrix (state measurement matrix)
    HyperPar.Psi = ( y(:,T)*State.xtN(ind,T)' ) / S11;
    System.H = [HyperPar.Psi zeros(d,p)];
    
    subplot(131)
    plot(k,logMarginal(k),'.b')
    hold on
    xlim([0 Niter])
    ylabel('Log marginal likelihood')
    xlabel('Iteration')
    
    subplot(132)
    semilogy(k,diag(HyperPar.Q(ind,ind)),'ob')
    hold on
    xlim([0 Niter])
    ylabel('State innov. variance')
    xlabel('State index')
    
    subplot(133)
    semilogy(k,diag(Qup_th),'ob')
    hold on
    xlim([0 Niter])
    ylabel('Param. innov. variance')
    xlabel('Param. index')
    
    drawnow
    
end

%% Pt 3 : Packing the results of the modal decomposition

xhat = State.xtN;

Modal.ym = xhat(1:2*M,:);
Modal.theta = xhat(2*M+(1:p),:);
Modal.omega = xhat(2*M+(1:p),:);
Modal.Am = zeros(M,N);
for m=1:M
    Modal.Am(m,:) = sqrt( xhat(2*m-1,:).^2 + xhat(2*m,:).^2 );
end

% Making sure that the modes are normalized
C = zeros(M,1);
for m=1:M
    C(m) = std(Modal.ym(2*m,:));
    Modal.ym((1:2)+2*(m-1),:) = Modal.ym((1:2)+2*(m-1),:)/C(m);
    Modal.Am(m,:) = Modal.Am(m,:)/C(m);
    HyperPar.Psi(:,(1:2)+2*(m-1)) = C(m)*HyperPar.Psi(:,(1:2)+2*(m-1));
end

HyperPar.p = p;

%--------------------------------------------------------------------------
function z_new = ffun(z_old,ord,p)

n = length(z_old)-p;
M = dffun_dz(z_old(n+(1:p)),ord,p);
z_new = z_old;
z_new(1:n) = M*z_old(1:n);

%--------------------------------------------------------------------------
function F = stm(z,ord,p)

n = length(z)-p;
M = dffun_dz(z(n+(1:p)),ord,p);
Z = dffun_dtheta(z(1:n),z(n+(1:p)),ord,p);
F = [M Z; zeros(p,n) eye(p)];

%--------------------------------------------------------------------------
function M = dffun_dz(theta,ord,p)

n = length(ord);
M = zeros(2*n*p);

for j=1:p
    for k=1:n
        ind = (1:2)  +2*(k-1) + 2*n*(j-1);
        M(ind,ind) = [ cos(ord(k)*theta(j)) sin(ord(k)*theta(j))
                      -sin(ord(k)*theta(j)) cos(ord(k)*theta(j))];
    end
end

%--------------------------------------------------------------------------
function Z = dffun_dtheta(z,theta,ord,p)

n = length(ord);
Z = zeros(2*n*p,p);

for j=1:p
    for k=1:n
        ind = (1:2) + 2*(k-1) + 2*n*(j-1);
        Z(ind,j) = ord(k)*[-sin(ord(k)*theta(j))  cos(ord(k)*theta(j));
                           -cos(ord(k)*theta(j)) -sin(ord(k)*theta(j))]*z(ind);
    end
end