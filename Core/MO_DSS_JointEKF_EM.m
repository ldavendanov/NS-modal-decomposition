function [Modal,logMarginal,HyperPar] = MO_DSS_JointEKF_EM(y,M,Niter)
%--------------------------------------------------------------------------
% Joint EKF estimator for Multiple-Output Diagonal State Space
% representation. This function estimates the 'M' modal components of a
% multicomponent non-stationary signal 'y' using the Joint EKF method.
%
% Created by : David Avendano - September 2020
%    Updated : David Avendano - November 2020
%--------------------------------------------------------------------------

%% Pt 1 : Initial set-up

n = 2*M;                    % Dimension of the modal and parameter vector
[d,N] = size(y);            % Size of the signal

% Set up initial values
% [Initial,HyperPar] = VAR_initialization(y(:,1:500*M),M);
[Initial,HyperPar] = STFT_initialization(y,M);

% Setting up the state space representation
System.ffun = @ffun;
System.F = @stm;
System.H = [HyperPar.Psi zeros(d,n)];

%% Pt 2 : Expectation-Maximization algorithm for state estimation and hyperparameter estimation

ind = 1:n;          % Indices of the state vector
logMarginal = zeros(1,Niter);

for k=1:Niter
    
    fprintf('EM algorithm - Iteration No. %3d ',k)
    
    % -- P2.1 : Expectation step - State estimation --
    [State,Covariances,logM,K] = KalmanFilter( y, System, HyperPar, Initial );
    [State,Covariances,J] = KalmanSmoother( State, Covariances, System );
    Covariances = LagOneCovSmoother( State, Covariances, System, K, J );
    logMarginal(k) = mean(logM);
    
    fprintf('Log-likelihood No. %2.4f\n',logMarginal(k))
    
    % -- P2.2 : Maximization step - Updating hyperparameters --
    
    % P2.2.1 : Update initial values
    Initial.x0 = State.xtN(:,1);
    Initial.P0 = Covariances.PtN(:,:,1);
    
    % P2.2.2 : Update state noise covariance
    % Individual covariances
    S11 = State.xtN(:,2:N)*State.xtN(:,2:N)' + sum( Covariances.PtN(:,:,2:N), 3 );
    S10 = State.xtN(:,2:N)*State.xtN(:,1:N-1)' + sum( Covariances.PttmN(:,:,2:N), 3 );
    S00 = State.xtN(:,1:N-1)*State.xtN(:,1:N-1)' + sum( Covariances.PtN(:,:,1:N-1), 3 );
    
    % Covariance update
    HyperPar.Q = ( S11 -  ( S10/S00 )*S10' )/(N-1);     % Covariance update equation
    HyperPar.Q(ind,ind+n) = 0;                          % Making sure that state and parameters are uncorrelated
    HyperPar.Q(ind+n,ind) = 0;
    
    % P2.2.3 : Update measurement noise covariance
    etN = y - System.H*State.xtN;
    HyperPar.R = etN*etN'/N + 1e-5*eye(d);
    for i=1:N
        HyperPar.R = HyperPar.R + ( System.H*Covariances.PtN(:,:,i)*System.H' )/N;
    end
    
    % P2.3.4 : Update mixing matrix (state measurement matrix)
    Psi = ( y(:,2:N)*State.xtN(ind,2:N)' ) / S11(ind,ind);
    System.H = [Psi zeros(d,n)];    
    
end

%% Pt 3 : Packing the results of the modal decomposition

xhat = State.xtN;

Modal.ym = xhat(1:2*M,:);
Modal.theta = xhat(2*M+1:4*M,:);
Modal.omega = zeros(M,N);
Modal.zeta = zeros(M,N);
Modal.Am = zeros(M,N);
for m=1:M
    Modal.omega(m,:) = atan2( Modal.theta(2*m,:), Modal.theta(2*m-1,:) );
    Modal.zeta(m,:) = -cos( angle( -Modal.theta(2*m,:) + 1i*Modal.theta(2*m-1,:) ) );
    Modal.Am(m,:) = sqrt( xhat(2*m-1,:).^2 + xhat(2*m,:).^2 );
end

%--------------------------------------------------------------------------
function z_new = ffun(z_old)

n = length(z_old);
M = dffun_dz(z_old(n/2+1:n));
z_new = z_old;
z_new(1:n/2) = M*z_old(1:n/2);

%--------------------------------------------------------------------------
function F = stm(z)

n = length(z);
M = dffun_dz(z(n/2+1:n));
Z = dffun_dtheta(z(1:n/2));
F = [M Z; zeros(n/2) eye(n/2)];

%--------------------------------------------------------------------------
function M = dffun_dz(theta)

n = length(theta);
M = zeros(n);
for k=1:n/2
    M(2*k-1,2*k-1) =  theta(2*k-1);
    M(2*k  ,2*k  ) =  theta(2*k-1);
    M(2*k-1,2*k  ) =  theta(2*k);
    M(2*k  ,2*k-1) = -theta(2*k);
end

%--------------------------------------------------------------------------
function Z = dffun_dtheta(z)

n = length(z);
Z = zeros(n);
for k=1:n/2
    Z(2*k-1,2*k-1) =  z(2*k-1);
    Z(2*k  ,2*k  ) = -z(2*k-1);
    Z(2*k-1,2*k  ) =  z(2*k);
    Z(2*k  ,2*k-1) =  z(2*k);
end

