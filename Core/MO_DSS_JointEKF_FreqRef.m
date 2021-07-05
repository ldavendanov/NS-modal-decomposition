function [Modal,logMarginal,HyperPar,Initial] = MO_DSS_JointEKF_FreqRef(y,omega_ref,Orders,HyperPar,Initial)
%--------------------------------------------------------------------------
% Joint EKF estimator for Multiple-Output Diagonal State Space
% representation. This function estimates the 'M' modal components of a
% multicomponent non-stationary signal 'y' using the Joint EKF method.
%
% Created by : David Avendano - September 2020
%    Updated : David Avendano - November 2020
%    Updated : David Avendano - April 2021
%--------------------------------------------------------------------------

%% Pt 1 : Initial set-up

M = numel(Orders);
n = 2*M;                    % Dimension of the modal and parameter vector
[d,N] = size(y);            % Size of the signal

% Set up initial values
if nargin<4
    [Initial,HyperPar] = STFT_initialization(y,M,Orders*InitialGuess.TargetFrequency);
    Initial.x0 = [Initial.x0(1:2*M); omega_ref(1)];
    Initial.P0 = Initial.P0(1:2*M+1,1:2*M+1);
    HyperPar.Q = InitialGuess.Variances(1)*eye(2*M+1);
    HyperPar.Q(end,end) = InitialGuess.Variances(2);
    HyperPar.R(d+1,d+1) = 1e-4;
end

% Setting up the state space representation
System.ffun = @(z)ffun(z,Orders);
System.F = @(z)stm(z,Orders);
System.H = [HyperPar.Psi zeros(d,1); zeros(1,2*M) 1];

%% Pt 2 : Expectation-Maximization algorithm for state estimation and hyperparameter estimation

y = [y; omega_ref]; % Extending output vector
T = 2:N;

% -- P2.1 : Expectation step - State estimation --
[State,Covariances,logM,K] = KalmanFilter( y, System, HyperPar, Initial );
[State,Covariances,J] = KalmanSmoother( State, Covariances, System );
logMarginal = mean(logM(T));

fprintf('Log-likelihood %2.4f\n',logMarginal)

%% Pt 3 : Packing the results of the modal decomposition

xhat = State.xtN;

Modal.ym = xhat(1:2*M,:);
Modal.theta = xhat(2*M+1,:);
Modal.omega = xhat(2*M+1,:);
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


%--------------------------------------------------------------------------
function z_new = ffun(z_old,ord)

n = length(z_old)-1;
M = dffun_dz(z_old(end),ord);
z_new = z_old;
z_new(1:n) = M*z_old(1:n);

%--------------------------------------------------------------------------
function F = stm(z,ord)

n = length(z)-1;
M = dffun_dz(z(end),ord);
Z = dffun_dtheta(z(1:n),z(end),ord);
F = [M Z; zeros(1,n) 1];

%--------------------------------------------------------------------------
function M = dffun_dz(theta,ord)

n = length(ord);
M = eye(2*n);
for k=1:n
    ind = (1:2)+2*(k-1);
    M(ind,ind) = [ cos(ord(k)*theta) sin(ord(k)*theta)
                  -sin(ord(k)*theta) cos(ord(k)*theta)];
end

%--------------------------------------------------------------------------
function Z = dffun_dtheta(z,theta,ord)

n = length(ord);
Z = zeros(2*n,1);
for k=1:n-1
    ind = (1:2)+2*(k-1);
    Z(ind) = ord(k)*[-sin(ord(k)*theta)  cos(ord(k)*theta);
                     -cos(ord(k)*theta) -sin(ord(k)*theta)]*z(ind);
end
