function [Modal,logMarginal,State,Covariances,Gain] = MO_DSS_JointEKF(y,M,estim,InitialValues,HyperPar)
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

% Default hyperparameters
switch nargin
    case 2
        [InitialValues,HyperPar] = VAR_initialization(y,M);
        estim = 'KF';
    case 3
        [InitialValues,HyperPar] = VAR_initialization(y,M);
    case 4
        [~,HyperPar] = VAR_initialization(y,M);
end

% Setting up the initial values
Initial.x0 = InitialValues.x0;
Initial.P0 = InitialValues.P0;

% Setting up the state space representation
System.ffun = @ffun;
System.F = @stm;
System.H = [HyperPar.Psi zeros(d,n)];

%% Pt 2 : Perform state estimation based on the provided input

switch estim
    case 'KF'       
        % -- Estimate with Kalman filter --
        [State,Covariances,logMarginal,K] = KalmanFilter( y, System, HyperPar, Initial );
        Gain.K = K;
        
    case 'KS'       
        % -- Estimate with Kalman filter and smoother --
        [State,Covariances,logMarginal,K] = KalmanFilter( y, System, HyperPar, Initial );
        [State,Covariances,J] = KalmanSmoother( State, Covariances, System );
        Gain.K = K; Gain.J = J;
        
    case 'KS+'
        % -- Estimate with Kalman filter and smoother, and then perform lag-one covariance smoother --
        [State,Covariances,logMarginal,K] = KalmanFilter( y, System, HyperPar, Initial );
        [State,Covariances,J] = KalmanSmoother( State, Covariances, System );
        Covariances = LagOneCovSmoother( State, Covariances, System, K, J );
        Gain.K = K; Gain.J = J;
        
end

%% Pt 3 : Packing the results of the modal decomposition

switch estim
    case 'KF'
        xhat = State.xtt;
        
    case {'KS','KS+'}
        xhat = State.xtN;

end

Modal.ym = xhat(1:2*M,:);
Modal.theta = xhat(2*M+1:4*M,:);
Modal.omega = zeros(M,N);
Modal.Am = zeros(M,N);
for m=1:M
    Modal.omega(m,:) = atan2( Modal.theta(2*m,:), Modal.theta(2*m-1,:) );
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

