function [Modal,HyperPar,State,Covariances,logMarginal] = DiagonalMvSS1(y,M,InitialValues,estim)
%--------------------------------------------------------------------------
% Estimation of the modal components of a non-stationary signal using the
% Joint EKF method 
% Created by : David Avendano - September 2020
%    Updated : David Avendano - November 2020
%--------------------------------------------------------------------------

%% Pt 0 : Initial set-up

n = 2*M;                    % Dimension of the modal and parameter vector
[d,N] = size(y);            % Size of the signal

% Setting up the state space representation
System.ffun = @ffun;
System.F = @stm;
System.H = [InitialValues.Psi zeros(d,n)];

% Setting up the initial guess of the noise structure
if nargin < 3
    CovStructure.R = eye(d);
    lambda = 1e-4;
    mu = 5e2*lambda;
    estim = 'KF';
    
elseif nargin < 4
    CovStructure.R = InitialValues.Variances(1)*eye(d);                     % Measurement noise variance
    lambda = InitialValues.Variances(2);                                    % Variance of the parameters
    mu = InitialValues.Variances(3);                                        % Variance of the mode trajectories
    estim = 'KF';
    
else
    CovStructure.R = InitialValues.Variances(1)*eye(d);                     % Measurement noise variance
    lambda = InitialValues.Variances(2);                                    % Variance of the parameters
    mu = InitialValues.Variances(3);                                        % Variance of the mode trajectories

end
CovStructure.Q = lambda*eye(2*n);
CovStructure.Q(1:n,1:n) = mu*eye(n);

% Setting up the initial values
Initial.x0 = InitialValues.x0;
Initial.P0 = 1e-8*eye(2*n);

%% Pt 1 : First run

% -- Kalman filter
[State,Covariances,logMarginal,K] = KalmanFilter( y, System, CovStructure, Initial );

% -- Kalman smoother
if strcmp( estim, 'KS' )
    [State,Covariances,J] = KalmanSmoother( State, Covariances, System );
    Covariances = LagOneCovSmoother( State, Covariances, System, K, J );
end

%% Pt 2 : Hyperparameter adjustment with Expectation-Maximization 

for k=1:10
    
    fprintf('EM iteration No. %2d - Log marginal %2.4f\n',k,mean(logMarginal))
    
    %-- 2.1 Updating hyperparameters (maximization step)
    % Initial values
    Initial.x0 = State.xtN(:,1);
    Initial.P0 = Covariances.PtN(:,:,1);
    
    % Covariances
    ind = 1:n;
    etN = y - System.H*State.xtN;
    CovStructure.R = etN*etN'/N;
    
    S11 = State.xtN(:,2:N)*State.xtN(:,2:N)' + sum( Covariances.PtN(:,:,2:N), 3 );
    S10 = State.xtN(:,2:N)*State.xtN(:,1:N-1)' + sum( Covariances.PttmN(:,:,2:N), 3 );
    S00 = State.xtN(:,1:N-1)*State.xtN(:,1:N-1)' + sum( Covariances.PtN(:,:,1:N-1), 3 );
    
    for i=1:N
        
        CovStructure.R = CovStructure.R + ( System.H*Covariances.PtN(:,:,i)*System.H' )/N;
        
    end
    
    Psi = y(:,2:N)*State.xtN(ind,2:N)'/S11(ind,ind);
    System.H = [Psi zeros(d,n)];
    CovStructure.Q = ( S11 -  ( S10/S00 )*S10' )/N;
        
    %-- 2.2 State vector estimates (expectation)
    [State,Covariances,logMarginal,K] = KalmanFilter( y, System, CovStructure, Initial );
    [State,Covariances,J] = KalmanSmoother( State, Covariances, System );
    Covariances = LagOneCovSmoother( State, Covariances, System, K, J );
    
end

HyperPar.Psi = Psi;
HyperPar.Q = CovStructure.Q;
HyperPar.R = CovStructure.R;

%% Pt 3 : Packing the results of the modal decomposition
switch estim
    case 'KF'
        xhat = State.xtt;
        
    case 'KS'
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
function H = smm(C)

[d,M] = size(C);
H = [C zeros(d,M)];

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