function [State,Covariances,logLikelihood,K] = KalmanFilter( y, System, CovStructure, Initial )

% Extracting information from the input
[d,N] = size(y);                                        % Dimension of the measurements
n = numel( Initial.x0 );                                % Dimension of the state vector

% Initializing the estimation matrices
xtt = zeros(n,N);       xttm = zeros(n,N);
Ptt = zeros(n,n,N);     Pttm = zeros(n,n,N);
ett = zeros(d,N);       ettm = zeros(d,N);
K = zeros(n,d,N);
logLikelihood = zeros(1,N);

% Setting up the initial values
xtt(:,1) = Initial.x0;
Ptt(:,:,1) = Initial.P0;

%-- Kalman filter ---------------------------------------------------------
H = System.H;                                           % State measurement matrix
for i=2:N
    
    % System matrices
    F = System.F( xtt(:,i-1) );                         % State transition matrix
    
    % Update equations
    xttm(:,i) = System.ffun( xtt(:,i-1) );
    ettm(:,i) = y(:,i) - H*xttm(:,i);
    Pttm(:,:,i) = F*Ptt(:,:,i-1)*F' + CovStructure.Q;
    
    % Calculate the Kalman gain
    Sigma_e2 = H*Pttm(:,:,i)*H' + CovStructure.R;
    logLikelihood(i) = -0.5*( log( det(Sigma_e2) ) + ettm(:,i)'*( Sigma_e2\ettm(:,i) ) );
    K(:,:,i) = Pttm(:,:,i)*H' / Sigma_e2;
    
    % Correction equations
    xtt(:,i) = xttm(:,i) + K(:,:,i) * ( y(:,i) - H*xttm(:,i)  );
    ett(:,i) = y(:,i) - H*xtt(:,i);
    Ptt(:,:,i) = ( eye(n) - K(:,:,i)*H )*Pttm(:,:,i);

end

% Packing output
State.xtt = xtt;
State.xttm = xttm;

Covariances.Ptt = Ptt;
Covariances.Pttm = Pttm;
