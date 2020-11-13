function [InitialValues,HyperPar] = VAR_initialization(y,M)

d = size(y,1);                                                              % Dimensions of the input signal
na = ceil(2*M/d);                                                           % Model order for requested number of modes
[A,R] = VARestimate(y,na);                                                  % Calculate Vector AR model for the given data

% Eigenvalue decomposition of the associated SS representation
[V,D] = eig([A; eye(d*(na-1),d*na)]);

% Initial value of the parameter vector
omega = angle(diag(D));                                                     % Natural frequencies
[omega,ind] = sort(abs(omega));
theta0 = zeros(2*M,1);
theta0(1:2:2*M) = cos(omega(1:2:2*M));
theta0(2:2:2*M) = sin(omega(1:2:2*M));

% Initial value of the state vector
x = V(1:d,ind)'*y(:,1);
x0 = zeros(2*M,1);
x0(1:2:end) = real( x(1:2:2*M) );
x0(2:2:end) = imag( x(1:2:2*M) );

% Packing output
InitialValues.x0 = [x0; theta0];
InitialValues.P0 = 1e-8*eye(4*M);
HyperPar.Q = 1e-4*eye(4*M);
HyperPar.R = R;
HyperPar.Psi = zeros(d,2*M);
HyperPar.Psi(:,1:2:end) = real(V(1:d,ind(1:2:2*M)));
HyperPar.Psi(:,2:2:end) = imag(V(1:d,ind(1:2:2*M)));