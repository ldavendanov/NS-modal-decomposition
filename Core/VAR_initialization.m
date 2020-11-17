function [InitialValues,HyperPar] = VAR_initialization(y,M,InitialGuess)

if ~isfield(InitialGuess,'Variances')
    Sigma = [1e-5 1e-5];
else
    Sigma = InitialGuess.Variances;
end

d = size(y,1);                                                              % Dimensions of the input signal
if ~isfield(InitialGuess,'TargetFrequencies')
   
    na = ceil(2*M/d);                                                           % Model order for requested number of modes
    [A,R] = VARestimate(y,na);                                                  % Calculate Vector AR model for the given data
    
    % Eigenvalue decomposition of the associated SS representation
    [V,D] = eig([A; eye(d*(na-1),d*na)]);
    omega = angle(diag(D));                                                     % Natural frequencies
    [~,ind] = sort(abs(omega));
    
    % Initial value of the parameter vector
    theta0 = zeros(2*M,1);
    theta0(1:2:2*M) = cos(omega(ind(1:2:2*M)));
    theta0(2:2:2*M) = sin(omega(ind(1:2:2*M)));
    
    % Initial value of the state vector
    x = V\[y(:,2); y(:,1)];
    x0 = zeros(2*M,1);
    x0(1:2:end) = real( x(ind(1:2:2*M)) );
    x0(2:2:end) = imag( x(ind(1:2:2*M)) );
    
    % Packing output
    InitialValues.x0 = [x0; theta0];
    InitialValues.P0 = 1e-6*eye(4*M);
    HyperPar.Q = Sigma(2)*eye(4*M);
    HyperPar.Q(1:2*M,1:2*M) = Sigma(1)*eye(2*M);
    HyperPar.R = R;
    HyperPar.Psi = zeros(d,2*M);
    HyperPar.Psi(:,1:2:end) = real(V(1:d,ind(1:2:2*M)));
    HyperPar.Psi(:,2:2:end) = imag(V(1:d,ind(1:2:2*M)));
    
else
    
    na = ceil(4*M/d);                                                           % Model order for requested number of modes
    [A,R] = VARestimate(y,na);                                                  % Calculate Vector AR model for the given data
    
    % Eigenvalue decomposition of the associated SS representation
    [V,D] = eig([A; eye(d*(na-1),d*na)]);
    omega = angle(diag(D));                                                     % Natural frequencies
    
    % Finding the modes of the target frequencies
    ind = zeros(1,M);
    for i=1:M
        delta = omega - InitialGuess.TargetFrequencies(i);
        [~,ind(i)] = min(abs(delta));
    end
    
    % Initial value of the parameter vector
    theta0 = zeros(2*M,1);
    theta0(1:2:2*M) = cos(InitialGuess.TargetFrequencies);
    theta0(2:2:2*M) = sin(InitialGuess.TargetFrequencies);
    
    % Initial value of the state vector
%     x = V\[y(:,2); y(:,1)];
    x0 = zeros(2*M,1);
%     x0(1:2:end) = real( x(ind(1:2:2*M)) );
%     x0(2:2:end) = imag( x(ind(1:2:2*M)) );

    % Packing output
    InitialValues.x0 = [x0; theta0];
    InitialValues.P0 = 1e-6*eye(4*M);
    HyperPar.Q = Sigma(2)*eye(4*M);
    HyperPar.Q(1:2*M,1:2*M) = Sigma(1)*eye(2*M);
    HyperPar.R = R;
    HyperPar.Psi = zeros(d,2*M);
    HyperPar.Psi(:,1:2:end) = real(V(1:d,ind(1:M)));
    HyperPar.Psi(:,2:2:end) = imag(V(1:d,ind(1:M)));
    
end
