function [InitialValues,HyperPar] = FFT_initialization(y,M,omega)

[d,N] = size(y);                                                            % Dimensions of the input signal

Y = fft(y,[],2)/N;
ff = 2*pi*(0:N-1)/N;

if nargin < 3
    
    % Finding peaks in the FFT
    [pks,locs] = findpeaks(max(So,[],2),'MinPeakDistance',10);
    [~,ind] = sort(pks,'descend');
    locs = locs(ind(1:M));
    locs = sort(locs,'ascend');
    omega = ff(locs);                                                           % Target frequencies
    
else
    
    locs = zeros(M,1);
    for i=1:numel(omega)
        dOm = abs( ff - omega(i) );
        [~,locs(i)] = min(dOm);
    end
    
end

% Initial value of the parameter vector
theta0 = zeros(2*M,1);
theta0(1:2:2*M) = cos(omega);
theta0(2:2:2*M) = sin(omega);

% Initial value of the state vector
x0 = zeros(2*M,1);
for i=1:M
    x0(2*i-1) = 1;
    x0(2*i) = 0;
end

Psi = zeros(d,2*M);
Psi(:,1:2:2*M) = real(Y(:,locs));
Psi(:,2:2:2*M) = imag(Y(:,locs));

% Packing output
InitialValues.x0 = [x0; theta0];
InitialValues.P0 = 1e-8*eye(4*M);
HyperPar.Q = 1e-4*eye(4*M);
HyperPar.R = diag(var(y,[],2));
HyperPar.Psi = Psi;
