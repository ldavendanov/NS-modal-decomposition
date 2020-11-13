function [InitialValues,HyperPar] = STFT_initialization(y,M)

[d,N] = size(y);                                                            % Dimensions of the input signal
Nfft = 512;                                                                 % Number of points of the STFT
delta = 2;
Nover = Nfft - delta;

Sxx = zeros(Nfft/2+1,floor((N-Nover)/delta),d);
for i=1:d
    [Sxx(:,:,i),ff] = spectrogram( y(i,:), gausswin(Nfft,6), Nover, Nfft );
end

So = max( mean(abs(Sxx),2), [], 3 );

[pks,locs] = findpeaks(So,'MinPeakDistance',10);
[~,ind] = sort(pks,'descend');
locs = locs(ind(1:M));
locs = sort(locs,'ascend');

% Initial value of the parameter vector
omega = ff(locs);                                                           % Target natural frequencies
theta0 = zeros(2*M,1);
theta0(1:2:2*M) = cos(omega);
theta0(2:2:2*M) = sin(omega);

% Initial value of the state vector
Psi = zeros(d,2*M);
for i=1:M
    g = [cos(omega(i)*(0:N-1));
         sin(omega(i)*(0:N-1))];
     
    Psi(:,(1:2)+2*(i-1)) = y*g'/N;
     
end
x0 = Psi'*y(:,1);

% Packing output
InitialValues.x0 = [x0; theta0];
InitialValues.P0 = 1e-8*eye(4*M);
HyperPar.Q = 1e-4*eye(4*M);
HyperPar.R = diag(var(y,[],2));
HyperPar.Psi = Psi;


