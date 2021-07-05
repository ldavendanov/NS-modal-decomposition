function [InitialValues,HyperPar] = STFT_initialization(y,M,omega)

[d,N] = size(y);                                                            % Dimensions of the input signal

Nfft = 512;                                                             % Number of points of the STFT
delta = 2;                                                              % Time jump between samples
Nover = Nfft - delta;                                                   % Number of overlapping samples

% Calculating the spectrogram
Sxx = zeros(Nfft/2+1,floor((N-Nover)/delta),d);
for i=1:d
    [Sxx(:,:,i),ff] = spectrogram( y(i,:), gausswin(Nfft,6), Nover, Nfft );
end

% Calculating the mean spectrum
So = squeeze( mean(abs(Sxx),2) );

if nargin < 3
    
    % Finding peaks in the mean spectrum
    [pks,locs] = findpeaks(max(So,[],2),'MinPeakDistance',10);
    [~,ind] = sort(pks,'descend');
    locs = locs(ind(1:M));
    locs = sort(locs,'ascend');
    omega = ff(locs);                                                           % Target natural frequencies
    
else
    
    locs = zeros(M,1);
    for i=1:M
        locs(i) = find( ff > omega(i), 1, 'first' );
    end
    
end

% Initial value of the parameter vector
theta0 = zeros(2*M,1);
theta0(1:2:2*M) = cos(omega);
theta0(2:2:2*M) = sin(omega);

% Initial value of the state vector
Psi = zeros(d,2*M);
for i=1:M
    
    Psi(:,(1:2)+2*(i-1)) = [real(So(locs(i),:))' imag(So(locs(i),:))'];
     
end

x0 = zeros(2*M,1);

% Packing output
InitialValues.x0 = [x0; theta0];
InitialValues.P0 = 1e-8*eye(4*M);
HyperPar.Q = 1e-4*eye(4*M);
HyperPar.R = diag(var(y,[],2));
HyperPar.Psi = Psi;
