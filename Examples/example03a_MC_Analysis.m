clear
close all
clc

addpath('..\Core\')
load('OptimizedHyperParams.mat','HyperPar')

%% Signal features

%-- Signal features
N = 4e3;                                                                    % Number of samples
fs = 5e2;                                                                   % Sampling frequency
t = 0:N-1;                                                                  % Time vector

Aii = [1 0.2; 0.8 0.4; -0.2 1.0];                                           % Amplitude coefficients
fii = [50 -20; 80 20; 120 20];                                              % Frequency coefficients
alpha = [0.5 0.5 0.25];                                                     % Cyclic frequency
Psi =  [1.0 0.8 0.5 0.5 0.2 0.8;
        0.5 0.5 0.8 1.0 0.4 0.6;
        0.4 0.9 0.1 0.2 1.0 0.0];
d = size(Psi,1);

for i=1:d
    Psi(i,:) = Psi(i,:)/norm(Psi(i,:));
end

%-- Instantaneous Amplitudes
IA = zeros(3,N);
IA(1,:) = Aii(1,1) + Aii(1,2)*cos(2*pi*alpha(1)*t/fs);
IA(2,:) = Aii(2,1) + Aii(2,2)*sin(2*pi*alpha(2)*t/fs);
IA(3,:) = max( Aii(3,1) + Aii(3,2)*abs( cos(2*pi*alpha(3)*t/fs) ), 0 );

%-- Instantaneous frequencies
IF = zeros(3,N);
IF(1,:) = fii(1,1) + fii(1,2)*sin( 2*pi*alpha(1)*t/fs );
IF(2,:) = fii(2,1) + fii(2,2)*sin( 2*pi*alpha(2)*t/fs );
IF(3,:) = fii(3,1) + fii(3,2)*sin( 2*pi*alpha(3)*t/fs );

%% Monte Carlo analysis
close all
clc

Nmc = 40;
snr = 0:10:60;
Nsnr = numel(snr);

% Initialization
M = 3;
IniGuess.Variances = [1e-4 1e-5];
logMarginal = zeros(Nmc,Nsnr,2);
NMSE = zeros(Nmc,M,Nsnr,2);
NMSE_IA = zeros(Nmc,M,Nsnr,2);
NMSE_IF = zeros(Nmc,M,Nsnr,2);

for j=1:Nsnr
    
    fprintf('SNR = %2d dB\n',snr(j))
    fprintf('Evaluating signal No. %2d\n',0)
    
    for k = 1:Nmc
        
        fprintf('\b\b\b%2d\n',k)
        
        % Pt.1 : Creating a signal realization --------------------------------
        
        % -- Modes
        ym = zeros(6,N);
        Phi0 = 2*pi*rand(3,1);
        for i=1:3
            ym(2*i-1,:) = IA(i,:).*cos( cumsum( 2*pi*IF(i,:)/fs ) + Phi0(i) );
            ym(2*i,:) = IA(i,:).*sin( cumsum( 2*pi*IF(i,:)/fs ) + Phi0(i) );
        end
        C = var(ym(2:2:end,:),[],2);
        
        % -- Measured signals
        y = Psi*ym;
        sigmaN = sqrt( min(var(y,[],2))*10.^( -snr(j)/10 ) );
        y = y + sigmaN*randn(d,N);
        
        % Pt.2 : Estimating the signal components with the diagonal SS method -
        
        % -- Estimate with manually set hyperparameters
        [InitialValues,HPar] = VAR_initialization(y(:,1:200*M),M,IniGuess);
        [Modal{1},lnM] = MO_DSS_JointEKF(y,M,'KS',InitialValues,HPar);
        logMarginal(k,j,1) = mean(lnM);
        
        NMSE(k,:,j,1) = var( HPar.Psi*Modal{1}.ym - y, [], 2 )./var( y, [], 2 );
        NMSE_IA(k,:,j,1) = var( Modal{1}.Am - IA./repmat(C,1,N), [], 2 )./var( IA./repmat(C,1,N), [], 2 );
        NMSE_IF(k,:,j,1) = var( Modal{1}.omega*fs/(2*pi) - IF, [], 2 )./var( IF, [], 2 );
        
        % -- Estimate EM-optimized hyperparameters
        HyperPar.R = HPar.R;
        [Modal{2},lnM] = MO_DSS_JointEKF(y,M,'KS',InitialValues,HyperPar);
        logMarginal(k,j,2) = mean(lnM);
        
        NMSE(k,:,j,2) = var( HyperPar.Psi*Modal{2}.ym - y, [], 2 )./var( y, [], 2 );
        NMSE_IA(k,:,j,2) = var( Modal{2}.Am - IA./repmat(C,1,N), [], 2 )./var( IA./repmat(C,1,N), [], 2 );
        NMSE_IF(k,:,j,2) = var( Modal{2}.omega*fs/(2*pi) - IF, [], 2 )./var( IF, [], 2 );
        
    end
end

%%
close all
clc

FName = 'Times New Roman';
FSize = 12;

figure('Position',[100 100 600 640])
for i=1:2
    subplot(311)
    semilogy(snr,squeeze(median( min( NMSE(:,:,:,i), [], 2 ) ))*100,'LineWidth',1.5)
    hold on
    grid on
    ylim([1e-1 1e2])
    set(gca,'YTick',10.^(-1:2))
    ylabel('NMSE [%]')
    set(gca,'FontName',FName,'FontSize',FSize)
    
    subplot(312)
    semilogy(snr,squeeze(median( min( NMSE_IA(:,:,:,i), [], 2 )))*100,'LineWidth',1.5)
    hold on
    grid on
    ylim([5e0 5e2])
    set(gca,'YTick',10.^(1:2))
    ylabel('NMSE - IA [%]')
    set(gca,'FontName',FName,'FontSize',FSize)
    
    subplot(313)
    semilogy(snr,squeeze(median( min( NMSE_IF(:,1,:,i), [], 2 )))*100,'LineWidth',1.5)
    hold on
    grid on
    ylim([1e-2 1e2])
    set(gca,'YTick',10.^(-2:2:2))
    ylabel('NMSE - IF [%]')
    set(gca,'FontName',FName,'FontSize',FSize)
end