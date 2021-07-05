clear
close all
clc

addpath('..\Core\')

% Pt.1 : Creating the simulation signal
N = 4e3;                                                                    % Number of samples
fs = 5e2;                                                                   % Sampling frequency
SNR = 5;                                                                   % Signal to noise ratio
m = 3;
Nmc = 100;

Methods = {'MO-Manual','MO-EM','SO-EM'};

IAhat = zeros(m,N,Nmc,3);
IFhat = zeros(m,N,Nmc,3);
yhat = zeros(m,N,Nmc,3);
y0 = zeros(m,N,Nmc);

fprintf('Simulation No. %2d\n',0)
for i = 1:Nmc
    
    fprintf('\b\b\b%2d\n',i)
    
    %-- Signal features
    [t,y,ym,IA,IF,Psi] = SignalGenerator(N,fs,SNR);
    y0(:,:,i) = Psi*ym;
    
    % Pt.2 : Estimating the signal components with the diagonal SS method
    
    % Initialization
    M = 3;
    IniGuess.Variances = [1e-4 1e-5];
    IniGuess.Freqs = 2*pi*( [50 80 120] + 4*randn(1,3) )/fs;
    
    % Estimate without EM optimization
    [InitialValues,HPar] = VAR_initialization(y(:,1:200*M),M,IniGuess);
    Modal = MO_DSS_JointEKF(y,M,'KS',InitialValues,HPar);
    IAhat(:,:,i,1) = Modal.Am;
    IFhat(:,:,i,1) = Modal.omega*fs/(2*pi);
    yhat(:,:,i,1) = HPar.Psi*Modal.ym;
    
    % Estimate after EM optimization - Multiple Output
    load(['OptimizedHyperParams_SNR',num2str(SNR)],'HyperPar')
    Modal = MO_DSS_JointEKF(y,M,'KS',InitialValues,HyperPar);
    IAhat(:,:,i,2) = Modal.Am;
    IFhat(:,:,i,2) = Modal.omega*fs/(2*pi);
    yhat(:,:,i,2) = HyperPar.Psi*Modal.ym;
    
    % Estimate after EM optimization - Single Output
    load(['OptimizedHyperParamsSO_SNR',num2str(SNR)],'HyperParSO')
    Modal = MO_DSS_JointEKF(y(1,:),M,'KS',InitialValues,HyperParSO);
    IAhat(:,:,i,3) = Modal.Am;
    IFhat(:,:,i,3) = Modal.omega*fs/(2*pi);
    yhat(1,:,i,3) = HyperParSO.Psi*Modal.ym;
    
end

%%
close all
clc

FName = 'Times New Roman';
FSize = 12;

C = zeros(M,1);
for i=1:M
    C(i,1) = std(ym(2*i,:));
end

clr = lines(3);
yl = [0.5 2; 0 3; 0 3];

for j=1:3
    figure('Position',[100 100 600 600])
    for i=1:m
        subplot(m,1,i)
        plot(t/fs,squeeze(IAhat(i,:,:,j)),':','Color',0.6*[1 1 1])
        hold on
        plot(t/fs,median(squeeze(IAhat(i,:,:,j)),2),'.','Color',clr(i,:))
        plot(t/fs,C(i)\IA(i,:),'--k','LineWidth',2)
        xlim([4 8]), ylim(yl(i,:))
        grid on
        xlabel('Time [s]')
        ylabel('IA (-)')
        set(gca,'FontName',FName,'FontSize',FSize)
    end
    
    set(gcf,'PaperPositionMode','auto')
    print(['Figures\Ex1_IAestimatesMC_',Methods{j},'_SNR',num2str(SNR)],'-dpng','-r300')
    print(['Figures\Ex1_IAestimatesMC_',Methods{j},'_SNR',num2str(SNR)],'-dmeta','-r300')
    
    figure('Position',[700 100 600 600])
    for i=1:m
%         subplot(m,1,i)
        plot(t/fs,squeeze(IFhat(i,:,:,j)),':','Color',0.6*[1 1 1])
        hold on
        plot(t/fs,median(squeeze(IFhat(i,:,:,j)),2),'.','Color',clr(i,:))
        plot(t/fs,IF(i,:),'--k','LineWidth',2)
    end
    xlim([4 8]), ylim([20 160])
    grid on
    xlabel('Time [s]')
    ylabel('IF (Hz)')
    set(gca,'FontName',FName,'FontSize',FSize)
    
    set(gcf,'PaperPositionMode','auto')
    print(['Figures\Ex1_IFestimatesMC_',Methods{j},'_SNR',num2str(SNR)],'-dpng','-r300')
    print(['Figures\Ex1_IFestimatesMC_',Methods{j},'_SNR',num2str(SNR)],'-dmeta','-r300')
    
end

%% Error analysis
close all
clc

T = 501:N;
rss_sss = zeros(m,Nmc,3);
rss_sssIA = zeros(m,Nmc,3);
rss_sssIF = zeros(m,Nmc,3);
for i=1:3
    errIA = repmat( diag(C)\IA, 1,1,Nmc ) - IAhat(:,:,:,i);
    errIF = repmat( IF, 1,1,Nmc ) - IFhat(:,:,:,i);
    err = y0 - yhat(:,:,:,i);
    
    rss_sssIA(:,:,i) = squeeze( sum( errIA(:,T,:).^2, 2)./sum( IA(:,T,:).^2, 2 ) );
    rss_sssIF(:,:,i) = squeeze( sum( errIF(:,T,:).^2, 2)./sum( IF(:,T,:).^2, 2 ) );
    rss_sss(:,:,i) = squeeze( sum( err(:,T,:).^2, 2)./sum( y0(:,T,:).^2, 2 ) );
end

figure('Position',[100 100 600 600])
subplot(311)
bar(log10(squeeze(median(rss_sss,2))'))
ylabel('1-R^2 - Signal')
set(gca,'XTickLabel',Methods)
ylim([-4 0])
grid on
set(gca,'FontName',FName,'FontSize',FSize)

subplot(312)
bar(log10(squeeze(median(rss_sssIA,2))'))
ylabel('1-R^2 - IA')
set(gca,'XTickLabel',Methods)
grid on
set(gca,'FontName',FName,'FontSize',FSize)
ylim([-4 0])

subplot(313)
bar(log10(squeeze(median(rss_sssIF,2))'))
ylabel('1-R^2 - IF')
ylim([-6 0])
set(gca,'XTickLabel',Methods)
grid on
set(gca,'FontName',FName,'FontSize',FSize)

set(gcf,'PaperPositionMode','auto')
print(['Figures\Ex1_ErrorPerformanceMC_SNR',num2str(SNR)],'-dpng','-r300')
print(['Figures\Ex1_ErrorPerformanceMC_SNR',num2str(SNR)],'-dmeta','-r300')

%% ------------------------------------------------------------------------
function [t,y,ym,IA,IF,Psi] = SignalGenerator(N,fs,SNR)

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

%-- Individual modes
ym = zeros(6,N);
Phi0 = 2*pi*rand(3,1);
for i=1:3
    ym(2*i-1,:) = IA(i,:).*cos( cumsum( 2*pi*IF(i,:)/fs ) + Phi0(i) );
    ym(2*i,:) = IA(i,:).*sin( cumsum( 2*pi*IF(i,:)/fs ) + Phi0(i) );
end

%-- Measured signals
y = Psi*ym;
sigmaN = sqrt( min(var(y,[],2))*10.^( -SNR/10 ) );
y = y + sigmaN*randn(d,N);

end