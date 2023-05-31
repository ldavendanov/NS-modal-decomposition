clear
close all
clc

% Load simulation data
load("URsystem_response.mat")
y = y(10001:20000,:);
Omega = Omega(10001:20000);
mnOmega = mean(Omega);
time = t(1:10000);

y = resample(y,1,5);
Omega = resample(Omega-mnOmega,1,5)+mnOmega;
time = time(1:5:end);
fs = fs/5;

[N,n] = size(y);

%% Test at different noise levels
close all
clc

addpath('..\..\Core\')

M = 4;
Orders = 1:M;
Niter = 100;

SNR = 0:5:60;

varY = min(var(y));
varN = 10.^(-SNR/10)*varY;

logML = zeros(size(SNR,2),2);
rss_sss = zeros(3,size(SNR,2),2);

for i=1:numel(SNR)

    %-- Add noise
    yN = y + sqrt(varN(i))*randn(size(y));

    %-- KF harmonic decomposition
    InitialGuess.TargetFrequency = OmegaRef/fs;
    InitialGuess.Variances = [varN(i) 1e-2 1e-6];
    [Modal,logMarginal,HyperPar,~] = MO_DSS_JointEKF_MultiHar_EM(yN',Orders,Niter,InitialGuess);

    %-- Error analysis
    yhat = (HyperPar.Psi*Modal.ym)';
    err = y - yhat;
    rss_sss(:,i,1) = var(err)./var(y);
    logML(i,1) = logMarginal(end);

    %-- Block diagonal SS
    InitialGuess.TargetFrequency = OmegaRef/fs*(1:M);
    InitialGuess.Variances = [varN(i) 1e-2 1e-6];
    [Modal,logMarginal,HyperPar,~] = MO_DSS_JointEKF_EM(yN',M,Niter,InitialGuess);

    %-- Error analysis
    yhat = (HyperPar.Psi*Modal.ym)';
    err = y - yhat;
    rss_sss(:,i,2) = var(err)./var(y);
    logML(i,2) = logMarginal(end);

end

save('UR_NoisePerformanceWhite','logML','rss_sss','SNR')

%% Colored noise
close all
clc

% Filter design for colored noise
f0 = 2; sigmaF = 2;
f = linspace(0,fs/2);
Hmag = exp(-(f-f0).^2/(2*sigmaF^2));
b = fir2(80,f*2/fs,Hmag);

figure
plot(f,Hmag.^2)

%%
addpath('..\..\Core\')

M = 4;
Orders = 1:M;
Niter = 100;

SNR = 0:5:60;

varY = min(var(y));
varN = 10.^(-SNR/10)*varY;

logML = zeros(size(SNR,2),2);
rss_sss = zeros(3,size(SNR,2),2);

for i=1:numel(SNR)

    %-- Add noise
    eta = randn(size(y));
    eta = filter(b,1,eta);
    yN = y + sqrt(varN(i))*eta./var(eta);

    %-- KF harmonic decomposition
    InitialGuess.TargetFrequency = OmegaRef/fs;
    InitialGuess.Variances = [varN(i) 1e-2 1e-6];
    [Modal,logMarginal,HyperPar,~] = MO_DSS_JointEKF_MultiHar_EM(yN',Orders,Niter,InitialGuess);

    %-- Error analysis
    yhat = (HyperPar.Psi*Modal.ym)';
    err = y - yhat;
    rss_sss(:,i,1) = var(err)./var(y);
    logML(i,1) = logMarginal(end);

    %-- Block diagonal SS
    InitialGuess.TargetFrequency = OmegaRef/fs*(1:M);
    InitialGuess.Variances = [varN(i) 1e-2 1e-6];
    [Modal,logMarginal,HyperPar,~] = MO_DSS_JointEKF_EM(yN',M,Niter,InitialGuess);

    %-- Error analysis
    yhat = (HyperPar.Psi*Modal.ym)';
    err = y - yhat;
    rss_sss(:,i,2) = var(err)./var(y);
    logML(i,2) = logMarginal(end);

    close all

end

save('UR_NoisePerformanceColored','logML','rss_sss','SNR')

%% Results
close all
clc

load('UR_NoisePerformanceWhite','rss_sss','SNR')
clr = lines(3);

figure('Position',[100 100 900 360])
t = tiledlayout('flow','TileSpacing','compact');
nexttile
for i=1:3
    semilogy(SNR,rss_sss(i,:,1)','Color',clr(i,:),'LineWidth',2)
    hold on
    semilogy(SNR,rss_sss(i,:,2)','--','Color',clr(i,:),'LineWidth',2)
end
grid on
xlabel('SNR (dB)')
ylabel('RSS/SSS')
text(1,1e0,'White noise','HorizontalAlignment','left','VerticalAlignment','top')

load('UR_NoisePerformanceColored','logML','rss_sss','SNR')
clr = lines(3);

nexttile
for i=1:3
    p(i) = semilogy(SNR,rss_sss(i,:,1)','Color',clr(i,:),'LineWidth',2);
    hold on
    p(i+3) = semilogy(SNR,rss_sss(i,:,2)','--','Color',clr(i,:),'LineWidth',2);
end
grid on
xlabel('SNR (dB)')
ylabel('RSS/SSS')
text(1,1e0,'Coloured noise','HorizontalAlignment','left','VerticalAlignment','top')

lgd = legend(p,{'Horz. KF-Harm','Vert. KF-Harm','Axial KF-Harm',...
    'Horz. KF-Diag','Vert. KF-Diag','Axial KF-Diag'});
lgd.Layout.Tile = 'north';
lgd.Orientation = 'horizontal';

set(gcf,'PaperPositionMode','auto')
print('Figures\UR_NoisePerformance','-r300','-dpng')
