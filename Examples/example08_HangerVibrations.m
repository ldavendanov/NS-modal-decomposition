clear
close all
clc

% Load signal
data = load("Data\transverse_and_longitudinal_data_2021-11-22_one_hour.mat");
y = [data.sensorL data.sensorT];
y = detrend(y);
y = y ./ std(y);

time = time2num(data.sensorTime);
fs0 = 1/( time(2)-time(1) );
fs = 8;

b = fir1(300,2*[0.01 fs/2]/fs0,'bandpass');
y = filtfilt(b,1,y);

% Nf = 2^13;
% subplot(211)
% pwelch(y,hann(Nf),3*Nf/4,Nf,fs0 )
% xlim([0 fs/2])
% 
y = resample(y,fs,fs0);
% 
% Nf = 2^12;
% subplot(212)
% pwelch(y,hann(Nf),3*Nf/4,Nf,fs )

N = size(y,1);
time = (0:N-1)/fs;

%%
close all
clc

Nf = 2^11;

figure
for i=1:2
    [Syy,ff,tt] = spectrogram(y(:,i),gausswin(Nf),Nf-20,Nf,fs);
    subplot(2,1,i)
    imagesc(tt,ff,20*log10(abs(Syy)))
    xlim([0 3600])
    cbar = colorbar;
    set(gca,'CLim',[10 24])
    axis xy
end

%% Pt.2 : Estimating the signal components with the diagonal SS method
close all
clc

% Initialization
M = 8;
Orders = 1:M;
IniGuess.Variances = [1e-2 1e-10];
IniGuess.TargetFrequency = 2*pi*[0.398 0.405 0.45]'/fs;
Niter = 20;
[~,~,HyperPar,Initial] = MO_DSS_JointEKF_MultiHar_EM(y(1:4e3,:)',Orders,Niter,IniGuess);
[Modal,logMarginal] = MO_DSS_JointEKF_MultiHar(y',Orders,HyperPar,Initial);

%% Pt.3 : Showing results
close all
clc

xl = [0.0 3600];

clr = lines(3);
FName = 'Times New Roman';
FSize = 12;

figure('Position',[100 100 600 500])
for m=1:M
    subplot(M/2,2,m)
    plot(time,Modal.ym(2*m,:),'Color',clr(1,:),'LineWidth',1.5)
    hold on
    plot(time,Modal.ym(2*m+2*M,:),'Color',clr(2,:),'LineWidth',1.5)
    plot(time,Modal.ym(2*m+4*M,:),'Color',clr(3,:),'LineWidth',1.5)
    xlim(xl)
    grid on
    if m==3, xlabel('Time [s]'), end
    ylabel(['IA Mode ',num2str(m)])
    set(gca,'FontName',FName,'FontSize',FSize)
    
end

figure('Position',[700 100 600 500])
plot(time,Modal.omega*fs/(2*pi),'Color',clr(1,:),'LineWidth',1.5)
xlim(xl)
grid on
xlabel('Time [s]')
ylabel('IF [Hz]')
set(gca,'FontName',FName,'FontSize',FSize)

% set(gcf,'PaperPositionMode','auto')
% print('Figures\Ex1_IFestimates','-dpng','-r300')

%%
close all
clc

figure
for i=1:2
    subplot(2,1,i)
    
    plot(time,y(:,i)')
    hold on
    plot(time,HyperPar.Psi(i,:)*Modal.ym)
    
    xlim(xl)
    
end

figure
for i=1:2
    subplot(2,1,i)
    
    pwelch(y(:,i)',hamming(1024),512,1024,fs)
    hold on
    pwelch(HyperPar.Psi(i,:)*Modal.ym,hamming(1024),512,1024,fs)
    
    
end

%%
close all
clc

figure
for i=1:3
    subplot(3,1,i)
    plot(time,(Modal.Am((1:M)+2*(i-1),:)))
end

%%
close all
clc

figure
imagesc(HyperPar.Psi(:,1:M)./max(abs(HyperPar.Psi(:,1:M))))