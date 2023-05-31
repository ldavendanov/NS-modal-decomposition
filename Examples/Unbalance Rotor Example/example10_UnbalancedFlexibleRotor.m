clear
close all
clc

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

Y = fft(y);
frec = fs*(0:size(Y,1)-1)/size(Y,1);

figure
subplot(211)
plot(time,y)

subplot(212)
semilogy(frec,abs(Y))

%% Decomposition using EMD
close all
clc

y_emd = [];

figure('Position',[100 100 900 360])
for i=1:3
    
    figure(1)
    [imf,res] = emd(y(:,i));

    y_emd = [y_emd, imf];

    subplot(1,3,i)
    hht(imf,fs,'FrequencyLimits',[0 10])
    grid on

end

A_imf = abs(hilbert(y_emd));
phi_imf = unwrap(angle(hilbert(y_emd)));
f_imf = diff([zeros(1,size(y_emd,2)); phi_imf])*fs/(2*pi);

figure
plot(time,f_imf)
hold on
plot(time,Omega/(2*pi))

%% Multivariate SSA
close all
clc

na = 8;
H = zeros(N-na,na*n);
tau = 1:N-na;
for i=1:na
    ind = (1:n)+n*(i-1);
    H(:,ind) = y(tau+i-1,:);
end

[ypca,S,V] = svd(H,'econ');
lambda = diag(S).^2;
ypca = ypca./std(ypca);

npc = find( cumsum(lambda)/sum(lambda) > 1-1e-4, 1, "first" );

figure
plot(cumsum(lambda)/sum(lambda))

figure
for i=1:npc
    plot(time(tau),ypca(:,i)+5*i)
    hold on
end

yh = hilbert(ypca);
Ah = abs(yh);
fh = diff([zeros(1,n*na); unwrap(angle(yh))])*fs/(2*pi);

mnfh = mean(fh);


%% SSdiagonal decomposition 
close all
clc

clear Modal

addpath('..\..\Core\')

M = 4;
Orders = 1:M;
Niter = 100;

%-- KF harmonic decomposition 
tic
InitialGuess.TargetFrequency = OmegaRef/fs;
InitialGuess.Variances = [1e-2 1e-2 1e-6];
[Modal{1},logMarginal{1},HyperPar{1},Initial{1}] = MO_DSS_JointEKF_MultiHar_EM(y',Orders,Niter,InitialGuess);
CT(1) = toc;

print('Figures\OptimizationKFHarm','-r300','-dpng')

%-- Block diagonal SS
tic
InitialGuess.TargetFrequency = OmegaRef/fs*(1:M);
InitialGuess.Variances = [1e-2 1e-2 1e-6];
[Modal{2},logMarginal{2},HyperPar{2},Initial{2}] = MO_DSS_JointEKF_EM(y',M,Niter,InitialGuess);
CT(2) = toc;

print('Figures\OptimizationKFdiag','-r300','-dpng')

%% Results
close all
clc

TLim = [0 20];

clr = lines(4);
LStyle = {'-','--','-.'};
cmap = OriBlue;
p = zeros(2,1);

figure('Position',[100 100 450 480])
for j=1:2
    for i=1:M
        p(j) = plot(time,Modal{j}.ym(2*i-1,:)+5*i,LStyle{j},'Color',clr(j,:),'LineWidth',1.5);
        hold on
    end
    grid on
end

for i=1:M
    p(3) = plot(time(tau),-ypca(:,2*i)+5*i,LStyle{3},'Color',clr(3,:),'LineWidth',1.5);
    hold on
end
grid on

xlim(TLim)
xlabel('Time (s)')
yticks([])
legend(p,{'KF-Harmonic','KF-Diagonal'},'Location','NorthOutside','Orientation','Horizontal')

for i=1:M
    text(0.98*TLim(1),5*i,['Mode ',num2str(i)],'HorizontalAlignment','right')
end

figure('Position',[600 100 450 480])
for i=1:M
    p(1) = plot(time,i*Modal{1}.omega*fs/(2*pi),'Color',clr(1,:),'LineWidth',1.5);
    hold on
end
for i=1:M
    p(2) = plot(time,Modal{2}.omega(i,:)*fs/(2*pi),'--','Color',clr(2,:),'LineWidth',1.5);
end
plot(time,Omega/(2*pi),':k')
grid on
xlim(TLim)
xlabel('Time (s)')
ylabel('Frequency (Hz)')

for i=1:M
    text(0.99*TLim(2),Modal{2}.omega(i,1)*fs/(2*pi),['Mode ',num2str(i)],'HorizontalAlignment','right','VerticalAlignment','top')
end
legend(p,{'KF-Harmonic','KF-Diagonal'},'Location','NorthOutside','Orientation','Horizontal')


figure('Position',[100 520 900 300])
for i=1:2
    subplot(2,1,i)
    imagesc(HyperPar{i}.Psi)
end

%% Results - Modal Components in time
close all
clc

TLim = [8 16];
YLim = [2 23];

clr = lines(4);
LStyle = {'-','-','-'};
cmap = OriBlue;

Method = {'KF-Harmonic','KF-Diagonal','SSA'};

figure('Position',[100 100 800 540])
for j=1:2
    subplot(1,3,j)
    for i=1:M
        p(j) = plot(time,Modal{j}.ym(2*i-1,:)+5*i,LStyle{j},'Color',clr(j,:),'LineWidth',1.5);
        hold on
    end
    grid on
    xlim(TLim)
    xlabel('Time (s)')
    yticks([])

    for i=1:M
        text(0.98*TLim(1),5*i,['Mode ',num2str(i)],'HorizontalAlignment','right')
    end
    title(Method{j})
    ylim(YLim)
end

subplot(133)
for i=1:M
    p(3) = plot(time(tau),-ypca(:,2*i)+5*i,LStyle{3},'Color',clr(3,:),'LineWidth',1.5);
    hold on
end
grid on

xlim(TLim)
xlabel('Time (s)')
yticks([])

for i=1:M
    text(0.98*TLim(1),5*i,['Mode ',num2str(i)],'HorizontalAlignment','right')
end
title(Method{3})
ylim(YLim)

set(gcf,'PaperPositionMode','auto')
print('Figures\UR_AnalysisResultsTime','-r300','-dpng')

%% Results - Instantaneous frequencies
close all
clc

figure('Position',[100 100 300 540])

for i=1:3
    p(3) = plot(time(tau),fh(:,2*i-1),':','Color',clr(3,:),'LineWidth',1.5);
    hold on
end
for i=1:M
    p(1) = plot(time,i*Modal{1}.omega*fs/(2*pi),'Color',clr(1,:),'LineWidth',1.5);
end
for i=1:M
    p(2) = plot(time,Modal{2}.omega(i,:)*fs/(2*pi),'-','Color',clr(2,:),'LineWidth',1.5);
end
    
xlim(TLim)
ylim([0 8])
grid on
xlabel('Time (s)')
ylabel('Frequency (Hz)')

legend(p,Method)

set(gcf,'PaperPositionMode','auto')
print('Figures\UR_AnalysisResultsIF','-r300','-dpng')

% plot(time,Omega/(2*pi),':k')
% grid on
% xlim(TLim)
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% 
% for i=1:M
%     text(0.99*TLim(2),Modal{2}.omega(i,1)*fs/(2*pi),['Mode ',num2str(i)],'HorizontalAlignment','right','VerticalAlignment','top')
% end
% legend(p,{'KF-Harmonic','KF-Diagonal'},'Location','NorthOutside','Orientation','Horizontal')



%% Reconstruction error
close all
clc

yhat{1} = HyperPar{1}.Psi * Modal{1}.ym;
yhat{2} = HyperPar{2}.Psi * Modal{2}.ym;
indT = 20:N;

RSS_SSS = var( y(indT,:) - yhat{1}(:,indT)' )./var(y(indT,:));
RSS_SSS(2,:) = var( y(indT,:) - yhat{2}(:,indT)' )./var(y(indT,:));

e{1} = (y-yhat{1}')./std(y);
e{2} = (y-yhat{2}')./std(y);
E{1} = fft(e{1})/size(e{1},1);
E{2} = fft(e{2})/size(e{2},1);

y_title = {{'Normalized error';'Horizontal (m/s^2)'},...
           {'Normalized error';'Vertical (m/s^2)'},...
           {'Normalized error';'Axial (m/s^2)'}};

figure('Position',[100 100 900 450])
for i=1:n
    subplot(3,2,2*i-1)
    plot(time,e{1}(:,i))
    hold on
    plot(time,e{2}(:,i))
    xlim(TLim)
    grid on
    ylabel(y_title{i})
    xlabel('Time (s)')

    subplot(3,2,2*i)
    semilogy(frec,abs(E{1}(:,i)))
    hold on
    semilogy(frec,abs(E{2}(:,i)))
    xlim([0 10]), 
    ylim([1e-5 1e0])
    grid on
    yticks([1e-4 1e-2 1e0])
    ylabel(y_title{i})
    xlabel('Frequency (Hz)')
end

set(gcf,'PaperPositionMode','auto')
print('Figures\UR_ErrorAnalysis','-r300','-dpng')

figure('Position',[100 560 900 240])
bar(log10(RSS_SSS'))
ylabel('$\log_{10}$ RSS/SSS','Interpreter','latex')
xticklabels({'Horizontal','Vertical','Axial'})
legend({'KF-Harmonic','KF-Diagonal'},'Location','northeastoutside')

set(gcf,'PaperPositionMode','auto')
print('Figures\UR_SummaryRSS','-r300','-dpng')




