clear
close all
clc

addpath('..\Core\')

%% Pt.1 : Loading ECG recording

Fs = 1e3;
fs = 500;
load('Data\ECG PTB db\patient284\s0543_rem.mat','val')

y = resample(val',fs,Fs);
y = y(5000+(1:16e3),1:6);
y = detrend(y);
y = y/diag(std(y));
[N,n] = size(y);
t = (0:N-1)/fs;

%% Select a signal segment
close all
clc

[Pyy,ff] = pwelch(y,hann(2^13),2^12,2^13,fs);
[~,ind] = max(Pyy,[],'all');
f0 = ff(ind);

figure
plot(t,y)

figure
plot(ff,10*log10(Pyy))
hold on
yl = get(gca,'YLim');
for i=1:40
    plot(i*f0*[1 1],yl,'--k')
end


%% Pt.2 : Estimating the signal components with the diagonal SS method
close all
clc

Y = y(1:2000,:) - mean(y(1:2000,:));
time = t(1:2000);

% Initialization
M = 50;
Orders = 1:M;
IniGuess.Variances = [5e-5 1e-5].^2;
IniGuess.TargetFrequency = 2*pi*f0/fs;
Niter = 30;
[Modal,logMarginal,HyperPar] = MO_DSS_JointEKF_MultiHar_EM(Y',Orders,Niter,IniGuess);

%% Pt.3 : Showing results
close all
clc

xl = [10.0 20.0];

clr = lines(3);
FName = 'Times New Roman';
FSize = 12;

figure('Position',[100 100 600 500])
for m=1:10
    subplot(10/2,2,m)
    plot(time,Modal.ym(2*m,:),'Color',clr(2,:),'LineWidth',1.5)
%     xlim(xl)
    grid on
    if m==M, xlabel('Time [s]'), end
    ylabel(['Mode ',num2str(m)])
    set(gca,'FontName',FName,'FontSize',FSize)
    
end

figure('Position',[700 100 600 500])
plot(time,Modal.omega*fs/(2*pi),'Color',clr(1,:),'LineWidth',1.5)
% xlim(xl)
grid on
xlabel('Time [s]')
ylabel('IF [Hz]')
set(gca,'FontName',FName,'FontSize',FSize)

%%
close all
clc

for i=1:6
    subplot(3,2,i)
    hht((diag(HyperPar.Psi(i,:))*Modal.ym)',fs,'FrequencyLimits',[0 50])
end

%%
close all
clc

xl = [10 20];

figure
for i=1:n
    subplot(n/2,2,i)
    plot(sqrt( HyperPar.Psi(i,1:2:end).^2 + HyperPar.Psi(i,2:2:end).^2 ))
    ylim([0 inf])
end

figure
for i=1:n
    psi = HyperPar.Psi(i,:);
    psi(2*(m-1)+(1:2)) = 0;
    
    subplot(n,1,i)
    plot(time,Y(:,i))
    hold on
    plot(time,psi*Modal.ym)
%     xlim(xl)
end

%%
close all
clc

xl = [12 14];
r = 2;

figure
for i=1:M
    
    ind = (1:2) + 2*(i-1);
    psi = zeros(1,2*M);
    psi(ind) = HyperPar.Psi(r,ind);
    
    subplot(M/2,2,i)
    plot(time,Y(:,r))
    hold on
    plot(time,psi*Modal.ym)
%     xlim(xl)
end

%%
close all
clc

figure
subplot(131)
imagesc(sqrt(HyperPar.Psi(:,1:2:end).^2 + HyperPar.Psi(:,2:2:end).^2))

subplot(132)
semilogy(diag(HyperPar.Q))

subplot(133)
imagesc(HyperPar.R)
