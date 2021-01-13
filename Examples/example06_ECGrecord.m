clear
close all
clc

addpath('..\Core\')

%% Pt.1 : Loading ECG recording

Fs = 1e3;
fs = 200;
load('Data\ECG PTB db\patient284\s0543_rem.mat','val')

y = resample(val',fs,Fs);
y = y(1:5e3,:);
y = detrend(y);
y = y/diag(std(y));
[N,n] = size(y);
t = (0:N-1)/fs;

pwelch(y,hamming(2^11),2^10,2^11,fs)

%% Pt.2 : Estimating the signal components with the diagonal SS method
close all
clc

% Initialization
M = 18;
Orders = 1:M;
IniGuess.Variances = [1e-3 1e-15];
IniGuess.TargetFrequency = 2*pi*1.27/fs;
Niter = 40;
[Modal,logMarginal,HyperPar] = MO_DSS_JointEKF_MultiHar_EM(y',Orders,Niter,IniGuess);

%% Pt.3 : Showing results
close all
clc

xl = [10.0 20.0];

clr = lines(3);
FName = 'Times New Roman';
FSize = 12;

figure('Position',[100 100 600 500])
for m=1:M
    subplot(M/2,2,m)
    plot(t,Modal.ym(2*m,:),'Color',clr(2,:),'LineWidth',1.5)
    xlim(xl)
    grid on
    if m==M, xlabel('Time [s]'), end
    ylabel(['Mode ',num2str(m)])
    set(gca,'FontName',FName,'FontSize',FSize)
    
end

figure('Position',[700 100 600 500])
plot(t,Modal.omega*fs/(2*pi),'Color',clr(1,:),'LineWidth',1.5)
% xlim(xl)
grid on
xlabel('Time [s]')
ylabel('IF [Hz]')
set(gca,'FontName',FName,'FontSize',FSize)

%%
close all
clc

xl = [10 20];

figure
for i=1:n
    subplot(n/3,3,i)
    plot(sqrt( HyperPar.Psi(i,1:2:end).^2 + HyperPar.Psi(i,2:2:end).^2 ))
    ylim([0 inf])
end

figure
for i=1:8
    psi = HyperPar.Psi(i,:);
%     psi(39:end) = 0;
    
    subplot(8,1,i)
    plot(t,y(:,i))
    hold on
    plot(t,psi*Modal.ym)
    xlim(xl)
end

%%
close all
clc

xl = [12 14];
r = 4;

figure
for i=1:M
    
    ind = (1:2) + 2*(i-1);
    psi = zeros(1,2*M);
    psi(ind) = HyperPar.Psi(r,ind);
    
    subplot(M/2,2,i)
    plot(t,y(:,r))
    hold on
    plot(t,psi*Modal.ym)
    xlim(xl)
end
