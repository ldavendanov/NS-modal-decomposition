clear
close all
clc

addpath('..\Core\')

%% Pt.1 : Loading signal from test rig

f_hss = 12;
[y,tacho,Fs] = LoadRigVibration('H',f_hss,20);
fs = Fs/8;
y = resample(y,fs,Fs);
N = size(y,1);
t = (0:N-1);

pwelch(y,hamming(2^12),2^11,2^12,fs)

%% Pt.2 : Estimating the signal components with the diagonal SS method
close all
clc

% Initialization
M = 24;
Orders = 1:M;
IniGuess.Variances = [1e-3 1e-10];
IniGuess.TargetFrequency = 2*pi*f_hss/fs;
Niter = 50;
[Modal,logMarginal,HyperPar] = MO_DSS_JointEKF_MultiHar_EM(y(1:4e3,:)',Orders,Niter,IniGuess);

% save('OptimizedHyperParams','HyperPar')

%% Pt.3 : Showing results
close all
clc

xl = [0.0 6.0];

clr = lines(3);
FName = 'Times New Roman';
FSize = 12;

figure('Position',[100 100 600 500])
for m=1:M
    subplot(M/4,4,m)
    plot(t(1:4e3)/fs,Modal.ym(2*m,:),'Color',clr(2,:),'LineWidth',1.5)
    xlim(xl)
    grid on
    if m==3, xlabel('Time [s]'), end
    ylabel(['IA Mode ',num2str(m)])
    set(gca,'FontName',FName,'FontSize',FSize)
    
end

figure('Position',[700 100 600 500])
plot(t(1:4e3)/fs,Modal.omega*fs/(2*pi),'Color',clr(1,:),'LineWidth',1.5)
xlim(xl)
grid on
xlabel('Time [s]')
ylabel('IF [Hz]')
set(gca,'FontName',FName,'FontSize',FSize)

% set(gcf,'PaperPositionMode','auto')
% print('Figures\Ex1_IFestimates','-dpng','-r300')

%% Plotting results - Hyperparameter estimates
close all
clc

q = diag(HyperPar.Q);

Qth = HyperPar.Q(2*M+1:end,2*M+1:end);

figure('Position',[100 100 600 600])
subplot(311)
plot(log10(diag(HyperPar.R)),'LineWidth',2)
xlabel('Ouput index')
ylabel('$\log_{10}\sigma_{\varepsilon_{i}}^2$','Interpreter','latex')
set(gca,'FontName',FName,'FontSize',FSize,'XTick',1:3)
% ylim([-5 -0.5]), xlim([0.5 3.5])
% legend({'Initial','After optimization'},'Orientation','horizontal')
grid on

subplot(312)
plot(1:2*M,log10( IniGuess.Variances(1) )*ones(1,2*M),'--','LineWidth',2)
hold on
plot(log10(q(1:2*M)),'LineWidth',2)
xlabel('State index')
ylabel('$\log_{10} \sigma_{u_{i}}^2$','Interpreter','latex')
set(gca,'FontName',FName,'FontSize',FSize)
% ylim([-6.2 -2.2]), xlim([0.5 2*M+0.5])
% legend({'Initial','After optimization'},'Orientation','horizontal')
grid on

subplot(313)
plot(1:2*M,log10( IniGuess.Variances(2) )*ones(1,2*M),'--','LineWidth',2)
hold on
plot(log10(q(2*M+1:end)),'LineWidth',2)
xlabel('Param. index')
ylabel('$\log_{10} \sigma_{v_{i}}^2$','Interpreter','latex')
set(gca,'FontName',FName,'FontSize',FSize)
legend({'Initial','After optimization'},'Orientation','horizontal')
ylim([-6.4 -4.4]), xlim([0.5 2*M+0.5])
grid on


%%
close all
clc

clr = lines(2);

figure('Position',[100 100 600 600])
for i=1:2
    
    psi = HyperPar.Psi(i,1:2:end).^2 + HyperPar.Psi(i,2:2:end).^2;
    subplot(2,1,i)
    semilogy( psi, 'Color', clr(2,:) )
    if i==3, xlabel('Output index'), end
    ylabel(['|\psi_',num2str(i),'|'])
    set(gca,'FontName','Times New Roman','FontSize',12)
    
end

%%
close all
clc

figure
for i=1:2
    subplot(2,1,i)
    
    plot(t(1:4e3)/fs,y(1:4e3,i)')
    hold on
    plot(t(1:4e3)/fs,HyperPar.Psi(i,:)*Modal.ym)
    
    xlim(xl)
    
end

figure
for i=1:2
    subplot(2,1,i)
    
    pwelch(y(1:4e3,i)',hamming(1024),512,1024,fs)
    hold on
    pwelch(y(1:4e3,i)'-HyperPar.Psi(i,:)*Modal.ym,hamming(1024),512,1024,fs)
    
    for k=1:M
        plot(k*f_hss*[1 1],[-80 0],'--k')
    end
    
end