clear
close all
clc

addpath('..\Core\')

%% Pt.1 : Creating the simulation signal
close all
clc

%-- Signal features
N = 4e3;                                                                    % Number of samples
fs = 5e2;                                                                   % Sampling frequency
t = 0:N-1;                                                                  % Time vector
SNR = 40;                                                                   % Signal to noise ratio

Aii = [1 0.2; 0.8 0.4; -0.2 1.0];                                           % Amplitude coefficients
fii = [50 -20 10];                                                          % Frequency coefficients
alpha = 0.5;                                                                % Cyclic frequency
Psi =  [1.0 0.8 0.5 0.5 0.2 0.8;
        0.5 0.5 0.8 1.0 0.4 0.6;
        0.4 0.9 0.1 0.2 1.0 0.0];
d = size(Psi,1);

for i=1:d
    Psi(i,:) = Psi(i,:)/norm(Psi(i,:));
end

%-- Instantaneous Amplitudes
IA = zeros(3,N);
IA(1,:) = Aii(1,1) + Aii(1,2)*cos(2*pi*alpha*t/fs);
IA(2,:) = Aii(2,1) + Aii(2,2)*sin(2*pi*alpha*t/fs);
IA(3,:) = max( Aii(3,1) + Aii(3,2)*abs( cos(2*pi*alpha*t/fs) ), 0 );

%-- Instantaneous frequencies
IF = fii(1,1) + fii(1,2)*sin( 2*pi*alpha(1)*t/fs ) + fii(1,3)*cos( 2*pi*alpha(1)*t/fs );

%-- Individual modes
ym = zeros(6,N);
Phi0 = 2*pi*rand(3,1);
for i=1:3
    ym(2*i-1,:) = IA(i,:).*cos( cumsum( i*2*pi*IF/fs ) + Phi0(i) );
    ym(2*i,:) = IA(i,:).*sin( cumsum( i*2*pi*IF/fs ) + Phi0(i) );
end

%-- Measured signals
y = Psi*ym;
sigmaN = sqrt( min(var(y,[],2))*10.^( -SNR/10 ) );
y = y + sigmaN*randn(d,N);

xl = [0 4];
FName = 'Times New Roman';
FSize = 12;

figure('Position',[100 100 600 600])
for i=1:d
    subplot(d,1,i)
    plot(t/fs,y(i,:))
    xlim(xl)
    grid on
    set(gca,'FontName',FName,'FontSize',FSize)
    ylabel(['$y_{',num2str(i),',t}$'],'Interpreter','latex')
    if i==3, xlabel('Time [s]'), end
end

set(gcf,'PaperPositionMode','auto')
print('Figures\Ex1_SignalSample','-dpng','-r300')

figure('Position',[100 100 600 600])
subplot(211)
plot(t/fs,IA)
ylim([0 2]), xlim(xl)
set(gca,'FontName',FName,'FontSize',FSize)
grid on
ylabel('IA')
xlabel('Time [s]')
legend({'Mode 1','Mode 2','Mode 3'},'Orientation','horizontal')

subplot(212)
plot(t/fs,IF)
xlim(xl), ylim([0 200])
grid on
set(gca,'FontName',FName,'FontSize',FSize)
ylabel('IF [Hz]')
xlabel('Time [s]')


%% Pt.2 : Estimating the signal components with the diagonal SS method
close all
clc

% Initialization
M = 3;
IniGuess.Variances = [1e-4 1e-5];
% IniGuess.TargetFrequencies = 2*pi*IF(:,1)/fs;

% Estimate without EM optimization
[InitialValues,HPar] = VAR_initialization(y(:,1:200*M),M,IniGuess);
[Modal{1},logMarginal{1}] = MO_DSS_JointEKF(y,M,'KS',InitialValues,HPar);

% Estimate after EM optimization
Niter = 50;
Orders = [1 2 3];
IniGuess.TargetFrequency = 2*pi*IF(1)/fs;
IniGuess.Variances = [1e-4 1e-5];
[Modal{2},logMarginal{2},HyperPar] = MO_DSS_JointEKF_MultiHar_EM(y,Orders,Niter,IniGuess);

% save('OptimizedHyperParams','HyperPar')

%% Pt.3 : Showing results
close all
clc

xl = [0.0 4.0];

C = zeros(M,1);
for i=1:M
    C(i,1) = std(ym(2*i,:));
    C(i,2) = std(Modal{1}.ym(2*i,:));
end

clr = lines(3);

figure('Position',[100 100 600 500])
for m=1:M
    subplot(M,1,m)
    plot(t/fs,IA(m,:)/C(m),'--k')
    hold on
    plot(t/fs,Modal{1}.Am(m,:),'Color',clr(1,:),'LineWidth',1.5)
    plot(t/fs,Modal{2}.Am(m,:),'Color',clr(2,:),'LineWidth',1.5)
    xlim(xl)
    ylim([0 3])
    grid on
    if m==3, xlabel('Time [s]'), end
    ylabel(['IA Mode ',num2str(i)])
    set(gca,'FontName',FName,'FontSize',FSize)
    
end

set(gcf,'PaperPositionMode','auto')
print('Figures\Ex1_IAestimates','-dpng','-r300')

figure('Position',[700 100 600 500])
for m=1:M
    plot(t/fs,m*IF,'--k')
    hold on
    plot(t/fs,Modal{1}.omega(m,:)*fs/(2*pi),'Color',clr(1,:),'LineWidth',1.5)
    plot(t/fs,m*Modal{2}.omega*fs/(2*pi),'Color',clr(2,:),'LineWidth',1.5)
    xlim(xl)
    grid on
    xlabel('Time [s]')
    ylabel('IF [Hz]')
    set(gca,'FontName',FName,'FontSize',FSize)
    text(0.1,m*IF(1),['Mode ',num2str(m)],'FontName',FName,'FontSize',FSize-1)
end
legend({'Original','Manual adjustment','EM optimization'},...
    'Location','northoutside','Orientation','horizontal')

set(gcf,'PaperPositionMode','auto')
print('Figures\Ex1_IFestimates','-dpng','-r300')

%% Plotting results - Hyperparameter estimates
close all
clc

q = diag(HyperPar.Q);

Qth = HyperPar.Q(2*M+1:end,2*M+1:end);

figure('Position',[100 100 600 600])
subplot(311)
plot(log10(diag(HPar.R)),'--','LineWidth',2)
hold on
plot(log10(diag(HyperPar.R)),'LineWidth',2)
xlabel('Ouput index')
ylabel('$\log_{10}\sigma_{\varepsilon_{i}}^2$','Interpreter','latex')
set(gca,'FontName',FName,'FontSize',FSize,'XTick',1:3)
ylim([-5 -0.5]), xlim([0.5 3.5])
legend({'Initial','After optimization'},'Orientation','horizontal')
grid on

subplot(312)
plot(1:2*M,log10( IniGuess.Variances(1) )*ones(1,2*M),'--','LineWidth',2)
hold on
plot(log10(q(1:2*M)),'LineWidth',2)
xlabel('State index')
ylabel('$\log_{10} \sigma_{u_{i}}^2$','Interpreter','latex')
set(gca,'FontName',FName,'FontSize',FSize)
ylim([-6.2 -2.2]), xlim([0.5 2*M+0.5])
legend({'Initial','After optimization'},'Orientation','horizontal')
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
for i=1:d
    
    psi1 = (Psi(i,1:2:end).*C(:,1)').^2 + (Psi(i,2:2:end).*C(:,1)').^2;
    psi2 = HyperPar.Psi(i,1:2:end).^2 + HyperPar.Psi(i,2:2:end).^2;
    psi3 = (HPar.Psi(i,1:2:end).*C(:,2)').^2 + (HPar.Psi(i,2:2:end).*C(:,2)').^2;
    subplot(d,1,i)
    plot( psi1, '--k' )
    hold on
    plot( psi3, 'Color', clr(1,:) )
    plot( psi2, 'Color', clr(2,:) )
    if i==3, xlabel('Output index'), end
    ylabel(['|\psi_',num2str(i),'|'])
    xlim([0.75 3.25])
    ylim([0 0.5])
    set(gca,'XTick',1:3)
    set(gca,'FontName','Times New Roman','FontSize',12)

    legend({'Original','VAR initialization','EM optimization'},'Orientation','horizontal')
    
end

set(gcf,'PaperPositionMode','auto')
print('Figures\Ex1_MixingMatrix','-dpng','-r300')