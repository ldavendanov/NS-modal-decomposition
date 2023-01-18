clear
close all
clc

addpath('C:\Users\ldav\Documents\MATLAB\wfdb-app-toolbox-0-10-0\mcode')
addpath('..\Core\')

%% Pt.1 : Loading ECG recording

remote_load = false;

if remote_load
    [val, Fs] = rdsamp('ptbdb/patient104/s0306lre');
    save('HealthyECG','val','Fs')
else
    load('HealthyECG','val','Fs')
end

leads = 7:12;
[y,fs,t,lead_names,RRinterval] = PreprocECG(val,Fs,leads);
omega_ref = (2*pi/fs)./RRinterval;
[N,n] = size(y);

%% Pt.2 : Model Order Selection
close all
clc

IniGuess.Variances = [1e-2 1e-10];
IniGuess.TargetFrequency = mean(omega_ref);
Niter = 100;                                                                % Maximum number of iterations
m = 2:2:50;                                                                 % Candidate orders

% Initialize computation matrices
logMarginal = zeros(Niter,numel(m));
RSS_SSS = zeros(6,numel(m));

% Initialization
for i=1:numel(m)
    
    M = m(i);
    Orders = (1:2:M);
    
    T = 1201:2200;
    [Modal,logMarginal(:,i),HyperPar] = MO_DSS_JointEKF_FreqRef_EM(y(T,:)',omega_ref(T),Orders,Niter,IniGuess);
    
    yhat = HyperPar.Psi*Modal.ym;
    err = y(T,:) - yhat';
    
    RSS_SSS(:,i) = sum( err.^2 ) ./ sum( y(T,:).^2 );
    
    close all
    
end

save('MOS_perf_FreqRef_1em2.mat','m','logMarginal','RSS_SSS')

%%
close all
clc

FName = 'Times New Roman';
FSize = 12;

load('MOS_perf_FreqRef_1em2.mat','m','logMarginal','RSS_SSS')

figure('Position',[100 100 900 300])
subplot(121)
semilogy(m,median(RSS_SSS*100),'-o','LineWidth',1.5)
hold on
semilogy(m,RSS_SSS*100,'--','Color',0.7*[1 1 1])
xlabel('Model order')
ylabel('RSS/SSS (%)')
set(gca,'FontName',FName,'FontSize',FSize)
grid on
% ylim([5e-3 1e-1])

subplot(122)
plot(m,logMarginal(end,:))
xlabel('Model order')
ylabel('Log marginal likelihood')
set(gca,'FontName',FName,'FontSize',FSize)
grid on
% ylim([18.5 21])

set(gcf,'PaperPositionMode','auto')
print('Figures\MOSresults','-dpng','-r300')

%% Pt.3 : Sensitivity to initial conditions
close all
clc

Niter = 100;                                                                 % Maximum number of iterations

% Initialize computation matrices
lambda = logspace(-6,2,20);
logMarginal = zeros(Niter,numel(lambda));
RSS_SSS = zeros(6,numel(lambda));

% Initialization
for i=1:numel(lambda)
    
    IniGuess.Variances = [lambda(i) 1e-10];
    IniGuess.TargetFrequency = mean(omega_ref);

    M = 12;
    Orders = (1:2:M);
    T = 1201:2200;
    
    [Modal,logMarginal(:,i),HyperPar] = MO_DSS_JointEKF_FreqRef_EM(y(T,:)',omega_ref(T),Orders,Niter,IniGuess);
    
    yhat = HyperPar.Psi*Modal.ym;
    err = y(T,:) - yhat';
    
    RSS_SSS(:,i) = sum( err.^2 ) ./ sum( y(T,:).^2 );
    
    close all
    
end

save('IC_perf_FreqRef.mat','lambda','logMarginal','RSS_SSS')

%%
close all
clc

FName = 'Times New Roman';
FSize = 12;

load('IC_perf_FreqRef.mat','lambda','logMarginal','RSS_SSS')

figure('Position',[100 100 900 300])
subplot(121)
loglog(lambda,median(RSS_SSS*100),'-o','LineWidth',1.5)
hold on
loglog(lambda,RSS_SSS*100,'--','Color',0.7*[1 1 1])
xlabel('Initial $\sigma_u^2$','Interpreter','latex')
ylabel('RSS/SSS (%)')
set(gca,'FontName',FName,'FontSize',FSize)
grid on
% ylim([3e-3 2e-2])

subplot(122)
semilogx(lambda,logMarginal(end,:))
xlabel('Initial $\sigma_u^2$','Interpreter','latex')
ylabel('Log marginal likelihood')
set(gca,'FontName',FName,'FontSize',FSize)
grid on
% ylim([20 22])

set(gcf,'PaperPositionMode','auto')
print('Figures\IC_performance','-dpng','-r300')

