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

leads = 10:12;
[y,fs,t,lead_names,RRinterval] = PreprocECG(val,Fs,leads);
omega_ref = (2*pi/fs)./RRinterval;


%% Pt.2 : Model Order Selection
close all
clc

IniGuess.Variances = [1e-4 1e-10];
IniGuess.TargetFrequency = mean(omega_ref);
Niter = 100;
T = 1201:2200;

logMarginal = zeros(Niter,23);
RSS_SSS = zeros(23,3);

% Initialization
for m = 2:2:46
    
    Orders = (1:2:m);    
    [Modal,logMarginal(:,m/2),HyperPar] = MO_DSS_JointEKF_MultiHar_Integrated_EM(y(T,:)',Orders,Niter,IniGuess);    
    yhat = HyperPar.Psi * Modal.ym;
    err = y(T,:) - yhat';
    
    RSS_SSS(m/2,:) = sum( err.^2 ) ./ sum( y(T,:).^2 );
    
    close all
    
end

save('MOS_perf','RSS_SSS','logMarginal')

%% MOS Results
close all
clc

FName = 'Times New Roman';
FSize = 12;

load('MOS_perf','RSS_SSS','logMarginal')

figure('Position',[100 100 900 360])
subplot(121)
semilogy(2:2:46,100*RSS_SSS,'-o')
grid on
set(gca,'FontName',FName,'FontSize',FSize)
xlabel('Model order')
ylabel('RSS/SSS (%)')


subplot(122)
plot(2:2:46,logMarginal(end,:))
grid on
set(gca,'FontName',FName,'FontSize',FSize)
xlabel('Model order')
ylabel('Log marginal likelihood')

%% Pt.3 : Sensitivity to initial values
close all
clc

Niter = 100;

logMarginal = zeros(Niter,10);
RSS_SSS = zeros(10,3);
lambda = logspace(-4,1,10);
T = 1201:2200;

% Initialization
for i=1:10
    
    IniGuess.Variances = [lambda(i) 1e-12];
    IniGuess.TargetFrequency = 2*pi*1.07/fs;
    
    Orders = (1:46);    
    [Modal,logMarginal(:,i),HyperPar] = MO_DSS_JointEKF_MultiHar_Integrated_EM(y(T,:)',Orders,Niter,IniGuess);    
    yhat = HyperPar.Psi * Modal.ym;
    err = y(T,:) - yhat';
    
    RSS_SSS(i,:) = sum( err.^2 ) ./ sum( y(T,:).^2 );
    
    close all
    
end

save('IC_perf','RSS_SSS','logMarginal','lambda')

%%
close all
clc

figure('Position',[100 100 900 360])
subplot(121)
loglog(lambda,100*RSS_SSS)
grid on

subplot(122)
loglog(lambda,logMarginal(end,:))
grid on