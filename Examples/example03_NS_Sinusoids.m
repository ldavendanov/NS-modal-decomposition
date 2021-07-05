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
SNR = 5;                                                                   % Signal to noise ratio

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
print('Figures\Ex1_SignalSample','-dmeta','-r300')

figure('Position',[100 100 600 600])
for i=1:d
    subplot(d,1,i)
    plot(t/fs,ym(2*i,:))
    xlim(xl)
    grid on
    set(gca,'FontName',FName,'FontSize',FSize)
    ylabel(['$q_{',num2str(i),',t}$'],'Interpreter','latex')
    if i==3, xlabel('Time [s]'), end
end

set(gcf,'PaperPositionMode','auto')
print('Figures\Ex1_OrgModesSample','-dpng','-r300')
print('Figures\Ex1_OrgModesSample','-dmeta','-r300')

figure('Position',[100 100 600 600])
subplot(211)
plot(t/fs,IA,'LineWidth',2)
ylim([0 2]), xlim(xl)
set(gca,'FontName',FName,'FontSize',FSize)
grid on
ylabel('IA')
xlabel('Time [s]')
legend({'Mode 1','Mode 2','Mode 3'},'Orientation','horizontal')

subplot(212)
plot(t/fs,IF,'LineWidth',2)
xlim(xl), ylim([0 200])
grid on
set(gca,'FontName',FName,'FontSize',FSize)
ylabel('IF [Hz]')
xlabel('Time [s]')
legend({'Mode 1','Mode 2','Mode 3'},'Orientation','horizontal')
print('Figures\Ex1_IAIFSample','-dpng','-r300')
print('Figures\Ex1_IAIFSample','-dmeta','-r300')

%% Pt.2a : Estimating the signal components with the diagonal SS method
close all
clc

% Initialization
M = 3;
IniGuess.Variances = [1e-4 1e-5];
IniGuess.Freqs = 2*pi*( [50 80 120] + 4*randn(1,3) )/fs;
% IniGuess.TargetFrequencies = 2*pi*IF(:,1)/fs;

% Estimate without EM optimization
[InitialValues,HPar] = VAR_initialization(y(:,1:200*M),M,IniGuess);
% InitialValues.x0(2*M+1:2:end) = cos(IniGuess.Freqs);
% InitialValues.x0(2*M+2:2:end) = sin(IniGuess.Freqs);
[Modal{1},logMarginal{1}] = MO_DSS_JointEKF(y,M,'KS',InitialValues,HPar);

% Estimate after EM optimization
Niter = 60;
IniGuess.Variances = [1e-6 1e-6];
IniGuess.Freqs = 2*pi*( [50 80 120] + 4*randn(1,3) )/fs;
[Modal{2},logMarginal{2},HyperPar,IniVal] = MO_DSS_JointEKF_EM(y,M,Niter,IniGuess);

save(['OptimizedHyperParams_SNR',num2str(SNR)],'HyperPar','IniVal')

%% Pt.2b : Estimating the signal components with the diagonal SS method (univariate)
close all
clc

% Estimate after EM optimization
Niter = 60;
IniGuess.Variances = [1e-5 1e-5];
IniGuess.Freqs = 2*pi*( [50 80 120] + 4*randn(1,3) )/fs;
[Modal{3},logMarginal{3},HyperParSO,IniValSO] = MO_DSS_JointEKF_EM(y(1,:),M,Niter,IniGuess);

save(['OptimizedHyperParamsSO_SNR',num2str(SNR)],'HyperParSO','IniValSO')

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
    plot(t/fs,Modal{3}.Am(m,:),'Color',clr(3,:),'LineWidth',1.5)
    xlim(xl)
    ylim([0 3])
    grid on
    if m==3, xlabel('Time [s]'), end
    ylabel(['IA Mode ',num2str(m)])
    set(gca,'FontName',FName,'FontSize',FSize)
    
end

set(gcf,'PaperPositionMode','auto')
print('Figures\Ex1_IAestimates','-dpng','-r300')
print('Figures\Ex1_IAestimates','-dmeta','-r300')

figure('Position',[700 100 600 500])
for m=1:M
    plot(t/fs,IF(m,:),'--k')
    hold on
    plot(t/fs,Modal{1}.omega(m,:)*fs/(2*pi),'Color',clr(1,:),'LineWidth',1.5)
    plot(t/fs,Modal{2}.omega(m,:)*fs/(2*pi),'Color',clr(2,:),'LineWidth',1.5)
    plot(t/fs,Modal{3}.omega(m,:)*fs/(2*pi),'Color',clr(3,:),'LineWidth',1.5)
    xlim(xl)
    grid on
    xlabel('Time [s]')
    ylabel('IF [Hz]')
    set(gca,'FontName',FName,'FontSize',FSize)
    text(0.1,IF(m,1),['Mode ',num2str(m)],'FontName',FName,'FontSize',FSize-1)
end
legend({'Original','MO - Manual','MO - EM','SO - EM'},...
    'Location','northoutside','Orientation','horizontal')

set(gcf,'PaperPositionMode','auto')
print('Figures\Ex1_IFestimates','-dpng','-r300')
print('Figures\Ex1_IFestimates','-dmeta','-r300')

%% Plotting results - Hyperparameter estimates
close all
clc

q = diag(HyperPar.Q);

Qth = HyperPar.Q(2*M+1:end,2*M+1:end);
R = Qth;
for i=1:2*M
    c = sqrt(Qth(i,i))*sqrt( diag( Qth )' );
    R(i,:) = R(i,:)./c;
end

figure('Position',[100 100 600 600])
subplot(311)
plot(log10(diag(HPar.R)),'--','LineWidth',2)
hold on
plot(log10(diag(HyperPar.R)),'LineWidth',2)
xlabel('Ouput index')
ylabel('$\log_{10}\sigma_{\varepsilon_{i}}^2$','Interpreter','latex')
set(gca,'FontName',FName,'FontSize',FSize,'XTick',1:3)
ylim([-5 -0.5]), xlim([0.5 3.5])
legend({'Initial','After optimization'},'Orientation','horizontal','Location','best')
grid on

subplot(312)
plot(1:2*M,log10( IniGuess.Variances(1) )*ones(1,2*M),'--','LineWidth',2)
hold on
plot(log10(q(1:2*M)),'LineWidth',2)
xlabel('State index')
ylabel('$\log_{10} \sigma_{u_{i}}^2$','Interpreter','latex')
set(gca,'FontName',FName,'FontSize',FSize)
ylim([-6.2 -2.2]), xlim([0.5 2*M+0.5])
legend({'Initial','After optimization'},'Orientation','horizontal','Location','best')
grid on

subplot(313)
plot(1:2*M,log10( IniGuess.Variances(2) )*ones(1,2*M),'--','LineWidth',2)
hold on
plot(log10(q(2*M+1:end)),'LineWidth',2)
xlabel('Param. index')
ylabel('$\log_{10} \sigma_{v_{i}}^2$','Interpreter','latex')
set(gca,'FontName',FName,'FontSize',FSize)
legend({'Initial','After optimization'},'Orientation','horizontal','Location','best')
ylim([-6.4 -4.4]), xlim([0.5 2*M+0.5])
grid on

set(gcf,'PaperPositionMode','auto')
print('Figures\Ex1_CovEstimates1','-dpng','-r300')


figure('Position',[700 100 600 500])
imagesc(R)
axis square
cbar = colorbar;
set(gca,'CLim',[-1 1])
xlabel('Param. index')
ylabel('Param. index')
cbar.Label.String = 'Correlation index';
set(gca,'FontName',FName,'FontSize',FSize)

set(gcf,'PaperPositionMode','auto')
print('Figures\Ex1_CovEstimates2','-dpng','-r300')

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
    p(2) = plot( psi3, 'Color', clr(1,:), 'LineWidth', 2 );
    hold on    
    p(3) = plot( psi2, 'Color', clr(2,:), 'LineWidth', 2 );
    p(1) = plot( psi1, '--k', 'LineWidth', 1 );
    if i==3, xlabel('Output index'), end
    ylabel(['|\psi_',num2str(i),'|'])
    xlim([0.75 3.25])
    ylim([0 0.5])
    set(gca,'XTick',1:3)
    set(gca,'FontName','Times New Roman','FontSize',12)

    legend(p,{'Original','VAR initialization','EM optimization'},'Orientation','horizontal')
    
end

set(gcf,'PaperPositionMode','auto')
print('Figures\Ex1_MixingMatrix','-dpng','-r300')
print('Figures\Ex1_MixingMatrix','-dmeta','-r300')

%% Mode estimation error
close all
clc

T = 501:N;
rss_sss = zeros(m,3);
rss_sssIA = zeros(m,3);
rss_sssIF = zeros(m,3);
for i=1:3
    errIA = diag(C(:,1))\IA - Modal{i}.Am;
    errIF = IF - Modal{i}.omega*fs/(2*pi);
    if i==1
        err = Psi*ym - HPar.Psi*Modal{i}.ym;
    elseif i==2
        err = Psi*ym - HyperPar.Psi*Modal{i}.ym;
    elseif i==3
        err = Psi(1,:)*ym - HyperParSO.Psi*Modal{i}.ym;
    end
    rss_sssIA(:,i) = sum( errIA(:,T).^2, 2)./sum( IA(:,T).^2, 2 );
    rss_sssIF(:,i) = sum( errIF(:,T).^2, 2)./sum( IF(:,T).^2, 2 );
    rss_sss(:,i) = sum( err(:,T).^2, 2)./sum( Psi*ym(:,T).^2, 2 );
end

subplot(311)
bar(1-rss_sss')
ylabel('R^2 - Signal')
ylim([0.8 1])

subplot(312)
bar(1-rss_sssIA')
ylabel('R^2 - IA')
ylim([0.6 1])

subplot(313)
bar(1-rss_sssIF')
ylabel('R^2 - IF')
ylim([0.998 1])

%% Other results
close all
clc

na = 10;
Phi = zeros(N-na,na*3);
tau = na+1:N;
for i=1:na
    ind = (1:3) + 3*(i-1);
    Phi(:,ind) = y(:,tau-i+1)';
end

[U,S,V] = svd(Phi,'econ');
lambda = diag(S).^2;

figure
plot(cumsum(lambda)/sum(lambda))

R = zeros(2*m,size(U,2));
for i=1:2*m
    D1 = sqrt(ym(i,tau)*ym(i,tau)');
%     D2 = sqrt(diag(U'*U));
    R(i,:) = ym(i,tau)*U/D1;
end

figure
imagesc(abs(R))

figure
for i=1:m
    subplot(3,1,i)
    ind = (1:2) + 2*(i-1);
    [r,idx0] = max(abs(R(ind,:)),[],2);
    [~,idx1] = max(r);
    idx = idx0(idx1);
    plot(t(tau),sqrt(lambda(idx))*U(:,idx))
    hold on
    plot(t(tau),ym(2*i,tau))
end

Uh = hilbert(U);
figure
for i=1:m
    subplot(3,1,i)
    ind = (1:2) + 2*(i-1);
    [r,idx0] = max(abs(R(ind,:)),[],2);
    [~,idx1] = max(r);
    idx = idx0(idx1);
    plot(t(tau),abs(sqrt(lambda(idx))*Uh(:,idx)))
    hold on
    plot(t(tau),IA(i,tau))
end

figure
for i=1:m
%     subplot(3,1,i)
    ind = (1:2) + 2*(i-1);
    [r,idx0] = max(abs(R(ind,:)),[],2);
    [~,idx1] = max(r);
    idx = idx0(idx1);
    plot( t(tau(2:end)), abs(diff(unwrap(angle(Uh(:,idx))))*fs/(2*pi)) )
    hold on
    plot(t(tau),IF(i,tau))
end
