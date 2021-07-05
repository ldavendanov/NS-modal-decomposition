clear
close all
clc

addpath('..\Core\')

%% Pt.1 : Loading ECG recording

fs = 100;
[val, Fs] = rdsamp('ptbdb/patient025/s0087lre');

% rdann('ptbdb/patient040/s0130lre','atr')

%%

y = resample(val,fs,Fs);
y = y(1:2e3,7:12);
y = detrend(y);
y = y/diag(std(y));
[N,n] = size(y);
t = (0:N-1)/fs;

%% Pt.1a : Plotting the signal
close all
clc

lead_names = {'V1','V2','V3','V4','V5','V6'};

FName = 'Times New Roman';
FSize = 12;
xl = [0 5]+6.5;
% xl = [0 20];
clr = lines(4);
delta = 8;

figure('Position',[100 100 900 900])
for i=1:n
    plot(t,y(:,i)+delta*(i-1),'LineWidth',1.5,'Color',clr(1,:))
    hold on
    text(xl(1)+0.01*diff(xl),delta*(i-1)-delta/4,lead_names{i},'FontName',FName,'FontSize',FSize)
end
xlim(xl)
grid on
set(gca,'YTickLabel',[])
set(gca,'FontName',FName,'FontSize',FSize)
xlabel('Time [s]')

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGpathology','-dpng','-r300')

%% Pt.2 : Estimating the signal components with the diagonal SS method
close all
clc

% Initialization
M = 46;
Orders = (1:M);
IniGuess.Variances = [1e-3 1e-8];
IniGuess.TargetFrequency = 2*pi*0.88/fs;
Niter = 100;
[Modal,logMarginal,HyperPar] = MO_DSS_JointEKF_MultiHar_Integrated_EM(y',Orders,Niter,IniGuess);

% [Modal,logMarginal,HyperPar] = MO_DSS_JointEKF_EM(y',M,Niter,IniGuess);

%% Pt.3 : Showing results
close all
clc

delta = 8;
lead = 4;

figure('Position',[100 100 900 900])
plot(t,y(:,lead),'Color',clr(1,:),'LineWidth',1.5)   
text(xl(1)-0.01*diff(xl),0,lead_names{lead},...
    'FontName',FName,'FontSize',FSize,'HorizontalAlignment','right')

hold on
for m=1:min(10,M)
    plot(t,Modal.ym(2*m,:)+delta*m,'Color',clr(2,:),'LineWidth',1.5)   
%     plot(t,Modal.ym(2*m-1,:)+delta*m,'Color',clr(2,:),'LineWidth',1.5)   
    
    text(xl(1)-0.01*diff(xl),delta*m,['m = ',num2str(m)],...
        'FontName',FName,'FontSize',FSize,'HorizontalAlignment','right')
end

xlim(xl)
grid on
if m==M, xlabel('Time [s]'), end
set(gca,'YTickLabel','')
set(gca,'FontName',FName,'FontSize',FSize)
xlabel('Time [s]')

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGpathologic_Modes','-dpng','-r300')

figure
plot(t,1./Modal.omega*fs/(2*pi))
xlim(xl)

%% Performance after EM optimization
close all
clc

err = y' - HyperPar.Psi*Modal.ym;

%- Reconstruction error plot
figure('Position',[100 100 900 300])
bar(100*var(err,[],2)./var(y)')
set(gca,'XTickLabel',lead_names)
set(gca,'FontName',FName,'FontSize',FSize)
ylabel({'Reconstruction error';'RSS/SSS (%)'})

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGpathologic_ErrorPerf','-dpng','-r300')

%- Optimized covariances plot
figure('Position',[100 400 900 360])
subplot(121)
bar(diag(HyperPar.Q(1:2:2*M,1:2:2*M))*1e2)
xlabel('Harmonic index')
ylabel('State innov. variance $\times 10^{-2}$','Interpreter','latex')
set(gca,'FontName',FName,'FontSize',FSize)
grid on

subplot(122)
imagesc(HyperPar.R*1e3)
axis square
cbar = colorbar('Location','northoutside');
set(gca,'CLim',max(max(abs(HyperPar.R)))*[-1 1]*1e3)
set(gca,'XTickLabel',lead_names,'YTickLabel',lead_names,'XTick',1:n,'YTick',1:n)
set(gca,'FontName',FName,'FontSize',FSize)
cbar.Label.String = 'Measurement noise covariance $\times 10^{-3}$';
cbar.Label.Interpreter = 'latex';

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGpathologic_Covariances','-dpng','-r300')

%- Optimized mixing matrix plot
figure('Position',[1000 100 900 360])
imagesc(HyperPar.Psi)
set(gca,'YTickLabel',lead_names,'YTick',1:n)
xlabel('Harmonic index')
set(gca,'FontName',FName,'FontSize',FSize)
set(gca,'CLim',max(max(abs(HyperPar.Psi)))*[-1 1])
cbar = colorbar;
cbar.Label.String = 'Amplitude of mixing matrix components';

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGpathologic_MixMatrix','-dpng','-r300')

%%
close all
clc

leads = [1 4];
M_max = min(M,4);
yl = [-1.5 1; -2 1.5; -1 1; -1 1];

figure('Position',[100 100 900 800])
for j=1:2
    
    subplot(M_max+1,2,j)
    plot(t,y(:,leads(j)),'Color',clr(1,:),'LineWidth',1.5)
    xlim(xl)
    set(gca,'XTickLabel',[])
    ylabel(lead_names{leads(j)})
    set(gca,'FontName',FName,'FontSize',FSize)
    grid on
    
    for i=1:M_max
        
        Psi = zeros(1,2*M);
        ind = (1:2)+2*(i-1);
        Psi(ind) = HyperPar.Psi(leads(j),ind);
        
        subplot(M_max+1,2,2*i-j+3)
        plot(t,Psi*Modal.ym,'Color',clr(1,:),'LineWidth',1.5)
        xlim(xl)
        ylim(yl(i,:))
        ylabel(['Comp. ',num2str(i)])
        if i<M_max
            set(gca,'XTickLabel',[])
        else
            xlabel('Time [s]')
        end
        set(gca,'FontName',FName,'FontSize',FSize)
        grid on
    end
end

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGpathologic_HarmonicComponents','-dpng','-r300')

% figure
% for i=1:n
%     subplot(n/2,2,i)
%     plot(sqrt( HyperPar.Psi(i,1:2:end).^2 + HyperPar.Psi(i,2:2:end).^2 ))
%     ylim([0 inf])
% end
% 
% figure
% for i=1:n
%     psi = HyperPar.Psi(i,:);
% %     psi(2*(m-1)+(1:2)) = 0;
%     
%     subplot(n/2,2,i)
%     plot(t,y(:,i))
%     hold on
%     plot(t,psi*Modal.ym)
%     xlim(xl)
% end

%%
close all
clc

r = 4;

figure('Position',[100 100 900 800])
plot3(zeros(1,N),t,y(:,r),'LineWidth',1.5,'Color',clr(2,:))
hold on
for i=1:min(10,M)
    
    ind = (1:2) + 2*(i-1);
    psi = zeros(1,2*M);
    psi(ind) = HyperPar.Psi(r,ind);
    
    plot3(i*ones(1,N),t,psi*Modal.ym,'Color',clr(1,:),'LineWidth',1.5)

end
ylim(xl)
view([60 50])
grid on
xlabel('Harmonics')
ylabel('Time [s]')
zlabel('Normalized amplitude')
set(gca,'FontName',FName,'FontSize',FSize)

%%
close all
clc

psi = HyperPar.Psi(r,ind);
psi = sqrt( psi(1:2:end).^2 + psi(2:2:end).^2 );
A = diag(psi)*Modal.Am;

figure('Position',[100 100 900 800])
subplot(4,1,4)
plot(t,y(:,r))
xlim(xl)
grid on
ylabel('Normalized amplitude')
xlabel('Time [s]')
set(gca,'FontName',FName,'FontSize',FSize)

subplot(4,1,1:3)
imagesc(t,1:M,A)
colorbar('Location','northoutside')
xlim(xl)
grid on
ylabel('Harmonics')
set(gca,'FontName',FName,'FontSize',FSize)

%%
close all
clc

figure
subplot(131)
imagesc(sqrt(HyperPar.Psi(:,1:2:end).^2 + HyperPar.Psi(:,2:2:end).^2))

subplot(132)
semilogy(diag(HyperPar.Q(1:2*M,1:2*M)))

subplot(133)
imagesc(HyperPar.R)

%%
close all
clc

imagesc(Modal.ym(:,100:end)*Modal.ym(:,100:end)')
axis square

%%
close all
clc

plot(Modal.ym(end,:))