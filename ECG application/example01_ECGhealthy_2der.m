clear
close all
clc

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
[N,n] = size(y);

%% Pt.2 : Estimating the signal components with the diagonal SS method
close all
clc

% Initialization
M = 34;
Orders = (1:M);

Optimize = true;

if Optimize
    IniGuess.Variances = [1e-4 1e-6];
    IniGuess.TargetFrequency = 2*pi*1.07/fs;
    Niter = 100;
    T = 1201:2200;
    
    [Modal,logMarginal,HyperPar,Initial] = MO_DSS_JointEKF_MultiHar_Integrated_EM(y(T,:)',Orders,Niter,IniGuess);
    
    save('Optimized_HyperPar_2der','HyperPar','Initial','IniGuess','logMarginal','Niter')
    
else
    load('Optimized_HyperPar_2der','HyperPar','IniGuess','logMarginal','Niter')
end

[Modal,logMarginal] = MO_DSS_JointEKF_MultiHar_Integrated(y',Orders,HyperPar,Initial);

%% Pt.3 : Showing results
close all
clc

delta = 8;
lead = 1;

xl = [10 20];
FName = 'Times New Roman';
FSize = 12;

clr = lines(4);

figure('Position',[100 100 900 900])
plot(t,y(:,lead),'Color',clr(1,:),'LineWidth',1.5)   
text(xl(1)-0.01*diff(xl),0,lead_names{lead},...
    'FontName',FName,'FontSize',FSize,'HorizontalAlignment','right')

hold on
for m=1:min(10,M)
    plot(t,Modal.ym(2*m,:)+delta*m,'Color',clr(2,:),'LineWidth',1.5)   
    plot(t,Modal.ym(2*m-1,:)+delta*m,'Color',clr(2,:),'LineWidth',1.5)   
    
    text(xl(1)-0.01*diff(xl),delta*m,['m = ',num2str(m)],...
        'FontName',FName,'FontSize',FSize,'HorizontalAlignment','right')
end

xlim(xl)
grid on
if m==M, xlabel('Time [s]'), end
set(gca,'YTickLabel','')
set(gca,'FontName',FName,'FontSize',FSize)
xlabel('Time [s]')

% set(gcf,'PaperPositionMode','auto')
% print('Figures\ECGcontrol_Modes','-dpng','-r300')

%% Performance after EM optimization
close all
clc

yhat = HyperPar.Psi*Modal.ym;
err = y' - yhat;

%- Reconstruction error plot
figure('Position',[100 100 900 300])
bar(100*var(err,[],2)./var(y)')
set(gca,'XTickLabel',lead_names)
set(gca,'FontName',FName,'FontSize',FSize)
ylabel({'Reconstruction error';'RSS/SSS (%)'})

% set(gcf,'PaperPositionMode','auto')
% print('Figures\ECGcontrol_ErrorPerf','-dpng','-r300')

%- Optimized covariances plot
figure('Position',[100 400 900 360])
subplot(121)
bar(diag(HyperPar.Q(1:2:2*M,1:2:2*M)))
xlabel('Harmonic index')
ylabel('State innov. variance','Interpreter','latex')
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

% set(gcf,'PaperPositionMode','auto')
% print('Figures\ECGcontrol_Covariances','-dpng','-r300')

%- Optimized mixing matrix plot
figure('Position',[1000 100 900 360])
imagesc(HyperPar.Psi)
set(gca,'YTickLabel',lead_names,'YTick',1:n)
xlabel('Harmonic index')
set(gca,'FontName',FName,'FontSize',FSize)
set(gca,'CLim',max(max(abs(HyperPar.Psi)))*[-1 1])
cbar = colorbar;
cbar.Label.String = 'Amplitude of mixing matrix components';

% set(gcf,'PaperPositionMode','auto')
% print('Figures\ECGcontrol_MixMatrix','-dpng','-r300')

%%
close all
clc

leads = [1 2];
M_max = min(M,4);
yl = [-1 1; -1 1; -1 1; -1 1];

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

% set(gcf,'PaperPositionMode','auto')
% print('Figures\ECGcontrol_HarmonicComponents','-dpng','-r300')

%%
close all
clc

r = 1;

figure('Position',[100 100 900 800])
plot3(zeros(1,N),t,y(:,r),'LineWidth',1.5,'Color',clr(2,:))
hold on
for i=1:min(12,M)
    
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

xl = [20 30];

psi = HyperPar.Psi(r,ind);
psi = sqrt( psi(1:2:end).^2 + psi(2:2:end).^2 );
A = diag(psi)*Modal.Am;
A = A(1:M-1,:);
Mmax = 28;

figure('Position',[100 100 900 800])
subplot(4,1,4)
plot(t,y(:,r))
xlim(xl)
grid on
ylabel('Normalized amplitude')
xlabel('Time [s]')
set(gca,'FontName',FName,'FontSize',FSize)

subplot(4,1,1:3)
imagesc(t,1:Mmax,A(1:Mmax,:))
axis xy
colorbar('Location','northoutside')
xlim(xl)
grid on
ylabel('Harmonics')
set(gca,'FontName',FName,'FontSize',FSize)
set(gca,'CLim',[0 0.6])

%%
close all
clc

figure
subplot(211)
plot(t,Modal.ym(end-1:end,:))
xlim(xl)

subplot(212)
plot(t,err)
xlim(xl)

figure
for i=1:n
    subplot(n,1,i)
    plot(t,y(:,i),t,yhat(i,:))
    xlim(xl)
end

Nf = 4096;

figure
for i=1:n
    subplot(n,1,i)
    pwelch([y(:,i) err(i,:)'],hanning(Nf),3*Nf/4,Nf,fs)
%     xlim(xl)
end

