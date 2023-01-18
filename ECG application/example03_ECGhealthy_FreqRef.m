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


%% Pt.1a : Plotting the signal
close all
clc

FName = 'Times New Roman';
FSize = 12;
xl = [10 20];
% xl = [1201 2200]/fs;
clr = lines(4);
delta = 8;

figure('Position',[100 100 900 400])
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
print('Figures\ECGcontrol','-dpng','-r300')

%% Pt.2 : Estimating the signal components with the diagonal SS method
close all
clc

% Initialization
M = 44;
Orders = (1:2:M);

Optimize = true;

if Optimize
    IniGuess.Variances = [1e-2 1e-10];
    IniGuess.TargetFrequency = mean(omega_ref);
    Niter = 100;
    [~,logMarginal,HyperPar,Initial] = MO_DSS_JointEKF_FreqRef_EM(y(1001:2000,:)',omega_ref(1001:2000),Orders,Niter,IniGuess);
    
    save('Optimized_HyperPar','HyperPar','Initial','IniGuess','logMarginal','Niter')
    
else
    load('Optimized_HyperPar','HyperPar','Initial','IniGuess','logMarginal','Niter')
end

[Modal,logMarginal] = MO_DSS_JointEKF_FreqRef(y',omega_ref,Orders,HyperPar,Initial);

% [Modal,logMarginal,HyperPar] = MO_DSS_JointEKF_EM(y',M,Niter,IniGuess);

%%
close all
clc

lead = 2;
delta = 3;
xl = [0 2];

z = zeros(N,M/2);
for i=1:M/2
    ind = (1:2) + 2*(i-1);
    z(:,i) = (HyperPar.Psi(lead,ind)*Modal.ym(ind,:))';
end

figure('Position',[100 100 300 700])
for i=1:6
    plot(t-12,z(:,i)+delta*(i-1),'Color',clr(1,:))
    hold on
    text(xl(1)+0.05,delta*(i-1)+1,['Cmp. ',num2str(i)],'FontName',FName,'FontSize',FSize-2)
    xlim(xl)
end
set(gca,'FontName',FName,'FontSize',FSize,'YTick',[])
ylim([-2 18])
xlabel('Time [s]')

set(gcf,'PaperPositionMode','auto')
print(['Figures\ECGhealthy_NSHmodes',lead_names{lead}],'-dpng','-r300')

%% Signal spectrogram
close all
clc

lead = 1;
xl = [10 20];

Nf = 1024;
g = gausswin(Nf,12);
[Syy,ff,tt] = spectrogram(y(:,lead),g,Nf-1,Nf,fs);
Syy = abs(Syy).^2 * ( var(g)/Nf );
t_idx = tt >= xl(1) & tt <= xl(2);

figure('Position',[100 100 900 320])
imagesc(tt(t_idx),ff,10*log10(Syy(:,t_idx)))
axis xy
grid on
ylim([0 60])
xlim(xl)
ylabel('Frequency [Hz]')
xlabel('Time [s]')
cbar = colorbar('Location','eastoutside');
cbar.Label.String = 'Spectrogram [dB]';
set(gca,'CLim',[-70 -10],'FontName',FName,'FontSize',FSize)

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGspectrogram','-dpng','-r300')

%% More on spectrogram
close all
clc

xl = [9.2 15];
t_idx = t >= xl(1) & t < xl(1)+Nf/fs;

figure
plot(t,y(:,1))
hold on
xlim(xl)
grid on

alpha = [4 6 8 10 12];
for i=1:numel(alpha)
    win = gausswin(Nf,alpha(i));
    plot(t(t_idx),win+i)
end

set(gca,'XTick',10:15)

%% Pt.3 : Showing results
close all
clc

for lead = 1:6
    
    figure('Position',[100 100 900 800])
    plot3(zeros(1,N),t,y(:,lead),'LineWidth',1.5,'Color',clr(2,:))
    hold on
    for i=1:min(10,M/2)
        
        ind = (1:2) + 2*(i-1);
        psi = zeros(1,M);
        psi(ind) = HyperPar.Psi(lead,ind);
        
        plot3(i*ones(1,N),t,psi*Modal.ym,'Color',clr(1,:),'LineWidth',1.5)
        
    end
    ylim(xl)
    zlim([-4 8])
    view([60 50])
    grid on
    xlabel('Harmonic index')
    ylabel('Time [s]')
    zlabel('Normalized amplitude')
    set(gca,'FontName',FName,'FontSize',FSize+2)
    
    set(gcf,'PaperPositionMode','auto')
    print(['Figures\ECGcontrol_Modes',lead_names{lead}],'-dpng','-r300')
    
end

%% Performance after EM optimization
close all
clc

yhat = HyperPar.Psi*Modal.ym;
err = y' - yhat;

figure
for i=1:6
    subplot(3,2,i)
    plot(t,y(:,i))
    hold on
    plot(t,yhat(i,:))
    xlim(xl)
end

%- Reconstruction error plot
figure('Position',[100 100 900 300])
bar(100*sum( err(:,1001:2000).^2, 2 )./sum(y(1001:2000,:).^2)')
set(gca,'XTickLabel',lead_names)
set(gca,'FontName',FName,'FontSize',FSize+2)
ylabel({'Reconstruction error';'RSS/SSS (%)'})

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGcontrol_ErrorPerf','-dpng','-r300')

%- Optimized covariances plot
figure('Position',[100 400 900 360])
subplot(121)
bar(diag(HyperPar.Q(1:2:M,1:2:M)))
xlabel('Harmonic index')
ylabel('State innov. variance ','Interpreter','latex')
set(gca,'FontName',FName,'FontSize',FSize+2)
grid on

subplot(122)
imagesc(HyperPar.R(1:n,1:n)*1e4)
axis square
cbar = colorbar('Location','northoutside');
% set(gca,'CLim',max(max(abs(HyperPar.R(1:n,1:n))))*[-1 1]*1e3)
set(gca,'XTickLabel',lead_names,'YTickLabel',lead_names,'XTick',1:n,'YTick',1:n)
set(gca,'FontName',FName,'FontSize',FSize+2)
cbar.Label.String = 'Measurement noise covariance $\times 10^{-4}$';
cbar.Label.Interpreter = 'latex';

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGcontrol_Covariances','-dpng','-r300')

%- Optimized mixing matrix plot
figure('Position',[1000 100 900 360])
imagesc(abs(HyperPar.Psi))
set(gca,'YTickLabel',lead_names,'YTick',1:n)
xlabel('Harmonic index')
set(gca,'FontName',FName,'FontSize',FSize+2)
% set(gca,'CLim',max(max(abs(HyperPar.Psi)))*[-1 1])
cbar = colorbar;
cbar.Label.String = 'Absolute amplitude';

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGcontrol_MixMatrix','-dpng','-r300')

%%
close all
clc

leads = [4 6];
M_max = min(M/2,4);
yl = [-1 1; -2 2; -1 1; -1 1];

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
        
        Psi = zeros(1,M);
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

delta = 8;

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

%%
close all
clc

xl = [10 20];
lead = 2;

psi = HyperPar.Psi(lead,ind);
psi = sqrt( psi(1:2:end).^2 + psi(2:2:end).^2 );
A = diag(psi)*Modal.Am;
Mmax = M/2;
t_idx = t>=xl(1) & t<=xl(2);

figure('Position',[100 100 900 720])
subplot(3,1,3)
plot(t,y(:,lead))
xlim(xl)
grid on
ylabel('Normalized amplitude')
xlabel('Time [s]')
set(gca,'FontName',FName,'FontSize',FSize)

subplot(3,1,1:2)
imagesc(t(t_idx),1:Mmax,A(1:Mmax,t_idx))
axis xy
cbar = colorbar('Location','northoutside');
cbar.Label.String = 'Relative Amplitude';
xlim(xl)
grid on
ylabel('Harmonic component index')
xlabel('Time [s]')
set(gca,'FontName',FName,'FontSize',FSize)

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGcontrol_HarmonicComponents','-dpng','-r300')

%%
close all
clc

m = 3;
Nf = 2048;
spectrogram(Modal.ym(2*m-1,:),gausswin(Nf,8),Nf-5,Nf,fs)
