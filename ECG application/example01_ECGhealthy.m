clear
close all
clc

addpath('..\Core\')

%% Pt.1 : Loading ECG recording

addpath('C:\Users\ldav\Documents\MATLAB\wfdb-app-toolbox-0-10-0\mcode')
[val, Fs] = rdsamp('ptbdb/patient104/s0306lre');

%%

fs = 100;
y = resample(val,fs,Fs);
y = y(1:8e3,7:12);
y = detrend(y);
y = y/diag(std(y));
[N,n] = size(y);
t = (0:N-1)/fs;

%% RR interval estimation
close all
clc

[pks,idx] = findpeaks(y(:,4),'MinPeakDistance',0.75,'MinPeakProminence',4);
locs = t(idx);

RR0 = diff(locs);
RRinterval = resample( diff(locs), locs(2:end), fs, 1, 1, 'spline' );

figure
subplot(211)
plot(t,y(:,4))
hold on
plot(locs,pks,'o')
grid on

subplot(212)
plot( locs(2)+(1:numel(RRinterval))/fs, RRinterval )
hold on
plot(locs(2:end),RR0,'o')
grid on

%% Remove baseline interference
close all
clc

idx_bl = idx - round(1*diff([0 idx'])/4)';
locs_bl = t(idx_bl);

f = zeros(size(y));
for i=1:n
    [baseline,gof,output] = fit( locs_bl', y(idx_bl,i), 'cubicinterp' );
    f(:,i) = feval(baseline,t);
end

figure
for i=1:n
    subplot(n,1,i)
    plot(t,y(:,i)-f(:,i))
end

figure
Nf = 4096;
for i=1:n
    subplot(n,1,i)
    pwelch([y(:,i) y(:,i)-f(:,i)],hanning(Nf),Nf/2,Nf,fs)
    xlim([0 5])
end

y0 = y-f;

%% Pt.1a : Plotting the signal
close all
clc

lead_names = {'V1','V2','V3','V4','V5','V6'};

FName = 'Times New Roman';
FSize = 12;
xl = [0 10]+10;
clr = lines(4);
delta = 8;

figure('Position',[100 100 900 900])
for i=1:n
    plot(t,y0(:,i)+delta*(i-1),'LineWidth',1.5,'Color',clr(1,:))
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
M = 46;
Orders = (1:2:M);

Optimize = true;

if Optimize
    IniGuess.Variances = [1e0 1e-12];
    IniGuess.TargetFrequency = 2*pi*1.07/fs;
    Niter = 100;
    [Modal,logMarginal,HyperPar,Initial] = MO_DSS_JointEKF_MultiHar_EM(y0(1001:2000,:)',Orders,Niter,IniGuess);
    
    save('Optimized_HyperPar','HyperPar','Initial','IniGuess','logMarginal','Niter')
    
else
    load('Optimized_HyperPar','HyperPar','IniGuess','logMarginal','Niter')
end

[Modal,logMarginal] = MO_DSS_JointEKF_MultiHar_Integrated(y0',Orders,HyperPar,Initial);

% [Modal,logMarginal,HyperPar] = MO_DSS_JointEKF_EM(y',M,Niter,IniGuess);


%% Pt.3 : Showing results
close all
clc

delta = 8;
lead = 4;

figure('Position',[100 100 900 900])
plot(t,y0(:,lead),'Color',clr(1,:),'LineWidth',1.5)   
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
print('Figures\ECGcontrol_Modes','-dpng','-r300')

%% Performance after EM optimization
close all
clc

err = y0' - HyperPar.Psi*Modal.ym;

%- Reconstruction error plot
figure('Position',[100 100 900 300])
bar(100*var(err,[],2)./var(y0)')
set(gca,'XTickLabel',lead_names)
set(gca,'FontName',FName,'FontSize',FSize)
ylabel({'Reconstruction error';'RSS/SSS (%)'})

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGcontrol_ErrorPerf','-dpng','-r300')

%- Optimized covariances plot
figure('Position',[100 400 900 360])
subplot(121)
bar(diag(HyperPar.Q(1:2:M,1:2:M))*1e2)
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
print('Figures\ECGcontrol_Covariances','-dpng','-r300')

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
print('Figures\ECGcontrol_MixMatrix','-dpng','-r300')

%%
close all
clc

leads = [1 4];
M_max = min(M,4);
yl = [-1 1; -1 1; -1 1; -1 1];

figure('Position',[100 100 900 800])
for j=1:2
    
    subplot(M_max+1,2,j)
    plot(t,y0(:,leads(j)),'Color',clr(1,:),'LineWidth',1.5)
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

set(gcf,'PaperPositionMode','auto')
print('Figures\ECGcontrol_HarmonicComponents','-dpng','-r300')

%%
close all
clc

r = 6;

figure('Position',[100 100 900 800])
plot3(zeros(1,N),t,y0(:,r),'LineWidth',1.5,'Color',clr(2,:))
hold on
for i=1:min(10,M)
    
    ind = (1:2) + 2*(i-1);
    psi = zeros(1,M);
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

r = 4;

psi = HyperPar.Psi(r,ind);
psi = sqrt( psi(1:2:end).^2 + psi(2:2:end).^2 );
A = diag(psi)*Modal.Am;
Mmax = 20;

figure('Position',[100 100 900 800])
subplot(4,1,4)
plot(t,y0(:,r))
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

%%
close all
clc

subplot(211)
plot(t,Modal.ym(M-1:M,:))

subplot(212)
plot(t,err)

figure
subplot(211)
pwelch(Modal.ym(M-1:M,:)')

subplot(212)
pwelch(err')

figure
for i=1:n
    subplot(3,2,i)
    plot(t,y0(:,i),t,y0(:,i)+err(i,:)')
    xlim(xl)
end