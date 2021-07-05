clear
close all
clc

addpath('..\Core\')
FName = 'Times New Roman';
FSize = 12;

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
[N,n] = size(y);
t = t-12;

%% Singular spectrum analysis

na = 20;
tau = na+1:N;

Phi = zeros(N-na,na*n);
for i=1:na
    ind = (1:n) + n*(i-1);
    Phi(:,ind) = y(tau-i+1,:);
end

[U,S,V] = svd(Phi,'econ');
Psi = S*V(:,1:n);
lambda = diag(S).^2;
Ncomp = find( cumsum(lambda) > 0.99*sum(lambda), 1, 'first');

plot(cumsum(lambda)/sum(lambda))

%% Plotting results
close all
clc

xl = [0 2];
lead = 2;
delta = 3;
clr = lines(4);

figure('Position',[100 100 300 700])
for i=1:6
    plot(t(tau),Psi(i,lead)*U(:,i)+delta*(i-1),'Color',clr(1,:))
    hold on
    text(xl(1)+0.05,delta*(i-1)+0.5,['Cmp. ',num2str(i)],'FontName',FName,'FontSize',FSize-2)
    xlim(xl)
end
set(gca,'FontName',FName,'FontSize',FSize,'YTick',[])
ylim([-2 17])
xlabel('Time [s]')

set(gcf,'PaperPositionMode','auto')
print(['Figures\ECGhealthy_SSAmodes',lead_names{lead}],'-dpng','-r300')

%% Discrete wavelet transform
close all
clc

lead = 2;
lvl = 5;
delta = 4; 

% [c,l] = wavedec(y(:,lead),lvl,'sym4');
[c,l] = wavedec(y(:,lead),lvl,'db6');
n_idx = [0; cumsum(l)];
n_idx = n_idx(1:end-1);
z = zeros(N,lvl+1);

for i=1:lvl+1
    c0 = zeros(size(c));
    c0(n_idx(i)+1:n_idx(i+1)) = c(n_idx(i)+1:n_idx(i+1));
    z(:,i) = waverec( c0, l, 'sym4' );    
end

figure('Position',[100 100 300 700])
for i=1:lvl+1
    plot(t,z(:,i)+delta*(i-1),'Color',clr(1,:))
    hold on
    text(xl(1)+0.05,delta*(i-1)+0.75,['Cmp. ',num2str(i)],'FontName',FName,'FontSize',FSize-2)
    xlim(xl)
end
set(gca,'FontName',FName,'FontSize',FSize,'YTick',[])
ylim([-2 22])
xlabel('Time [s]')

set(gcf,'PaperPositionMode','auto')
print(['Figures\ECGhealthy_DWTmodes',lead_names{lead}],'-dpng','-r300')

%% Empirical mode decomposition
close all
clc

lead = 2;

delta = 5;
[z,residual] = emd(y(:,lead),'Interpolation','pchip');
z = [z residual];
n_comp = 6;

figure('Position',[100 100 300 700])
for i=1:n_comp
    plot(t,z(:,i)+delta*(i-1),'Color',clr(1,:))
    hold on
    text(xl(1)+0.05,delta*(i-1)+1,['Cmp. ',num2str(i)],'FontName',FName,'FontSize',FSize-2)
    xlim(xl)
end
set(gca,'FontName',FName,'FontSize',FSize,'YTick',[])
ylim([-4 28])
xlabel('Time [s]')

set(gcf,'PaperPositionMode','auto')
print(['Figures\ECGhealthy_EMDmodes',lead_names{lead}],'-dpng','-r300')
