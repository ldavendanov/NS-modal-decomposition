clear
close all
clc

m = [20 0.5];               % Mass (kg)
k = [2000 6000 500];             % Stiffness (N/m)
c = 0.01*[1 0.1];                % Damping (N*s/m)


omegaN = [sqrt(k(1)/sum(m)) sqrt(k(2)/sum(m))];

OmegaRef = 2*pi*1.4;
T = 300;
fs = 200;
z0 = zeros(9,1);
t = 0:1/fs:T;

[t,z] = ode45( @(t,z)UnbalancedRotorODE4dof(t,z,OmegaRef), t, z0 );

%% Calculate accelerations from displacements and velocities
clc

y = zeros(numel(t),3);
for i=1:numel(t)
    dz = UnbalancedRotorODE4dof(t(i),z(i,:)',OmegaRef);
    y(i,:) = dz(5:7);
end
Omega = z(:,8);

save('URsystem_response','y','t','z','OmegaRef','Omega','fs')

%% Simulation results
close all
clc

Tplot = [200 240];

figure('Position',[100 100 900 600])
subplot(411)
plot(t,y(:,1)*1e3)
grid on
xlabel('Time (s)')
ylabel({'Stator horizontal';'acceleration (m/s^2)'})
xlim(Tplot)

subplot(412)
plot(t,y(:,2)*1e3)
grid on
xlabel('Time (s)')
ylabel({'Stator vertical';'acceleration (m/s^2)'})
xlim(Tplot)

subplot(413)
plot(t,y(:,3))
grid on
xlabel('Time (s)')
ylabel({'Rotor radial';'acceleration (m/s^2)'})
xlim(Tplot)

subplot(414)
plot(t,Omega/(2*pi))
hold on
plot([t(1) t(end)],OmegaRef*[1 1]/(2*pi),'-k')
text(Tplot(1),OmegaRef/(2*pi),'\Omega_d','HorizontalAlignment','left','VerticalAlignment','bottom')
grid on
xlabel('Time (s)')
ylabel('Rotor speed (Hz)')
xlim(Tplot)

print('Figures\UR_SampleResponseTime','-r300','-dpng')

figure('Position',[100 600 900 200])
plot(t,z(:,9))
grid on
xlabel('Time (s)')
ylabel('Applied torque (N/m)')
xlim(Tplot)

%% Frequency domain analysis
close all
clc

clr = lines(3);
yl = 1e-6*[1e-2 1e6];

Y = fft(y(2002:end,:));
Y = Y / size(Y,1);
ff = fs*(0:size(Y,1)-1)/size(Y,1);

figure('Position',[100 100 900 300])
for i=1:8
    semilogy(OmegaRef/(2*pi)*[1 1]*i,yl,'--k')
    hold on
    if i>1
        text(OmegaRef/(2*pi)*i,yl(2),[num2str(i),'\Omega_d'],'HorizontalAlignment','center','VerticalAlignment','bottom')
    else
        text(OmegaRef/(2*pi)*i,yl(2),'\Omega_d','HorizontalAlignment','center','VerticalAlignment','bottom')
    end
end

for i=1:2
    semilogy(omegaN(i)/(2*pi)*[1 1],yl,'-.','Color',0.25*[1 1 1])
    hold on
    text(omegaN(i)/(2*pi),yl(2),['\omega_{n',num2str(i),'}'],'HorizontalAlignment','center','VerticalAlignment','bottom','Color',0.25*[1 1 1])
end

p(1) = semilogy(ff,abs(Y(:,1)),'LineWidth',1.5,'Color',clr(1,:));
p(2) = semilogy(ff,abs(Y(:,2)),'LineWidth',1.5,'Color',clr(2,:));
p(3) = semilogy(ff,abs(Y(:,3)),'LineWidth',1.5,'Color',clr(3,:));
xlim([0 18]), ylim(yl)
xlabel('Frequency (Hz)')
ylabel('Acceleration magnitude (m/s^2)')
grid on

SensorLabel = {'Stator horizontal acceleration';'Stator vertical acceleration';'Rotor radial acceleration'};
legend(p,SensorLabel)

print('Figures\UR_SampleSignalFreq','-r300','-dpng')

%% Time-frequency analysis
close all
clc

cmap = OriBlue;

nf = 2^11;
nover = nf-4;
alpha = 4;
[Syy{1},ff,tt] = spectrogram(y(4001:12000,1),gausswin(nf,alpha),nover,nf,fs);
Syy{2} = spectrogram(y(4001:12000,2),gausswin(nf,alpha),nover,nf,fs);
Syy{3} = spectrogram(y(4001:12000,3),gausswin(nf,alpha),nover,nf,fs);

figure('Position',[100 100 900 540])
for i=1:3
    subplot(1,3,i)
    imagesc(tt,ff,20*log10(abs(Syy{i})))
    ylim([0 15])
    axis xy
    colormap(cmap)
    grid on
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title(SensorLabel{i})
    cbar = colorbar; 
%     set(gca,'CLim',[-120 40])
    cbar.Label.String = 'Power (dB)';
    hold on
    for j=1:8
        plot(t(4001:12000)-t(4001),j*Omega(4001:12000)/(2*pi),'-k')
    end
end

print('Figures\UR_SampleSignal_TF','-r300','-dpng')