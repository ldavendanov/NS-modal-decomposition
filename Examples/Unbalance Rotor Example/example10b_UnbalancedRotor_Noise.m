clear
close all
clc

OmegaRef = 2*pi*2;
T = 300;
fs = 60;
z0 = zeros(7,1);
t = 0:1/fs:T;

[t,z] = ode45( @(t,z)UnbalancedRotorODE(t,z,OmegaRef), t, z0 );

%% Calculate accelerations from displacements and velocities
clc

y = zeros(numel(t),2);
for i=1:numel(t)
    dz = UnbalancedRotorODE(t(i),z(i,:)',OmegaRef);
    y(i,:) = dz(4:5);
end

%% Simulation results
close all
clc

figure
subplot(511)
plot(t,y(:,1)*1e3)
grid on
xlabel('Time (s)')
ylabel('Engine acceleration (m/s^2)')
xlim([160 180])

subplot(512)
plot(t,y(:,2))
grid on
xlabel('Time (s)')
ylabel('Rotor radial acceleration (m/s^2)')
xlim([160 180])

subplot(513)
plot(t,z(:,6))
hold on
plot([t(1) t(end)],OmegaRef*[1 1])
grid on
xlabel('Time (s)')
ylabel('Rotor speed (rad/s)')
xlim([160 180])

subplot(514)
plot(t,OmegaRef-z(:,6))
grid on
xlabel('Time (s)')
ylabel('Speed error (rad/s)')
xlim([160 180])

subplot(515)
plot(t,z(:,7))
grid on
xlabel('Time (s)')
ylabel('Applied torque (N/m)')
xlim([160 180])


Nf = 4096;
[Pyy,frec] = pwelch(z(:,1:2),hann(Nf),Nf/2,Nf,fs);

figure
plot(frec,10*log10(Pyy))
yl = get(gca,'YLim');
hold on
for i=1:7
    plot(i*[1 1]*OmegaRef/(2*pi),yl,'--k')
end

figure
for i=1:2
    subplot(1,2,i)
    spectrogram(y(:,i),gausswin(512,4),510,512,fs)
end

%% Optimize harmonic decomposition
close all
clc

Y = y(16001:17000,:);
time = t(1:1000);

addpath('..\Core\')


Orders = 1:6;
Niter = 50;
InitialGuess.TargetFrequency = 2/fs;
InitialGuess.Variances = [1e-1 1e-3];
[Modal,logMarginal,HyperPar,Initial] = MO_DSS_JointEKF_MultiHar_EM(Y',Orders,Niter,InitialGuess);

%%
close all
clc

SNR = 12:2:30;
logL = zeros(size(SNR));
ym = zeros(2*numel(Orders),size(Y,1),numel(SNR));
fN = zeros(size(Y,1),numel(SNR));
for i=1:numel(SNR)

    sigmaY2 = var(Y);
    sigmaN2 = sigmaY2*10.^(-SNR(i)/10);

    Yn = Y + randn(size(Y))*diag(sigmaN2);

    [Modal,logMarginal] = MO_DSS_JointEKF_MultiHar(Yn',Orders,HyperPar,Initial);
    logL(i) = sum(logMarginal);
    ym(:,:,i) = Modal.ym;
    fN(:,i) = Modal.omega*fs/(2*pi);
end

%%
close all
clc

figure
bar(SNR,logL)

figure
for i=1:4
    subplot(5,1,i)
    plot(time,squeeze(ym(2*i,:,:)))
    ylim(6*[-1 1])
end

subplot(515)
plot(time,abs(fN))

%%
function dz = UnbalancedRotorODE(~,z,OmegaRef)

% System settings
m = [10 0.4];               % Mass (kg)
k = [1000 200];             % Stiffness (N/m)
c = 0.01*[1 0.1];                % Damping (N*s/m)
l2 = 0.2;                   % Rotor radius before deformation (m)

% System variables
Omega = z(6);               % Rotor speed (rad/s)
theta = z(3);               % Rotor angle (rad)
d2 = l2 + z(2);             % Deformed rotor radius (m)
Fc = m(2)*d2*Omega^2;       % Centrifugal force of rotor (N)

% Rotor speed control
Cp = 4e1;                   % Proportional coefficient
Ci = 1e2;                   % Integral coefficient
Tau = z(7);                 % Applied torque
err = OmegaRef - Omega;

% System matrices
M = [            sum(m) m(2)*sin(theta) m(2)*d2*cos(theta);
        m(2)*sin(theta)            m(2)                  0;
     m(2)*d2*cos(theta)               0          m(2)*d2^2];

f = [ Fc*sin(theta) - 2*m(2)*cos(theta)*z(5)*z(6) - k(1)*z(1) - c(1)*z(4); 
      Fc - k(2)*z(2) - c(2)*z(5);
      Tau - 2*m(2)*z(5)*z(6)*d2];

dx = M\f;
dz = [z(4:6); dx; Ci*err - Cp*dx(3)];

end