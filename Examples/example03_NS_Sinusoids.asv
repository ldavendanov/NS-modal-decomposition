clear
close all
clc

addpath('..\Core\')

%% Pt.1 : Creating the simulation signal
close all
clc

%-- Signal features
N = 6e3;                                                                    % Number of samples
fs = 1e3;                                                                   % Sampling frequency
t = 0:N-1;                                                                  % Time vector

Aii = [1 0.2; 0.8 0.4; -0.2 1.0];                                           % Amplitude coefficients
fii = [50 -20; 80 20; 120 20];                                              % Frequency coefficients
alpha = [0.5 0.5 0.25];                                                     % Cyclic frequency
Psi =  [1.0 0.8 0.5 0.5 0.2 0.8;
        0.5 0.5 0.8 1.0 0.4 0.6;
        0.1 0.4 0.9 0.2 0.8 0.2];
for i=1:3
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

figure
subplot(411)
plot(t,IA)
ylim([0 2])

subplot(412)
plot(t,IF)

subplot(413)
plot(t,ym)

subplot(414)
plot(t,y)


%% Pt.2 : Estimating the signal components with the diagonal SS method
close all
clc

M = 3;
Niter = 10;
[Modal,logMarginal,HyperPar] = MO_DSS_JointEKF_EM(y,M,Niter);

%% Pt.3 : Showing results
close all
clc

figure('Position',[100 100 1200 600])
for i=1:M
    
    subplot(3,3,3*i-2)
    plot(t,ym(2*i,:)/std(ym(2*i,:)))
    hold on
    plot(t,Modal.ym(2*i,:)/std(Modal.ym(2*i,:)))
    grid on
    xlabel('Time [s]')
    ylabel(['Mode ',num2str(i)])
    
    subplot(3,3,3*i-1)
    plot(t,IA(i,:)/max(IA(i,100:end)))
    hold on
    plot(t,Modal.Am(i,:)/max(Modal.Am(i,100:end)))
    grid on
    xlabel('Time [s]')
    ylabel(['IA Mode ',num2str(i),' [--]'])
    
    subplot(3,3,3*i)
    plot(t,IF(i,:))
    hold on
    plot(t,Modal.omega(i,:)*fs/(2*pi))
    grid on
    xlabel('Time [s]')
    ylabel(['IF Mode ',num2str(i),' [Hz]'])
    
end