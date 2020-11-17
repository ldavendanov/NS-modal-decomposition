clear
close all
clc

addpath('..\Core\')

%% Pt.1 : Creating the simulation signal
close all
clc

%-- Signal features
N = 6e3;                                                                    % Number of samples
fs = 5e2;                                                                   % Sampling frequency
t = 0:N-1;                                                                  % Time vector
SNR = 60;                                                                   % Signal to noise ratio

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
sigmaN = sqrt( min(var(y,[],2))*10.^( -SNR/10 ) );
y = y + sigmaN*randn(3,N);

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
Niter = 60;
IniGuess.Variances = [1e-6 1e-6];
IniGuess.TargetFreqs = 2*pi*IF(:,1)/fs;
[Modal,logMarginal,HyperPar] = MO_DSS_JointEKF_EM(y,M,Niter,IniGuess);

%% Pt.3 : Showing results
close all
clc

xl = [0.0 4.0];
clr = lines(2);

figure('Position',[100 100 1200 600])
for i=1:M
    
    subplot(3,3,3*(M-i+1)-2)
    plot(t/fs,ym(2*i-1,:)/std(ym(2*i,:)))
    hold on
    plot((t)/fs,Modal.ym(2*i,:)/std(Modal.ym(2*i,:)))
    grid on
    xlabel('Time [s]')
    ylabel(['Mode ',num2str(i)])
    set(gca,'FontName','Times New Roman','FontSize',12)
    xlim(xl)
    
    subplot(3,3,3*(M-i+1)-1)
    plot(t/fs,IA(i,:)/max(IA(i,100:end)))
    hold on
    plot(t/fs,Modal.Am(i,:)/max(Modal.Am(i,100:end)))
    grid on
    xlabel('Time [s]')
    ylabel(['IA Mode ',num2str(i),' [--]'])
    set(gca,'FontName','Times New Roman','FontSize',12)
    xlim(xl)
    
    subplot(1,3,3)
    plot(t/fs,IF(i,:),'Color',clr(1,:))
    hold on
    plot(t/fs,Modal.omega(i,:)*fs/(2*pi),'Color',clr(2,:))
    grid on
    xlabel('Time [s]')
    ylabel(['IF Mode ',num2str(i),' [Hz]'])
    set(gca,'FontName','Times New Roman','FontSize',12)
    xlim(xl)
    
end

Qth = HyperPar.Q(2*M+1:end,2*M+1:end);
R = Qth;
for i=1:2*M
    c = sqrt(Qth(i,i))*sqrt( diag( Qth )' );
    R(i,:) = R(i,:)./c;
end

figure
imagesc(R)
colorbar
set(gca,'CLim',[-1 1])

%% Plotting results - Hyperparameter estimates
close all
clc

q = diag(HyperPar.Q);

figure('Position',[100 100 1200 400])
subplot(131)
bar(log10(q(1:2*M)))
xlabel('Index')
ylabel('$\log_{10} \sigma_{z_{i}}^2$','Interpreter','latex')
title('State innov. variance')
set(gca,'FontName','Times New Roman','FontSize',12)

subplot(132)
bar(log10(q(2*M+1:end)))
xlabel('Index')
ylabel('$\log_{10} \sigma_{\theta_{i}}^2$','Interpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',12)
title('Parameter innov. variance')

subplot(133)
bar(log10(diag(HyperPar.R)))
xlabel('Index')
ylabel('$\log_{10}\sigma_{y_{i}}^2$','Interpreter','latex')
title('Measurement noise variance')
set(gca,'FontName','Times New Roman','FontSize',12)

figure('Position',[100 500 800 400])
imagesc(abs(HyperPar.Psi))
xlabel('Mode index')
ylabel('Output index')
set(gca,'XTick',1:2*M,'YTick',1:3,'CLim',[0 inf])
set(gca,'FontName','Times New Roman','FontSize',12)
colorbar
title('Mixing matrix - Absolute value')

R = HyperPar.Q(2*M+1:end,2*M+1:end);
for i=1:2*M
    R(i,:) = R(i,:)./diag( HyperPar.Q(2*M+1:end,2*M+1:end) )';
end

figure
imagesc(R)
colorbar
set(gca,'CLim',[-1 1])

%%
close all
clc

figure
for i=1:3
    
    psi1 = Psi(i,1:2:end).^2 + Psi(i,2:2:end).^2;
    psi2 = HyperPar.Psi(i,1:2:end).^2 + HyperPar.Psi(i,2:2:end).^2;
    subplot(3,1,i)
    plot( psi1/norm(psi1) )
    hold on
    plot( psi2/norm(psi2) )
    ylim([0 1])
    xlabel('Output index')
    ylabel(['|\psi_',num2str(i),'|'])
    set(gca,'FontName','Times New Roman','FontSize',12)

    
end