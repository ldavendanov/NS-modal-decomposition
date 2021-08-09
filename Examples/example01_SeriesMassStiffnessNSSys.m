clear
close all
clc

addpath('..\Core\')

%% Pt 1 : Creating a NS system --------------------------------------------
close all
clc

%-- Simulation parameters
T = 100;                                                                    % Simulation period (seconds)
fs = 50;                                                                    % Sampling frequency (Hz)
t = 0:1/fs:T;                                                               % Time vector
N = numel(t);                                                               % Number of samples
n = 3;                                                                      % Number of DOFs

%-- Frozen modal analysis
[fN,zeta,Phi] = FrozenModal(T,fs);                                          % Calculate the frozen modal analysis

figure('Position',[100 100 1200 450])
subplot(121)
plot(t,fN)
grid on
xlabel('Time [s]')
ylabel('Natural frequency [Hz]')

subplot(122)
plot(t,100*zeta)
grid on
xlabel('Time [s]')
ylabel('Damping ratio [%]')

figure('Position',[1300 100 600 600])
for i=1:n
    subplot(n,1,i)
    imagesc( t, 1:n, squeeze(real(Phi(:,i,:))) )
    colorbar
end

%% Pt 2 : Make a simulation of the NS system ------------------------------
close all
clc

% Simulate the system's response
x0 = zeros(2*n,1);                                                          % Initial condition
F = randn(1,N);                                                             % Force excitation
[~,x] = ode45(@(t,x)MassSpringSystem(t,x,F,fs),t,x0);                       % Simulate the system response
y = [zeros(1,n); diff(x(:,4:6))]';                                          % Acceleration signal

% Calculate the spectrogram
Nf = 512;
[Sxx,ff,tt] = spectrogram(y(1,:),gausswin(Nf),Nf-2,Nf,fs);

figure('Position',[100 100 1200 640])
subplot(311)
plot(t,y)
ylabel('Acceleration [m/s^2]')
legend({'x_1','x_2','x_3'},'Location','north','Orientation','horizontal')
grid on

subplot(3,1,[2 3])
imagesc(tt,ff,10*log10(abs(Sxx)))
axis xy
hold on
plot(t,fN,'r')
xlim([0 T])
xlabel('Time [s]')
ylabel('Frequency [Hz]')
grid on
legend('Frozen natural frequencies')

%% Pt 3 : Modal decomposition based on system matrices --------------------
close all
clc

TT = [0 T];

% Projecting to obtain the modal decomposition
yN = zeros(n,N);
for i=1:N
    yN(:,i) = Phi(:,:,i)'*y(:,i);
end

figure('Position',[100 100 900 600])
for i=1:n
    subplot(n,1,i)
    plot(t,real(yN(i,:)))
    xlim(TT)
    xlabel('Time [s]')
    ylabel(['Mode No. ',num2str(i)])
    grid on
end

%% Pt 3 : KF-based modal decomposition of the response --------------------
close all
clc

% Setting up the model parameters
M = 3;                                                                      % Number of modes
variances = [1 1e-1 1e-2];                                                  % KF variances

% Setting up initial values of the state vector
x0 = zeros(4*n,1);
for i=1:n
    x0(2*n+2*i-1) = cos(2*pi*fN(i,1)/fs);
    x0(2*n+2*i) = sin(2*pi*fN(i,1)/fs);
end

% Initialize computation matrices
Am = cell(3,1);                                                             % IA estimates
ym = cell(3,1);                                                             % Mode trajectories
fm = cell(3,1);                                                             % IF estimates

TT = [20 40];

for j=1:n    
    
    [ym{j},Am{j},omega] = DiagonalSS1(y(j,:),x0,M,variances);
    fm{j} = omega*fs/(2*pi);
    
    clr = lines(M);
    
    figure('Position',[100 100 1200 600])
    for i=1:M
        subplot(M,3,3*i-2)
        plot(t,ym{j}(i,:),'Color',clr(i,:))
        xlim(TT)
        grid on
        ylabel(['Mode ',num2str(i)])
        xlabel('Time [s]')
        
        subplot(M,3,3*i-1)
        plot(t,Am{j}(i,:),'Color',clr(i,:))
        xlim(TT)
        grid on
        ylabel(['Amplitude mode ',num2str(i)])
        xlabel('Time [s]')
    end
    
    subplot(1,3,3)
    for i=1:M
        plot(t,fm{j}(i,:),'Color',clr(i,:))
        hold on
        plot(t,fN(i,:),'--','Color',clr(i,:))
        grid on
        ylabel('Frequency [Hz]')
        xlabel('Time [s]')
    end
    xlim(TT)
    ylim([0 fs/2])
    
    title(['Estimate based on output ',num2str(j)])
    
end

%% Comparison of the frozen and KF-based modes ----------------------------
close all
clc

TT = [0 50];

figure('Position',[100 100 900 600])
for i=1:M
    subplot(M,1,i)
    for j=1:n
        plot(t,ym{j}(i,:)/norm(ym{j}(i,:)))
        hold on        
    end
    plot(t,real(yN(i,:)/norm(real(yN(i,:)))))
    xlim(TT)
end


%% Pt 4 : Modal decomposition based on KF with multiple output ------------
close all
clc

% Setting up the measurement matrix
Psi = zeros(n,2*n);
Psi(:,1:2:end) = real(Phi(:,:,1));
Psi(:,2:2:end) = imag(Phi(:,:,1));

% Setting up the initial state vector
x0 = zeros(4*n,1);
for i=1:n
    x0(2*n+2*i-1) = cos(2*pi*fN(i,1)/fs);
    x0(2*n+2*i) = sin(2*pi*fN(i,1)/fs);
end

% Setting the KF variances
variances = [1e-4 1e-4 1e-6];

Niter = 100;
IniVal.Psi = Psi;
IniVal.x0 = x0;
IniVal.Variances = variances;
[Modal,logMarginal,HyperPar,Initial] = MO_DSS_JointEKF_EM(y,M,Niter,IniVal);

%% Plotting results
close all
clc

TT = [0 100];

clr = lines(M);

figure('Position',[100 100 1200 600])
for i=1:M
    subplot(M,3,3*i-2)
    plot(t,Modal.ym(2*i,:),'Color',clr(i,:))
    xlim(TT)
    grid on
    ylabel(['Mode ',num2str(i)])
    xlabel('Time [s]')
    
    subplot(M,3,3*i-1)
    plot(t,Modal.Am(i,:),'Color',clr(i,:))
    xlim(TT)
    grid on
    ylabel(['Amplitude mode ',num2str(i)])
    xlabel('Time [s]')
end

subplot(1,3,3)
for i=1:M
    plot(t,omega(i,:)*fs/(2*pi),'Color',clr(i,:))
    hold on
    plot(t,Modal.omega(i,:)*fs/(2*pi),'--','Color',clr(i,:))
    grid on
    ylabel('Frequency [Hz]')
    xlabel('Time [s]')
end
xlim(TT)
ylim([0 fs/2])

title('Estimate based on the three outputs')

%%
close all
clc

XLim = [40 50];

figure('Position',[100 100 900 800])

for i=1:3
    
    subplot(3,2,2*i-1)
    plot(t,y(i,:))
    hold on
    plot(t,HyperPar.Psi(i,:)*Modal.ym)
    xlim(XLim)
    
    subplot(3,2,2*i)
    plot(t,fN(i,:))
    hold on
    plot(t,omega(i,:)*fs/(2*pi))
    xlim(XLim)
end

%%
close all
clc

TT = [20 80];

figure('Position',[100 100 900 600])
for i=1:M
    subplot(M,1,i)
    plot(t,Modal.ym(2*i,:)/norm(Modal.ym(2*i,:)))
    hold on
    plot(t,real(yN(i,:)/norm(real(yN(i,:)))))
    xlim(TT)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%--- Other functions ---------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = StiffnessMatrix(k)
% Builds the stiffness matrix of a series spring-mass system

K = [k(1)+k(2)     -k(2)     0;
         -k(2) k(2)+k(3) -k(3);
             0     -k(3)  k(3)];

end

function dx = MassSpringSystem( t, x, f, fs )

[M,K,C] = SystemMatrices( t );
n = size(M,1);
A = [zeros(n) eye(n);
     -M\K -M\C];
B = [zeros(n,1); -M\[0 0 1]'];

k = max(round(t*fs),1);
dx = A*x + B*f(k);

end

function [M,K,C] = SystemMatrices( t )

%-- System properties
m = [1 1 1];                % Masses (kg)
k0 = 1000*[2 2 2];          % Stiffness - constant component
k1 = 1000*[0  0.5 0];       % Stiffness - cosine component
k2 = 1000*[0 -0.2 0];       % Stiffness - sine component
Omega = 2*pi*1/20;           % Period of stiffness variation
alpha = 2e-4;
beta = 2e-4;

%-- Building system matrices
M = diag(m);
K0 = StiffnessMatrix(k0);
K1 = StiffnessMatrix(k1);
K2 = StiffnessMatrix(k2);

K = K0 + K1*cos( Omega*t ) + K2*sin( Omega*t );
C = alpha*K0 + beta*M;

end

function [fN,zeta,Phi] = FrozenModal(T,fs)

n = 3;                                                                      % Number of DOFs
N = T*fs+1;                                                                 % Number of samples
t = (0:N-1)/fs;                                                             % Time vector

% Initializing computation matrices
fN = zeros(2*n,N);
zeta = zeros(2*n,N);
V = zeros(2*n,2*n,N);

% Calculating the frozen eigenanalysis
for i=1:N
    [M,K,C] = SystemMatrices( t(i) );                                       % Extract system matrices
    A = [zeros(n) eye(n); -M\K -M\C];                                       % Calculate the SS matrix
    [V(:,:,i),E] = eig(A);                                                  % Eigenvalue decomposition
    fN(:,i) = abs(diag(E))/(2*pi);                                          % Natural frequencies
    zeta(:,i) = -cos(angle(diag(E)));                                       % Damping ratios
end

% Extracting and sorting modes
fN = fN(1:2:end,:);
zeta = zeta(1:2:end,:);
[~,ind] = sort(mean(fN,2));
fN = fN(ind,:);
zeta = zeta(ind,:);

% Fixing the frozen mode shape matrix
Phi = zeros(n,n,N);
for i=1:N
    for j=1:n
        Phi(:,j,i) = V(1:n,2*ind(j)-1,i);
        Phi(:,j,i) = Phi(:,j,i)*sign(real(Phi(n,j,i)))/abs(Phi(n,j,i));
    end
end

end