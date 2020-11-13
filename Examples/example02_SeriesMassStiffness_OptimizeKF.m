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

%% Pt 3 : Optimize the KF with multiple output ----------------------------
close all
clc

% Performing the calculation
M = 3;
Niter = 10;
[Modal,logMarginal,HyperPar] = MO_DSS_JointEKF_EM(y,M,Niter);

%%
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
    plot(t,Modal.omega(i,:)*fs/(2*pi),'Color',clr(i,:))
    hold on
    plot(t,fN(i,:),'--','Color',clr(i,:))
    grid on
    ylabel('Frequency [Hz]')
    xlabel('Time [s]')
end
xlim(TT)
ylim([0 fs/2])

title('Estimate based on the three outputs')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%--- Other functions ---------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = StiffnessMatrix(k)
% Builds the stiffness matrix of a series spring-mass system

K = [k(1)+k(2)     -k(2)     0;
         -k(2) k(2)+k(3) -k(3);
             0     -k(3)  k(3)];

end

%--------------------------------------------------------------------------
function dx = MassSpringSystem( t, x, f, fs )

[M,K,C] = SystemMatrices( t );
n = size(M,1);
A = [zeros(n) eye(n);
     -M\K -M\C];
B = [zeros(n,1); -M\[0 0 1]'];

k = max(round(t*fs),1);
dx = A*x + B*f(k);

end

%--------------------------------------------------------------------------
function [M,K,C] = SystemMatrices( t )

%-- System properties
m = [1 1 1];                % Masses (kg)
k0 = 1000*[2 2 2];          % Stiffness - constant component
k1 = 1000*[0  0.8 0];       % Stiffness - cosine component
k2 = 1000*[0 -0.4 0];       % Stiffness - sine component
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

%--------------------------------------------------------------------------
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