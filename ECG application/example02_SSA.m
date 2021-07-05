clear
close all
clc

addpath('..\Core\')

%% Pt.1 : Loading ECG recording

fs = 100;
[val, Fs] = rdsamp('ptbdb/patient104/s0306lre');

y = resample(val,fs,Fs);
y = y(1:4e3,7:12);
y = detrend(y);
y = y/diag(std(y));
[N,n] = size(y);
t = (0:N-1)/fs;

%% Pt.2 : Calculating SSA

na = 800;
tau = na+1:N;
Phi = zeros(N-na,n*na);
for i=1:na
    Phi(:,(1:n)+n*(i-1)) = y(tau-i,:);
end

[U,S,V] = svd(Phi,'econ');
lambda = diag(S).^2;
Npc = find(cumsum(lambda)/sum(lambda)>0.99,1,'first');

%% Pt.3 : Plotting results
close all
clc

Nf = 1024;
[Puu,f] = pwelch(U,hanning(Nf),3*Nf/4,Nf,fs);

figure
plot(cumsum(lambda)/sum(lambda))

figure
for i=1:12
    subplot(4,3,i)
    plot(U(:,i))
end

figure
for i=1:12
    subplot(4,3,i)
    plot(f,Puu(:,i))
end
    

%%
close all
clc

Uh = hilbert(U);
IA = abs(Uh);
IF = abs( diff(unwrap(angle(Uh)))*(fs/(2*pi)) );

figure
subplot(211)
plot(tau/fs,IA(:,2:11))
xlim([8 10])

subplot(212)
plot(tau(2:end)/fs,IF(:,2:11))
% ylim([0 20])
xlim([8 10])

%%
close all
clc

Nf = 2048;
g = gausswin(Nf,4);
[Syy,ff,tt] = spectrogram(y(:,1),g,Nf-2,Nf,fs);
Syy = abs(Syy).^2 * ( var(g)/Nf );

figure
subplot(411)
plot(t,y(:,1))
xlim([20 22])
grid on

subplot(4,1,2:4)
imagesc(tt,ff,10*log10(Syy))
hold on
plot(tau(2:end)/fs,IF(:,1:10),'k')
ylim([0 20])
xlim([20 22])
axis xy
grid on
