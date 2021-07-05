clear
close all
clc

addpath('..\Core\')

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

%% Projection decomposition
close all
clc

z = zeros(N,100);
for i=1:50
    ind = (1:2) + 2*(i-1);
    z(:,ind) = [cos(cumsum(i*omega_ref')) sin(cumsum(i*omega_ref'))];
end

theta = z\y;

figure
plot((theta(1:2:end,:).^2+theta(2:2:end,:).^2))