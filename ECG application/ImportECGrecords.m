clear
close all
clc

% Control patient
[sig, Fs, tm] = rdsamp('ptbdb/patient104/s0306lre');


%%
close all
clc

plot(tm,sig)