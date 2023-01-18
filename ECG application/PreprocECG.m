function [y,fs,t,lead_names,RRinterval] = PreprocECG(val,Fs,leads)

lead_names = {'i','ii','iii','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'};
lead_names = lead_names(leads);

fs = 250;
y = resample(val,fs,Fs);
y = y(1:8e3,leads);
y = detrend(y);
y = y/diag(std(y));
[N,~] = size(y);
t = (0:N-1)/fs;

%% RR interval estimation
close all
clc

[~,idx] = findpeaks(y(:,5),'MinPeakDistance',0.9,'MinPeakProminence',4);
locs = t(idx);

[RRinterval,T] = resample( diff(locs), idx(2:end), 1, 1, 1, 'spline' );

% figure
% plot(t,y(:,2))
% hold on
% plot(locs,y(idx,2),'o')

%% Remove baseline interference

idx_bl = idx - round(1*diff([0 idx'])/4)';
locs_bl = t(idx_bl);

f = zeros(size(y));
for i=1:6
    [baseline,~,~] = fit( locs_bl', y(idx_bl,i), 'cubicinterp' );
    f(:,i) = feval(baseline,t);
end

y = y-f;
y = y(T,:);
t = t(T);

y = y - repmat( mean(y), numel(T), 1 );