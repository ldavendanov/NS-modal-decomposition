function [y,tacho,Fs] = LoadRigVibration(state,f_hss,rig_load)

BaseFolder = 'C:\Users\ldav\OneDrive - Syddansk Universitet\SDU Lectures\2020 Autumn - Modeling and design of rotating machines\Matlab Exercises\';

% -- Loading data according to input
if rig_load == 0
    i = f_hss - 10 + 1;
else
    i = f_hss - 10 + 12;
end
switch state
    case 'H'
        tacho = load([BaseFolder,'Data\Healthy\H',num2str(i,'%02d'),'_',num2str(f_hss,'%02d'),'Hz_',num2str(rig_load,'%02d'),'L\channel1']);
        data_vert = load([BaseFolder,'Data\Healthy\H',num2str(i,'%02d'),'_',num2str(f_hss,'%02d'),'Hz_',num2str(rig_load,'%02d'),'L\channel2']);
        data_horz = load([BaseFolder,'Data\Healthy\H',num2str(i,'%02d'),'_',num2str(f_hss,'%02d'),'Hz_',num2str(rig_load,'%02d'),'L\channel3']);
    case 'D'
        tacho = load([BaseFolder,'Data\Damage_OuterRing\D',num2str(i,'%02d'),'_',num2str(f_hss,'%02d'),'Hz_',num2str(rig_load,'%02d'),'L\channel1']);
        data_vert = load([BaseFolder,'Data\Damage_OuterRing\D',num2str(i,'%02d'),'_',num2str(f_hss,'%02d'),'Hz_',num2str(rig_load,'%02d'),'L\channel2']);
        data_horz = load([BaseFolder,'Data\Damage_OuterRing\D',num2str(i,'%02d'),'_',num2str(f_hss,'%02d'),'Hz_',num2str(rig_load,'%02d'),'L\channel3']);
end

% -- Extracting desired outputs
tacho = tacho.Data;
y = [data_horz.Data data_vert.Data];
Fs = 1/data_horz.Header.xIncrement;

% -- Some signal conditioning
y = detrend(y);