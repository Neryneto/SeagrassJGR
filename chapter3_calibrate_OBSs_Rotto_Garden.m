path='C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\GardenIsland_exp1\data\sawhorse';

OBS_Bot = ncreadall (fullfile(path,'\OBS_BOT\obs_3194_V0.nc'));
OBS_Mid = ncreadall (fullfile(path,'\OBS_MID\obs_2997_V0.nc')); 
OBS_Top = ncreadall (fullfile(path,'\OBS_TOP\obs_2535_V0.nc')); 

[Sig, TStr, Raw] = xlsread ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\GardenIsland_exp1\data\SSC\pump_samples\GardenIslandPumpSamples.xlsx','Sheet1','F3:M32');
time_pump_bot = Sig(:,1)+datenum('30-Dec-1899');
[Sig, TStr, Raw] = xlsread ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\GardenIsland_exp1\data\SSC\pump_samples\GardenIslandPumpSamples.xlsx','Sheet1','F33:M61');
time_pump_mid = Sig(:,1)+datenum('30-Dec-1899');
%%

% for i=1:length(time_pump_bot)
[OBS_Bot.time_5min, OBS_Bot.avg_5min] = create_obs_5_min_interval(OBS_Bot, 0.242, 50);
OBS_Bot.avg_5min(OBS_Bot.avg_5min<0 | OBS_Bot.avg_5min>70)=nan;
[OBS_Mid.time_5min, OBS_Mid.avg_5min] = create_obs_5_min_interval(OBS_Mid, 0.242, 50);
OBS_Mid.avg_5min(OBS_Mid.avg_5min<0 | OBS_Mid.avg_5min>70)=nan;
[OBS_Top.time_5min, OBS_Top.avg_5min] = create_obs_5_min_interval(OBS_Top, 0.0025, 60);
