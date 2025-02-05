%open Aquadopp
caminho = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\2.Part_2\Quadpod';

AQD11898_P2 = ncreadall(fullfile(caminho,'AQD11898\AQD11898_V0_Part2'));
for i = 1:size(AQD11898_P2.VEL_XE_RAW, 2)
    [AQD11898_P2.Despike_vel_x(:,i), AQD11898_P2.Despike_vel_y(:,i), AQD11898_P2.Despike_vel_z(:,i)] = func_despike_phasespace3d_3var_ADP (AQD11898_P2.VEL_XE_RAW(:,i),...
        AQD11898_P2.VEL_YN_RAW(:,i), AQD11898_P2.VEL_ZU_RAW(:,i), 0);
end

[AQD11898_P2.vel_maj, AQD11898_P2.vel_min, AQD11898_P2.Mean_Current...
    AQD11898_P2.Current_theta, AQD11898_P2.Wave_theta,...
    AQD11898_P2.Despike_vel_maj, AQD11898_P2.Despike_vel_min,...
    AQD11898_P2.freq_urms, AQD11898_P2.ubr, AQD11898_P2.ubj...
    AQD11898_P2.Suu_ADP, AQD11898_P2.Svv_ADP] = rotate_ADP_max_variance (...
    AQD11898_P2.VEL_XE_RAW, AQD11898_P2.VEL_YN_RAW, AQD11898_P2.Despike_vel_x, AQD11898_P2.Despike_vel_y,...
    AQD11898_P2.CORR_XE_RAW,AQD11898_P2.CORR_YN_RAW,60,1);

samples_fft = 1024;
AQD11898_P2.time_ft = AQD11898_P2.TIME_RAW(1:samples_fft:end); AQD11898_P2.time_ft(end) = [];

% Waves
lf=.04;                 % Hz - low frequency cutoff
maxfac=200;             % Maximum value of factor scaling pressure to waves
minspec=0.03;           % m^2/Hz - minimum spectral level for computing direction and spreading
Ndir=0;                 %deg - direction offset (includes compass error and misalignment of cable probe relative to case
                        % the offset for the Aquadopp Profiler is 0

parms=[lf maxfac minspec Ndir];
dt=1;
hp=0.9;
hv=0;
AQD11898_P2.z = 0.9-0.096:-0.03:0;

for i=1:round((size((AQD11898_P2.PRESSURE_RAW),1))/samples_fft)
    AQD11898_P2.Wave.meanDepth_ADP(i)=nanmean(AQD11898_P2.PRESSURE_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft))+.9;
    AQD11898_P2.Wave.time(i) = nanmean(AQD11898_P2.TIME_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft));

    aux2 = PcorFFTFun(AQD11898_P2.PRESSURE_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft),1,samples_fft,256,...
        AQD11898_P2.Wave.meanDepth_ADP(i),.9,1/25,1/2.5,1,'all','off','off');
    
    [AQD11898_P2.Wave.Hm0_ADP(i),AQD11898_P2.Wave.Tm01_ADP(i),AQD11898_P2.Wave.Tm02_ADP(i),AQD11898_P2.Wave.Tp_ADP(i),...
        AQD11898_P2.Wave.fp_ADP(i),AQD11898_P2.Wave.f_ADP(i,:),AQD11898_P2.Wave.Syy_ADP(i,:)]=WaveSpectraFun(aux2,...
        .9,samples_fft,256,AQD11898_P2.Wave.meanDepth_ADP(i)+.9,.9,1/25,1/2.5,1,-3,'off','off','off','off');
    
        clear aux aux2
        
        if sum (isnan(nanmean (AQD11898_P2.Despike_vel_x(i*samples_fft - (samples_fft - 1):i*samples_fft,4),2))) < 100
         
        auxu = AQD11898_P2.Despike_vel_x(i*samples_fft - (samples_fft - 1):i*samples_fft,4)';
        if isnan(auxu)==1; auxu(1)=0; end; if isnan(auxu(end))==1; auxu(end)=0; end;
        X = ~isnan(auxu);
        Y = cumsum(X-diff([1 X])/2); vu = interp1(1:nnz(X),auxu(X),Y); clear auxu X Y
        
        auxv = AQD11898_P2.Despike_vel_y(i*samples_fft - (samples_fft - 1):i*samples_fft,4)';
        if isnan(auxv)==1; auxv(1)=0; end; if isnan(auxv(end))==1; auxv(end)=0; end;
        X = ~isnan(auxv);
        Y = cumsum(X-diff([1 X])/2); vv = interp1(1:nnz(X),auxv(X),Y); clear auxv X Y
        
        auxp = AQD11898_P2.PRESSURE_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft,1)';
        if isnan(auxp)==1; auxp(1)=0; end; if isnan(auxp(end))==1; auxp(end)=0; end;
        X = ~isnan(auxp);
        Y = cumsum(X-diff([1 X])/2); vp = interp1(1:nnz(X),auxp(X),Y); clear auxv X Y
        
        [AQD11898_P2.Wave.Su(:,i),AQD11898_P2.Wave.Sp(:,i),AQD11898_P2.Wave.Dir(:,i),...
            AQD11898_P2.Wave.Spread(:,i),AQD11898_P2.Wave.F(:,i),AQD11898_P2.Wave.dF(:,i)] = ...
            wds(vu',vv',vp',dt,100,hp,hv,parms);
        else
%         [AQD11898_P2.Wave.ubr_sherwood(i),AQD11898_P2.Wave.Tbr_sherwood(i),AQD11898_P2.Wave.Snn_sherwood(:,i)]=ubspecdat(AQD11898_P2.Wave.meanDepth_ADP(i)+0.9,AQD11898_P2.Wave.Spp_ADP(i,:),AQD11898_P2.Wave.f_ADP(1,:));
        
        AQD11898_P2.Wave.Su(1:41,i)=nan;
        AQD11898_P2.Wave.Sp(1:41,i)=nan;
        AQD11898_P2.Wave.Dir(1:41,i)=nan;
        AQD11898_P2.Wave.Spread(1:41,i)=nan;
               
    end
        
end
   
AQD11898_P2.Wave.Dir = wrapTo360 (AQD11898_P2.Wave.Dir);
AQD11898_P2.Wave_Urms(1:5,:)=nan;AQD11898_P2.Wave_Urms(771:end,:)=nan;


for i=1:length(AQD11898_P2.Despike_vel_x)/300
    AQD11898_P2.Despike_vel_x_5min (i,:) = nanmean (AQD11898_P2.Despike_vel_x(i*300-299:i*300,:));
    AQD11898_P2.Despike_vel_y_5min (i,:) = nanmean (AQD11898_P2.Despike_vel_y(i*300-299:i*300,:));
AQD11898_P2.time_5min (i,:) = (AQD11898_P2.TIME_SENSOR_RAW(i*300-299));
end
AQD11898_P2.Current_mean_5min = sqrt ((AQD11898_P2.Despike_vel_x_5min.^2)+(AQD11898_P2.Despike_vel_y_5min.^2));

for i=1:length(AQD11898_P2.Despike_vel_x)/600
    AQD11898_P2.Despike_vel_x_10min (i,:) = nanmean (AQD11898_P2.Despike_vel_x(i*600-599:i*600,:));
    AQD11898_P2.Despike_vel_y_10min (i,:) = nanmean (AQD11898_P2.Despike_vel_y(i*600-599:i*600,:));
AQD11898_P2.time_10min (i,:) = (AQD11898_P2.TIME_SENSOR_RAW(i*600-599:i*600));
end
AQD11898_P2.time_10min(:,2:end)=[];
AQD11898_P2.Current_mean_10min = sqrt ((AQD11898_P2.Despike_vel_x_10min.^2)+(AQD11898_P2.Despike_vel_y_10min.^2));

for i=1:length(AQD11898_P2.Despike_vel_x)/900
    AQD11898_P2.Despike_vel_x_15min (i,:) = nanmean (AQD11898_P2.Despike_vel_x(i*900-899:i*900,:));
    AQD11898_P2.Despike_vel_y_15min (i,:) = nanmean (AQD11898_P2.Despike_vel_y(i*900-899:i*900,:));
AQD11898_P2.time_15min (i,:) = (AQD11898_P2.TIME_SENSOR_RAW(i*900-899));
end
AQD11898_P2.Current_mean_15min = sqrt ((AQD11898_P2.Despike_vel_x_15min.^2)+(AQD11898_P2.Despike_vel_y_15min.^2));

AQD11898_P2.Current_Mean = sqrt(AQD11898_P2.Despike_vel_maj.^2+AQD11898_P2.Despike_vel_min.^2);
%% Read FLNTU
caminho_flntu = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\2.Part_2\Quadpod\FLNTU'
file = fileread([caminho_flntu '\3194_3.raw']);

lines = regexp (file, '\n', 'split');

[mmddyy,rest] = strtok(lines,char(9));
%Creates vector time
formatIn='m/dd/yy';
for i = 2:length(lines)-1
    try
        OBSday(i-1,:)=datevec(mmddyy(i),formatIn);
    catch
        OBSday(i-1,1:6)=nan;
    end
end

[hhmm,rest2] = strtok(rest,char(9));
formatIn2='HH:MM:SS';
for i = 2:length(lines)-1
    try
        OBShour(i-1,:)=datevec(hhmm(i),formatIn2);
    catch
        OBShour(i-1,:)=nan;
    end
end

FLNTU.time = struct;
FLNTU.time=datenum([OBSday(:,1:3) OBShour(:,4:6)]);

fid = fopen([caminho_flntu '\3194_3.raw'],'r');
var_dir=textscan(fid, '%s %s %d %d %d %d %d', 'Headerlines',1);
FLNTU.raw = struct;
FLNTU.raw = var_dir{1,6};
FLNTU.raw(1)=[];

window = 10;flntu.final = [];
for i=1:size(FLNTU.raw,1)-window
    FLNTU.aux(i,2) =nanmin(FLNTU.raw(i:i+(window-1),1));
end
FLNTU.final = FLNTU.aux(:,2);
FLNTU.NTU = 0.0242 * (double(FLNTU.final) - 50);
FLNTU.NTU(FLNTU.NTU<0)=0;
%% Read Aquascat
cam_Aquascat = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\2.Part_2\Aquascat';
pp=read_aquatec_v4(fullfile(cam_Aquascat,'20190320220000.aqa'))

%% Read Vector
cam_vector = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\2.Part_2\Quadpod\VEC1670';
cd (cam_vector)    
% Get_Nortek_ADV2NetCDF (cam_vector,'ADVQU201','ADV1670_P2_V0','',0,0,.1,.2,0,60,58,8)
ADV1670_P2 = ncreadall (fullfile(cam_vector,'\ADV1670_V0'))
  
[ADV1670_P2.Despike_vel_x(1,:), ADV1670_P2.Despike_vel_y(1,:), ADV1670_P2.Despike_vel_z(1,:)] = func_despike_phasespace3d_3var (ADV1670_P2.VEL_X(1,:),...
    ADV1670_P2.VEL_Y(1,:), ADV1670_P2.VEL_Z(1,:), 0);

for i=1:length(ADV1670_P2.VEL_X)/(8192)
    ADV1670_P2.Despike_vel_x_15min (i) = nanmean(ADV1670_P2.Despike_vel_x(i*8192-8191:i*8192));
    ADV1670_P2.Despike_vel_y_15min (i) = nanmean(ADV1670_P2.Despike_vel_y(i*8192-8191:i*8192));
    ADV1670_P2.Despike_vel_z_15min (i) = nanmean(ADV1670_P2.Despike_vel_z(i*8192-8191:i*8192));
end
ADV1670_P2.Despike_VEL_15min = sqrt((ADV1670_P2.Despike_vel_x_15min.^2)+(ADV1670_P2.Despike_vel_y_15min.^2));
ADV1670_P2.time_15_min = ADV1670_P2.TIME(1:8*60*15:end);ADV1670_P2.time_15_min(end)=[];

[ADV1670_P2.vel_maj, ADV1670_P2.vel_min , ADV1670_P2.Current_Mean,...
    ADV1670_P2.Current_theta, ADV1670_P2.Wave_theta, ...
    ADV1670_P2.Despike_cur_maj, ADV1670_P2.Despike_cur_min,...
    ADV1670_P2.fx, ADV1670_P2.ubr, ADV1670_P2.ubr_spect,...
    ADV1670_P2.Suu, ADV1670_P2.Svv] = rotate_continuous_ADV_max_variance (ADV1670_P2.VEL_X,ADV1670_P2.VEL_Y,...
    ADV1670_P2.Despike_vel_x,ADV1670_P2.Despike_vel_y,8);


ADV1670_P2.time_30s = ADV1670_P2.TIME(1:240:end);ADV1670_P2.time_30s(end)=[];
for i = 1:round(size(ADV1670_P2.PRES,2)/8192)-1
    if i==619 || i==620
        
        ADV1670_P2.Wave.Hm0(i)=nan;
        ADV1670_P2.Wave.Tm01(i)=nan;
        ADV1670_P2.Wave.Tm02(i)=nan;
        ADV1670_P2.Wave.Tp(i)=nan;
        ADV1670_P2.Wave.fp(i)=nan;
        ADV1670_P2.Wave.f(i,:)=nan;
        ADV1670_P2.Wave.SA_1(i,:)=nan;
    else
       
        ADV1670_P2.Wave.MeanDepth(i)=nanmean(ADV1670_P2.PRES(i*8192-8191:i*8192));
        aux=PcorFFTFun(ADV1670_P2.PRES(i*8192-8191:i*8192)',8,1024,2^10,ADV1670_P2.Wave.MeanDepth(i)+.4,0.4,1/25,1/2.5,1,'all','off','off');
        
        [ADV1670_P2.Wave.Hm0(i),ADV1670_P2.Wave.Tm01(i),ADV1670_P2.Wave.Tm02(i),ADV1670_P2.Wave.Tp(i),ADV1670_P2.Wave.fp(i),ADV1670_P2.Wave.f(i,:),...
            ADV1670_P2.Wave.SA_1(i,:)]=WaveSpectraFun(aux,8,1024,2^10,aux+.4,.4,.1,.8,1,-3,'off','off','off','off');
       
        clear aux
    end
end

ADV1670_P2.time_15_min=ADV1670_P2.TIME(1:8192:end);ADV1670_P2.time_15_min(end)=[];
%%
clearvars -except ADV1670_P2 AQD11898_P2 FLNTU

save ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\2.Part_2\Quadpod\Quadpod_instruments.mat',...  
    'ADV1670_P2','AQD11898_P2','FLNTU','-v7.3')
    
