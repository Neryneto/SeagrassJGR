%% Read, filter, calculate Urms values and save Quadpod instruments

% Quadpod:
%     ADV 1670 (08/03 12:00 to 15/03 8:45)
%     Vectrino (09/03 14:40 to 13/03 15:59)
%     RBR

%% Read Aquadopp in Tripod to use wave angle
load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments.mat')

%% Vector ADV 1670

caminho_VEC = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\VEC1670';
% Get_Nortek_ADV2NetCDF ([caminho_VEC],'QD_ADV01','ADV1670_V0','',0,0,.1,.2,0,60,58,8)
ADV1670 = ncreadall (fullfile(caminho_VEC,'\ADV1670_V0'))

for i=1:floor(length(ADV1670.VEL_X)./11717)
    [ADV1670.Despike_vel_x(1,i*11717-11716:i*11717), ADV1670.Despike_vel_y(1,i*11717-11716:i*11717), ADV1670.Despike_vel_z(1,i*11717-11716:i*11717),...
    ] = func_despike_phasespace3d_3var (ADV1670.VEL_X(1,i*11717-11716:i*11717),...
    ADV1670.VEL_Y(1,i*11717-11716:i*11717), ADV1670.VEL_Z(1,i*11717-11716:i*11717), 0);
end

% Get Wave Urms velocities and angle, current angle and current mean
[ADV1670.vel_maj, ADV1670.vel_min , ADV1670.Current_Mean,...
    ADV1670.Current_theta, ADV1670.Wave_theta, ...
    ADV1670.Despike_cur_maj, ADV1670.Despike_cur_min,...
    ADV1670.fx, ADV1670.ubr, ADV1670.ubr_spect,...
    ADV1670.Suu, ADV1670.Svv] = rotate_continuous_ADV_max_variance (ADV1670.VEL_X,ADV1670.VEL_Y,...
    ADV1670.Despike_vel_x,ADV1670.Despike_vel_y,8);

ADV1670.time_ft = ADV1670.TIME_SENSOR(1:1024:end); ADV1670.time_ft(end) = [];

% Read Waves
for i = 1:round(size(ADV1670.PRES,2)/8192)
    if i>=544
        ADV1670.Wave.Hm0(i,1)=nan;
        ADV1670.Wave.Tm01(i,1)=nan;
        ADV1670.Wave.Tm02(i,1)=nan;
        ADV1670.Wave.Tp(i,1)=nan;
        ADV1670.Wave.fp(i,1)=nan;
        ADV1670.Wave.f(i,:)=nan;
        ADV1670.Wave.SA_1(i,:)=nan;
    else
        ADV1670.meanDepth_15min (i) = nanmean(ADV1670.PRES(i*8192-8191:i*8192));
        aux=PcorFFTFun(ADV1670.PRES(i*8192-8191:i*8192)',8,1024,1024,...
            ADV1670.meanDepth_15min (i)+.4,0.4,1/25,1/2.5,1,'all','off','off');
                
        [ADV1670.Wave.Hm0(i,1),ADV1670.Wave.Tm01(i,1),ADV1670.Wave.Tm02(i,1),ADV1670.Wave.Tp(i,1),ADV1670.Wave.fw(i,1),ADV1670.Wave.f(i,:),...
            ADV1670.Wave.SA_1(i,:)]=WaveSpectraFun(aux,8,1024,512,ADV1670.meanDepth_15min (i)+.4,.4,.1,.8,1,-3,'off','off','off','off');
        clear aux
    end
end

if length (ADV1670.time_ft) > length(ADV1670.Wave.Hm0)
    ADV1670.time_ft(end)= [];
end

%% Calculate Reynolds Stresses
for i=1:length(ADV1670.Current_Mean)
    if sum(isnan(ADV1670.vel_maj(i*8192-8191:i*8192)))>=8192 ||...
            sum(isnan(ADV1670.vel_min(i*8192-8191:i*8192)))>=8192/2
        ADV1670.Benilov.upup_bar(i)=nan;
        ADV1670.Benilov.vpvp_bar(i)=nan;
        ADV1670.Benilov.wpwp_bar(i)=nan;
        ADV1670.Benilov.upwp_bar(i)=nan;
        ADV1670.Benilov.uwave_wwave_bar(i) = nan;
        
        ADV1670.BrickerMonismith.upwp(i)=nan;
        ADV1670.BrickerMonismith.upup_bar (i)=nan;
        ADV1670.BrickerMonismith.vpvp_bar (i)=nan;
        ADV1670.BrickerMonismith.wpwp_bar (i)=nan;
        ADV1670.BrickerMonismith.upwp_bar (i)=nan;
        ADV1670.BrickerMonismith.uwave_wwave_bar (i)=nan;
        
        ADV1670.YoungWebster.uu(i)=nan;
        ADV1670.YoungWebster.vv(i)=nan;
        ADV1670.YoungWebster.ww(i)=nan;
        ADV1670.YoungWebster.uw(i)=nan;
        ADV1670.YoungWebster.uw_wave(i)=nan;
        
    else
        
        % ADV Reynolds Stress Tensor from Benilov
        [ADV1670.Benilov.upup_bar(i), ADV1670.Benilov.vpvp_bar(i), ADV1670.Benilov.wpwp_bar(i),...
            ADV1670.Benilov.upwp_bar(i), ADV1670.Benilov.uwave_wwave_bar(i), ADV1670.Benilov.Suu(:,i),...
            ADV1670.Benilov.S_uwave_uwave(:,i), ADV1670.Benilov.Supup(:,i)] = ...
            calculate_Benilov (ADV1670.vel_maj(i*8192-8191:i*8192),...
            ADV1670.vel_min(i*8192-8191:i*8192),...
            ADV1670.Despike_vel_z(i*8192-8191:i*8192)',...
            ADV1670.PRES(i*8192-8191:i*8192)',...
            0.4,...
            8,...
            1,...
            1026);
        
        % ADV Reynolds Stress Tensor from Bricker and Monismith
        [ADV1670.BrickerMonismith.upup_bar(i), ADV1670.BrickerMonismith.vpvp_bar(i), ADV1670.BrickerMonismith.wpwp_bar(i),...
            ADV1670.BrickerMonismith.upwp_bar(i), ADV1670.BrickerMonismith.uwave_wwave_bar(i)] = ...
            calculate_BrickerMonismith (ADV1670.vel_maj(i*8192-8191:i*8192),...
            ADV1670.vel_min(i*8192-8191:i*8192),...
            ADV1670.Despike_vel_z(i*8192-8191:i*8192)',...
            ADV1670.PRES(i*8192-8191:i*8192)',...
            8);
        
        % ADV Reynolds Stress Tensor from Young and Webster
        [ADV1670.YoungWebster.uu(i), ADV1670.YoungWebster.vv(i), ADV1670.YoungWebster.ww(i),...
            ADV1670.YoungWebster.uw(i), ADV1670.YoungWebster.uw_wave(i)] = calculate_youngWebster (...
            ADV1670.vel_maj(i*8192-8191:i*8192),...
            ADV1670.vel_min(i*8192-8191:i*8192),...
            ADV1670.Despike_vel_z(i*8192-8191:i*8192)',...
            ADV1670.PRES(i*8192-8191:i*8192),...
            8);
    end
    
end
ADV1670.Benilov.TKE = 0.5.*(ADV1670.Benilov.upup_bar.^2 + ADV1670.Benilov.vpvp_bar.^2 + ADV1670.Benilov.wpwp_bar.^2);
ADV1670.BrickerMonismith.TKE = 0.5.*(ADV1670.BrickerMonismith.upup_bar.^2 + ADV1670.BrickerMonismith.vpvp_bar.^2 + ADV1670.BrickerMonismith.wpwp_bar.^2);
ADV1670.YoungWebster.TKE = 0.5.*(ADV1670.YoungWebster.uu.^2 + ADV1670.YoungWebster.vv.^2 + ADV1670.YoungWebster.ww.^2);
%% RBR_MATT
caminho_RBR = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\054322_20190315_1106'
fileID = fopen(fullfile(caminho_RBR,'054322_20190315_1106_data.txt'));
arq_txt=textscan(fileID, '%s %f','HeaderLines',1,'Delimiter',',');

formatIn='yyyy-mm-dd HH:MM:SS';
date_RBR=datenum(arq_txt{1,1},formatIn);

ntu_raw = arq_txt{1,2};
ntu_raw(ntu_raw>8 | ntu_raw<0)=nan;
for i=2:length(ntu_raw)
    if abs(ntu_raw(i)-ntu_raw(i-1))>5 && abs(ntu_raw(i)-ntu_raw(i+1))>5
        ntu_raw(i) = nan;
    end
end
ntu=movmean(ntu_raw,12,'omitnan');
RBR_Matt = [date_RBR ntu];

%% VECTRINO
caminho_Vectrino = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Vectrino_Rottnest'
dirlist = dir(fullfile(caminho_Vectrino,'\**\*.mat'))

for K = 1: size(dirlist,1)
    tmp = load(fullfile(dirlist(K).folder, dirlist(K).name));
    Vectrinos(K).Data = tmp.Data;
    Vectrinos(K).Config = tmp.Config;
    
    for i = 1:size(Vectrinos(K).Data.Profiles_VelX, 2)
        [Vectrino(K).Despike_vel_x(:,i), Vectrino(K).Despike_vel_y(:,i), Vectrino(K).Despike_vel_z(:,...
            i)] = func_despike_phasespace3d_3var_ADP (Vectrinos(K).Data.Profiles_VelX(:,i),...
            Vectrinos(K).Data.Profiles_VelY(:,i), Vectrinos(K).Data.Profiles_VelZ1(:,i), 0);
        
        Vectrino(K).Time = Vectrinos(K).Data.Profiles_HostTimeMatlab' + 8/24;
    end
        Vectrino(K).InitialTime = tmp.Data.Profiles_hostStartDate;
    Vectrino(K).InitialTime_Matlab = tmp.Data.Profiles_hostFirstRecordMatlab + 8/24;

    for j=1:31
        [l c] = min(abs(AQD11898.time_ft-Vectrino(K).InitialTime_Matlab));
        w = complex(Vectrino(K).Despike_vel_x(:,j),Vectrino(K).Despike_vel_y(:,j));
        w (w~=isfinite(w)==0);
        wr = w.*exp(-1i.*AQD11898.Wave_theta(c)*pi/180);
        Vectrino(K).vel_maj (:,j) = real(wr);
        Vectrino(K).vel_min (:,j) = imag(wr);
        
%         Vectrino(K).Uwrms(1,j) = UcUrms(Vectrino(K).vel_maj (:,j));
        Vectrino(K).ubr (1,j) = sqrt(2*(nanvar(Vectrino(K).vel_maj (:,j)-...
        nanmean(Vectrino(K).vel_maj (:,j))) + nanvar(Vectrino(K).vel_min (:,j)-...
        nanmean(Vectrino(K).vel_min (:,j)))));
        
    end
    xtmp = nanmean(Vectrino(K).vel_maj(:,10:15),2);
    L = length(xtmp);
    aux = xtmp; if isnan(xtmp(1))==1; aux(1)=0; end; if isnan(xtmp(end))==1; aux(end)=0; end;
    X = ~isnan(aux);
    Y = cumsum(X-diff([1; X])/2); eta = interp1(1:nnz(X),aux(X),Y);
    [Vectrino(K).Sxx,Vectrino(K).fx]=calculatespectrum(25,eta,2.^nextpow2(length(eta))/2); clear L aux X Y eta
    
    deltaf = Vectrino(K).fx(3)-Vectrino(K).fx(2);
    
    ytmp = nanmean(Vectrino(K).vel_min(:,10:15),2);
    L = length(ytmp);
    aux = ytmp; if isnan(ytmp(1))==1; aux(1)=0; end; if isnan(ytmp(end))==1; aux(end)=0; end;
    X = ~isnan(aux);
    Y = cumsum(X-diff([1; X])/2); eta = interp1(1:nnz(X),aux(X),Y);
    [Vectrino(K).Syy,~]=calculatespectrum(25,eta,2.^nextpow2(length(eta))/2); clear L aux X Y eta
    
    Vectrino(K).ubj = sqrt(2*(Vectrino(K).Sxx  + Vectrino(K).Syy).*deltaf);
    Vectrino(K).folder = dirlist(K).folder;
    Vectrino(K).name = dirlist(K).name;
    Vectrino(K).BottomDistance = mode(round(tmp.Data.BottomCheck_BottomDistance*100)/100);
    Vectrino(K).Temperature = tmp.Data.Profiles_Temperature;
    Vectrino(K).SNR1 = tmp.Data.Profiles_SNRBeam1;
    Vectrino(K).SNR2 = tmp.Data.Profiles_SNRBeam2;
    Vectrino(K).SNR3 = tmp.Data.Profiles_SNRBeam3;
    Vectrino(K).SNR4 = tmp.Data.Profiles_SNRBeam4;
    Vectrino(K).Corr1 = tmp.Data.Profiles_CorBeam1;
    Vectrino(K).Corr2 = tmp.Data.Profiles_CorBeam2;
    Vectrino(K).Corr3 = tmp.Data.Profiles_CorBeam3;
    Vectrino(K).Corr4 = tmp.Data.Profiles_CorBeam4;
    Vectrino(K).Velz2 = tmp.Data.Profiles_VelZ2;
               
    clear tmp l c
end

T = struct2table(Vectrino); % convert the struct array to a table
sortedT = sortrows(T, 'InitialTime_Matlab'); % sort the table by 'DOB'
clear Vectrino
Vectrino = table2struct(sortedT);

%cleaning SNR<15 or Correlation<75
for K = 1:length(Vectrino)
    for i=1:numel(Vectrino(K).SNR1)
        if Vectrino(K).SNR1(i)<15 || Vectrino(K).Corr1(i)<75
            Vectrino(K).Despike_vel_x(i) =nan;
        end
        if Vectrino(K).SNR2(i)<15 || Vectrino(K).Corr2(i)<75
            Vectrino(K).Despike_vel_y(i) =nan;
        end
        if Vectrino(K).SNR3(i)<15 || Vectrino(K).Corr3(i)<75
            Vectrino(K).Despike_vel_z(i) =nan;
        end
    end
Vectrino(K).avgx=nanmean(Vectrino(K).Despike_vel_x);    
Vectrino(K).avgy=nanmean(Vectrino(K).Despike_vel_y);
Vectrino(K).MeanVel=sqrt((Vectrino(K).avgx.^2)+(Vectrino(K).avgy.^2));    
Vectrino(K).Despike_Vel=sqrt((Vectrino(K).Despike_vel_x.^2)+(Vectrino(K).Despike_vel_y.^2));  
end

% clearvars -except ADV1670 RBR_Matt Vectrino Wave_Quadpod
save ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Quadpod_instruments.mat',...
    'ADV1670','Vectrino','RBR_Matt','-v7.3')
