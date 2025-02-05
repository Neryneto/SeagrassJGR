%% Read, filter, calculate Urms values and save Tripod instruments

% Tripod:
%     Aquadopp 11898 (8 to 15/03)

%% Tripod instruments Aquadopp AQD 11898 (Tripod near Quadpod)

caminho_AQD = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898';

AQD11898 = ncreadall (fullfile(caminho_AQD,'\AQD11898_V0'));

for i = 1:size(AQD11898.VEL_XE_RAW, 2)
    [AQD11898.Despike_vel_x(:,i), AQD11898.Despike_vel_y(:,i), AQD11898.Despike_vel_z(:,i)] = func_despike_phasespace3d_3var_ADP (AQD11898.VEL_XE_RAW(:,i),...
        AQD11898.VEL_YN_RAW(:,i), AQD11898.VEL_ZU_RAW(:,i), 0);
end
AQD11898.z = 0.1+0.104+[1:37].*0.05+0.1

[AQD11898.vel_maj, AQD11898.vel_min, AQD11898.Mean_Current...
    AQD11898.Current_theta, AQD11898.Wave_theta,...
    AQD11898.Despike_vel_maj, AQD11898.Despike_vel_min,...
    AQD11898.freq_urms, AQD11898.ubr, AQD11898.ubj...
    AQD11898.Suu_ADP, AQD11898.Svv_ADP] = rotate_ADP_max_variance (...
    AQD11898.VEL_XE_RAW, AQD11898.VEL_YN_RAW, AQD11898.Despike_vel_x, AQD11898.Despike_vel_y,...
    AQD11898.CORR_XE_RAW,AQD11898.CORR_YN_RAW,60,1);

samples_fft = 1024;
AQD11898.Wave.time = AQD11898.TIME_RAW(1:samples_fft:end); AQD11898.Wave.time(end) = [];
AQD11898.time_ft = AQD11898.TIME_RAW(1:samples_fft:end); AQD11898.time_ft(end) = [];

AQD11898.freestream = nanmean (AQD11898.Mean_Current(:,30:35),2);
BURSTS = floor(size(AQD11898.Despike_vel_x,1)/1024);
  
AQD11898.Despike_vel = sqrt((AQD11898.Despike_vel_x.^2)+(AQD11898.Despike_vel_y.^2));
heights = AQD11898.z(30:35);

% Waves
lf=.04;                 % Hz - low frequency cutoff
maxfac=200;             % Maximum value of factor scaling pressure to waves
minspec=0.03;           % m^2/Hz - minimum spectral level for computing direction and spreading
Ndir=0;                 %deg - direction offset (includes compass error and misalignment of cable probe relative to case
                        % the offset for the Aquadopp Profiler is 0

parms=[lf maxfac minspec Ndir];
dt=1;
hp=0.14;
hv=1.81-.14;

for i=1:round((size((AQD11898.PRESSURE_RAW),1))/samples_fft)-1
    try
        AQD11898.meanDepth_15min(i)=nanmean(AQD11898.PRESSURE_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft));
        AQD11898.time_15min(i) = nanmean(AQD11898.TIME_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft));

        
        aux2 = PcorFFTFun(AQD11898.PRESSURE_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft),1,samples_fft,256,...
            AQD11898.meanDepth_15min(i)+.15,.15,1/25,1/2.5,1,'all','off','off');
        
        [AQD11898.Wave.Hm0_ADP(i),...
            AQD11898.Wave.Tm01_ADP(i), AQD11898.Wave.Tm02_ADP(i), AQD11898.Wave.Tp_ADP(i),...
            AQD11898.Wave.fp_ADP(i),AQD11898.Wave.f_ADP(i,:), AQD11898.Wave.Spp_ADP(i,:)] = WaveSpectraFun(aux2,...
            1,samples_fft,256,AQD11898.meanDepth_15min(i)+.15,...
            .15,1/25,1/2.5,1,-3,'off','off','off','off');
        
        clear aux aux2
        
        if sum (isnan(nanmean (AQD11898.Despike_vel_x(i*samples_fft - (samples_fft - 1):i*samples_fft,30:35),2))) < 100
         
        auxu = nanmean (AQD11898.Despike_vel_x(i*samples_fft - (samples_fft - 1):i*samples_fft,30:35),2)';
        if isnan(auxu)==1; auxu(1)=0; end; if isnan(auxu(end))==1; auxu(end)=0; end;
        X = ~isnan(auxu);
        Y = cumsum(X-diff([1 X])/2); vu = interp1(1:nnz(X),auxu(X),Y); clear auxu X Y
        
        auxv = nanmean (AQD11898.Despike_vel_y(i*samples_fft - (samples_fft - 1):i*samples_fft,30:35),2)';
        if isnan(auxv)==1; auxv(1)=0; end; if isnan(auxv(end))==1; auxv(end)=0; end;
        X = ~isnan(auxv);
        Y = cumsum(X-diff([1 X])/2); vv = interp1(1:nnz(X),auxv(X),Y); clear auxv X Y
        
        auxp = AQD11898.PRESSURE_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft,1)';
        if isnan(auxp)==1; auxp(1)=0; end; if isnan(auxp(end))==1; auxp(end)=0; end;
        X = ~isnan(auxp);
        Y = cumsum(X-diff([1 X])/2); vp = interp1(1:nnz(X),auxp(X),Y); clear auxv X Y
        
        [AQD11898.Wave.Su(i,:),AQD11898.Wave.Sp(i,:),AQD11898.Wave.Dir(i,:),AQD11898.Wave.Spread,AQD11898.Wave.F,AQD11898.Wave.dF] = wds(vu',vv',vp',dt,100,hp,hv,parms);
        end
        [AQD11898.Wave.ubr_sherwood(i),AQD11898.Wave.Tbr_sherwood(i),AQD11898.Wave.Snn_sherwood(:,i)]=ubspecdat(AQD11898.meanDepth_15min(i)+0.15,AQD11898.Wave.Spp_ADP(i,:),AQD11898.Wave.f_ADP(1,:));
        catch
        AQD11898.Wave.Su(i,:)=nan;
        AQD11898.Wave.Sp(i,:)=nan;
        AQD11898.Wave.Dir(i,:)=nan;
        AQD11898.Wave.Spread(i,:)=nan;
               
    end
        
end
   
AQD11898.Wave.Dir = wrapTo360 (AQD11898.Wave.Dir);
clearvars -except AQD11898

save ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments','AQD11898' ,'-v7.3')