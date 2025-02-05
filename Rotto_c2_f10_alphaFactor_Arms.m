load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments')
% load ('C:\Users\22371812\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\2.Part_2\Quadpod\Quadpod_instruments.mat');
% load ('C:\Users\22371812\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\Sawhorse_instruments.mat')
% load('C:\Users\22371812\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd.mat')
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Quadpod_instruments.mat');


%% Test different alphar
% plot(AQD11898_P2.Wave.time,AQD11898_P2.ubr(1:end-1,27))
% % plot(AQD11898_P2.Wave.time,AQD11898_P2.Mean_Current(1:end-1,1))
% % plot(AQD11898_P2.Wave.time,AQD11898_P2.Mean_Current(1:end-1,27))
% plot(AQD11898_P2.time_15min,AQD11898_P2.Mean_Current(:,1))
% plot(AQD11898_P2.time_15min,AQD11898.Mean_Current(:,20))
% % plot(AQD11898.Wave.time,AQD11898.Wave.Tm01_ADP./10)
% % plot(AQD11898_P2.Wave.time,AQD11898_P2.Wave.Tm01_ADP/10)
% % plot(AQD11898_P2.Wave.time,AQD11898_P2.Wave.Tp_ADP/10)
% % plot(AQD11898.Wave.time,AQD11898.Wave.Tp_ADP./10)
% plot(RBR.ST9.time_wave,RBR.ST9.ubr_LWT)
% hold on
% plot(AQD11898.time_15min,AQD11898.ubr(:,1))
% plot(ADV1670.time_ft,ADV1670.ubr)
% for i=1:122
% plot(Vectrino(i).InitialTime_Matlab,nanmean(Vectrino(i).ubr),'.k')
% plot(Vectrino(i).InitialTime_Matlab,violin_vectrino_vector_urms(i),'.g')
% end
% %% Calculate \alpha^R
% for idx=1:length(AQD11898_P2.Wave.meanDepth_ADP)
% [AQD11898_P2.Wave.Lp(idx),AQD11898_P2.Wave.kp(idx),AQD11898_P2.Wave.sigmatm02(idx)]=disper(AQD11898_P2.Wave.meanDepth_ADP(idx)+0.9,AQD11898_P2.Wave.Tm02_ADP(idx));
% end
% AQD11898_P2.Wave.Arms=AQD11898_P2.Urms(:,1).*AQD11898_P2.Wave.Tm02_ADP./(2*pi);
% 
% phij = (cosh (AQD11898_P2.Wave.kp.*AQD11898_P2.z(5)))./(cosh (AQD11898_P2.Wave.kp.*AQD11898_P2.z(1)));
% phij(phij<1 | phij>1.025)=nan;
% alphar = sqrt((phij'.^2).*(AQD11898_P2.ubr(:,1))./AQD11898_P2.ubr(:,5));
% 
% alphaR=AQD11898_P2.ubr(:,1)./AQD11898_P2.ubr(:,5);
% Arms_RBR_tm02 = RBR.ST9.ubr_LWT.*RBR.ST9.Tm02./(2*pi);
% Arms_RBR_tm01 = RBR.ST9.ubr_LWT.*RBR.ST9.Tm01./(2*pi);
% Arms_RBR_tp = RBR.ST9.ubr_LWT.*RBR.ST9.Tp./(2*pi); Arms_RBR_tp(Arms_RBR_tp>5)=nan;
% AQD11898_P2.ubr(1,:)=[];
% Arms_aquadopp_p2_tm01 = AQD11898_P2.ubr(:,27).*AQD11898_P2.Wave.Tm01_ADP'./(2*pi); Arms_aquadopp_p2_tm01(Arms_aquadopp_p2_tm01>0.5)=nan;
% Arms_aquadopp_p2_tm02 = AQD11898_P2.ubr(:,27).*AQD11898_P2.Wave.Tm02_ADP'./(2*pi);Arms_aquadopp_p2_tm02(Arms_aquadopp_p2_tm02>0.5)=nan;
% Arms_aquadopp_p2_tp = AQD11898_P2.ubr(:,27).*AQD11898_P2.Wave.Tp_ADP'./(2*pi);

% Arms_aquadopp_tm01 = AQD11898.ubr(:,27).*AQD11898.Wave.Tm01_ADP'./(2*pi); Arms_aquadopp_tm01(Arms_aquadopp_tm01>0.5)=nan;
% Arms_aquadopp_tm02 = AQD11898.ubr(:,27).*AQD11898.Wave.Tm02_ADP'./(2*pi);Arms_aquadopp_tm02(Arms_aquadopp_tm02>0.5)=nan;
% Arms_aquadopp_tp = AQD11898.ubr(:,27).*AQD11898.Wave.Tp_ADP'./(2*pi);

% plot(RBR.ST9.time_wave, Arms_RBR_tm01)
% hold on
% plot(RBR.ST9.time_wave, Arms_RBR_tm02)
% plot(RBR.ST9.time_wave, Arms_RBR_tp)
% plot(AQD11898_P2.Wave.time, Arms_aquadopp_p2_tm01)
% % plot(AQD11898_P2.Wave.time, Arms_aquadopp_tm02)
% plot(AQD11898_P2.Wave.time, Arms_aquadopp_p2_tp)
% plot(AQD11898_P2.time_ft, AQD11898_P2.ubr(:,27))
% hold on
% plot(ADV1670_P2.time_15_min,ADV1670_P2.ubr)
% plot(ADV1670.time_ft,ADV1670.ubr)
% plot(AQD11898.time_15min, AQD11898.ubr(:,1))

%%
avg_vectrino_rms = [];

for ii=1:length(Vectrino);
    avg_vectrino_rms = nanmean(Vectrino(ii).ubr(11:15));

    [~,idx] = min(abs(AQD11898.time_15min-Vectrino(ii).InitialTime_Matlab));
    [~,idx_vector] = min(abs(ADV1670.time_ft-Vectrino(ii).InitialTime_Matlab));

%     if Vectrino(ii).BottomDistance<0.4 && Vectrino(ii).InitialTime_Matlab>737494.5
    if Vectrino(ii).BottomDistance<0.1 && Vectrino(ii).InitialTime_Matlab>737494.5
        violin_vectrino_vector_urms (ii) = avg_vectrino_rms/ADV1670.ubr(idx_vector);
        % violin_vectrino_aqd_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,1));
        violin_vectrino_aqd_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,1));
        
        Arms_aquadopp_tm01 (ii) = AQD11898.ubr(idx,1).*AQD11898.Wave.Tm01_ADP(idx)'./(2*pi); 
        Arms_aquadopp_tm02 (ii) = AQD11898.ubr(idx,1).*AQD11898.Wave.Tm02_ADP(idx)'./(2*pi);
        Arms_aquadopp_tp(ii) = AQD11898.ubr(idx,1).*AQD11898.Wave.Tp_ADP(idx)'./(2*pi);

    else
        violin_vectrino_vector_urms (ii) = nan;
        violin_vectrino_aqd_urms (ii) = nan;
        
         Arms_aquadopp_tm01 (ii) = nan;
        Arms_aquadopp_tm02 (ii) = nan;
        Arms_aquadopp_tp(ii) = nan;
    end
end
%%
% plot (violin_vectrino_aqd_Arms,violin_vectrino_aqd_urms,'.','markersize',10)
% plot (violin_vectrino_vector_Arms,violin_vectrino_vector_urms,'.','markersize',10)
f10=figure(10)

a=Arms_aquadopp_tm02./(1./sqrt(512));
violin_vectrino_aqd_urms(violin_vectrino_aqd_urms<0.41)=nan;
violin_vectrino_aqd_urms(61)=0.3222;
valid_idx = ~isnan(a) & ~isnan(violin_vectrino_aqd_urms);

% Keep only valid data (non-NaN rows)
a_valid = a(valid_idx);
violin_vectrino_aqd_urms_valid = violin_vectrino_aqd_urms(valid_idx);
% [h, p] = ttest(a_valid', violin_vectrino_aqd_urms_valid');

plot(a, violin_vectrino_aqd_urms,'.b','markersize',16)
hold on
plot([min(a):0.1:max(a)+1],[min(a):0.1:max(a)+1].*(-0.0968)+0.8246,'r');

% plot(violin_vectrino_aqd_Arms,0.8289.*violin_vectrino_aqd_Arms+0.04999,'r')
xlabel('$A_{\infty,R}/N^{-1/2}$','Interpreter','Latex','fontsize',14)
ylabel('$\alpha_{w,R}$ [ ]','Interpreter','Latex','fontsize',14)
grid on
box on

xlim([1 4.5])
ylim([0 0.8])
% [h, p] = ttest(a, violin_vectrino_aqd_urms);

text(1.016,0.23,'$R^2=0.52, p=0.003, \nu=14$','Interpreter','Latex','fontsize',14)
text(1.016,0.13,'$\alpha_{w,R} = -0.097*A_{\infty,R}/N^{(-1/2)}+0.82$',...
    'Interpreter','Latex','fontsize',14)



%%
print(f10,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\Figure10_alpha_Arms.tiff','-dtiff','-r300')
