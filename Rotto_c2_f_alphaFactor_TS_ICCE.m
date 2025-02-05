% BrickerMonismith = load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\BrickerMonismith.mat');
% BrickerMonismith_vectrino = load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\BrickerMonismith_vectrino.mat');
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Quadpod_instruments.mat');
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments')
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_lev.mat')
load ('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\Sawhorse_instruments.mat')
%% Calculate spectra

Vec_h=[];Svel = [];

heights = 0.05:0.01:max(Vec_h)+0.01;heights=heights';

for num=1:122 %each profile
    Vec_h(num,1) = nanmean(Vectrino(num).BottomDistance);
end
for i = 1:floor(length(ADV1670.PRES)/8192)
    ADV1670.press_15min (i) = nanmean(ADV1670.PRES(i*8192-8191:i*8192));
end
%% Calculate ratio Vectrino/Free stream and Aquadopp/free stream
avg_vectrino = [];
avg_vectrino_rms = [];
ratio_aqd_freestream_cur = [];
ratio_aqd_freestream_urms = [];

ratio_vectrino_freestream_cur = [];
ratio_vectrino_freestream_urms = [];

avg=struct;
% 
% for idx=1:length(AQD11898.meanDepth_15min)
% [AQD11898.Wave.Lp(idx),AQD11898.Wave.kp(idx),AQD11898.Wave.sigmatm02(idx)]=disper(AQD11898.meanDepth_15min(idx),AQD11898.Wave.Tm02_ADP(idx));
% end
% clear idx
% plot(AQD11898.Wave.time,   AQD11898.Wave.sigmatm02);
% hold on
% AQD11898.Wave.Arms=AQD11898.Wave.ubr_sherwood./AQD11898.Wave.sigmatm02;
% plot(AQD11898.Wave.time,   AQD11898.Wave.Arms);
le_vectrino=[];
for ii=1:length(Vectrino);
    avg_vectrino = nanmean(Vectrino(ii).MeanVel(11:15));
    avg_vectrino_rms = nanmean(Vectrino(ii).ubr(11:15));
    
    [~,idx] = min(abs(AQD11898.time_ft-Vectrino(ii).InitialTime_Matlab));
    [~,idx_vector] = min(abs(ADV1670.time_ft-Vectrino(ii).InitialTime_Matlab));
    [~,idx_le] = min(abs(AQP11266.time_15min-Vectrino(ii).InitialTime_Matlab));
le_vectrino(ii) = le(idx_le);
    
    ratio_vectrino_freestream_cur (ii) = avg_vectrino/AQD11898.freestream(idx);
    ratio_vectrino_freestream_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,30:35));
%     ratio_aqd_freestream_cur (ii,:) = AQD11898.Mean_Current(idx,:)/AQD11898.freestream(idx);
%     ratio_aqd_freestream_urms (ii,:) = AQD11898.ubr(idx,:)/nanmean(AQD11898.ubr(idx,30:35));
    ratio_aqd_freestream_cur (ii,:) = AQD11898.Mean_Current(idx,:)/nanmean(AQD11898.Mean_Current(idx,30:35));
    ratio_aqd_freestream_urms (ii,:) = AQD11898.ubr(idx,:)/nanmean(AQD11898.ubr(idx,30:35));
hm0_aqd(ii)=AQD11898.Wave.Hm0_ADP(idx);
    
    if Vectrino(ii).BottomDistance<0.4 && Vectrino(ii).InitialTime_Matlab>737494.5
        violin_vectrino_vector_cur (ii) = avg_vectrino/ADV1670.Current_Mean(idx_vector);
        violin_vectrino_vector_urms (ii) = avg_vectrino_rms/ADV1670.ubr(idx_vector);
        
        violin_vectrino_aqd_cur (ii) = avg_vectrino/AQD11898.freestream(idx);
        violin_vectrino_aqd_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,30:35));
        violin_vectrino_aqd_Arms (ii) = AQD11898.ubr(idx,35).*AQD11898.Wave.Tp_ADP(idx)./(2*pi);
        %         plot(AQD11898.time_ft(idx),violin_vectrino_aqd_Arms(ii),'.r')
    else
        violin_vectrino_vector_cur (ii) = nan;
        violin_vectrino_vector_urms (ii) = nan;
        
        violin_vectrino_aqd_cur (ii) = nan;
        violin_vectrino_aqd_urms (ii) = nan;
        %         violin_vectrino_aqd_Arms (ii) = nan;
    end
end

%% test
% 
% for i=1:122
% plot(i,Vectrino(i).BottomDistance,'.k')
% hold on
% end
% 
% plot(violin_vectrino_aqd_Arms)
% plot(le*10)

%%
top_canopy_col = [221/256 96/256 49/256];
bot_canopy_col = [33/256 104/256 105/256];
bottom_col = [2/255 48/255 71/255];
free_stream_col = [255/255 183/255 3/255];
map_can = load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\bathycolormap.mat');

% map_bare = colorcet('L18', 'N', 64,'reverse',1);

f6=figure(6)
set(f6,'Units','centimeters','Position',[7.2496 -9.2251 13.1586 25.9997])

s5=subplot(2,2,1)
hold on
for i=35:45
     plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
plot(nanmean(ratio_aqd_freestream_urms(35:45,:)),AQD11898.z,'s',... 
         'color',bottom_col)
plot([0 1],[0.43 0.43],'--k')
plot([0 1],[nanmean(le_vectrino(35:45)) nanmean(le_vectrino(35:45))],'-k')
set(gca,'xtick',[0.25:0.25:1],'yscale','log','yticklabel',[],'fontsize',14)   
ylim([0.04 1.5])
xlim([0 1]) 
a=(nanmean(violin_vectrino_aqd_Arms(35:45)))
title([{'$$A_{\infty}^{rms} = 1.0 m/s$$' }], 'Interpreter','latex')
% title([num2str (nanmean(violin_vectrino_aqd_Arms(35:45)))])
grid on
box on
text(0.7,2.25,'Velocity profiles','fontsize',16)
text(0.55,-0.06,'Normalised velocity [ ]','fontsize',14)
set(gca,'xtick',[0.25:0.25:1],'yscale','log','ytick',[0.05 0.1 0.5 1],'yticklabel',[{'5','10','50','100'}], 'fontsize',12)  
yylabel=get(gca,'ytick')
ylim([0.04 1.5])
xlim([0 1]) 
ylabel('Elevation above bed [cm]','fontsize',14)
grid on
box on


s7=subplot(2,2,2)
hold on
for i=61:71
   p2=  plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
p4=plot(nanmean(ratio_aqd_freestream_urms(61:71,:)),AQD11898.z,'s',... 
         'color',bottom_col)
p5=plot([0 1],[0.43 0.43],'--k')
p6=plot([0 1],[nanmean(le_vectrino(61:71)) nanmean(le_vectrino(61:71))],'-k')

set(gca,'xtick',[0.25:0.25:1],'yscale','log','yticklabel',[],'fontsize',12)   
ylim([0.04 1.5])
xlim([0 1]) 
b=(nanmean(violin_vectrino_aqd_Arms(61:71)))

title(['$$A_{\infty}^{rms} = 0.55 m/s$$'], 'Interpreter','latex')
grid on
box on



l1=legend([p2 p5 p6],[{'$\hat{u}_{rms}/u_{rms,\infty}$','$l_v$','$l_{v,e}$'}],'interpreter','latex','fontsize',14,'orientation','vertical')
       %%
s9=subplot(2,2,3)
violin_vectrino_vector_cur(violin_vectrino_vector_cur>0.3)=nan;
[counts,centers] = hist(violin_vectrino_vector_cur,0:0.05:0.8);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',flip(map_can.bathycolormap),'showMM',2,...
    'ylabel','Attenuation parameter (\alpha [  ])');
% histogram(violin_vectrino_vector_cur,0:0.03:0.8,'normalization','probability')
xlabel('$\alpha_{c} [ \ ]$','Interpreter','Latex','fontsize',14)
ylim([0 1])
set(gca, 'YTick',0:.1:1,'xtick',[], 'XTickLabel',[],'fontsize',14)
xlim([0.5 1.5])
grid on
box on

s10=subplot(2,2,4)
violin_vectrino_vector_urms(violin_vectrino_vector_urms>0.95)=nan;
[counts,centers] = hist(violin_vectrino_vector_urms,0:0.05:1);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',flip(map_can.bathycolormap),'showMM',2,'ylabel',' ');
% histogram(violin_vectrino_vector_urms,0:0.05:0.8,'normalization','probability')
ylim([0 1])
set(gca, 'YTick',0:.1:1,'yticklabel',[],'xtick',[], 'XTickLabel',[],'fontsize',14)
xlim([0.5 1.5])
xlabel ('$\alpha^{rms} [ \ ]$','Interpreter','Latex','fontsize',14)
% caxis([0 maximo])
grid on
box on
text(-0.3,1.1,'Attenuation histograms','fontsize',15)



%%
s5.Position = [0.1322    0.5838    0.26    0.34];
s7.Position = [0.42 0.5838 0.26 0.34];
s9.Position = [0.1300    0.1100    0.3347    0.3412]
s10.Position = [0.5703    0.1100    0.3347    0.3412];

l1.Position = [0.6950    0.7810    0.2975    0.1444];

%%
print(f6,'C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Plots\chapter2\Figure12_alpha_TS_ICCE.tiff','-dtiff','-r300')

