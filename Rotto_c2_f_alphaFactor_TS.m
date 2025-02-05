% BrickerMonismith = load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\functions\BrickerMonismith.mat');
% BrickerMonismith_vectrino = load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\functions\BrickerMonismith_vectrino.mat');
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Quadpod_instruments.mat');
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments')
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
for ii=1:length(Vectrino);
    avg_vectrino=nanmean(Vectrino(ii).MeanVel(11:15));
    avg_vectrino_rms=nanmean(Vectrino(ii).ubr(11:15));
    
    [~,idx] = min(abs(AQD11898.time_ft-Vectrino(ii).InitialTime_Matlab));
    ratio_vectrino_freestream_cur (ii) = avg_vectrino/AQD11898.freestream(idx);
    ratio_vectrino_freestream_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,30:35));
    ratio_aqd_freestream_cur (ii,:) = AQD11898.Mean_Current(idx,:)/AQD11898.freestream(idx);
    ratio_aqd_freestream_urms (ii,:) = AQD11898.ubr(idx,:)/nanmean(AQD11898.ubr(idx,30:35));
    
    if Vectrino(ii).BottomDistance<00.1 && Vectrino(ii).InitialTime_Matlab>737494.5
        violin_vectrino_freestream_cur (ii) = avg_vectrino/AQD11898.freestream(idx);
        violin_vectrino_freestream_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,30:35));
    else
        violin_vectrino_freestream_cur (ii) = nan;
        violin_vectrino_freestream_urms (ii) = nan;
    end
    if Vectrino(ii).BottomDistance<0.4 && Vectrino(ii).InitialTime_Matlab>737494.5
        violin_vectrino_aqd_cur (ii) = avg_vectrino/nanmean(AQD11898.Mean_Current(idx,1:3));
        violin_vectrino_aqd_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,1:3));
    else
        violin_vectrino_aqd_cur (ii) = nan;
        violin_vectrino_aqd_urms (ii) = nan;
    end
    
end


%%
top_canopy_col = [221/256 96/256 49/256];
bot_canopy_col = [33/256 104/256 105/256];
bottom_col = [2/255 48/255 71/255];
free_stream_col = [255/255 183/255 3/255];
map_can = colorcet('L5', 'N', 64,'reverse',1);
map_bare = colorcet('L18', 'N', 64,'reverse',1);

f6=figure(6)
set(f6,'Units','centimeters','Position',[0 -5 18 25])

s1=subplot(2,8,1)
hold on
for i=23:34
plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col, 'markeredgecolor',bottom_col);
end
plot(nanmean(ratio_aqd_freestream_cur(23:34,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
plot(nanmean(ratio_aqd_freestream_urms(23:34,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
set(gca,'xtick',[0:0.4:1],'yscale','log','fontsize',10)  
ylim([0.04 1.5])
xlim([0 0.9])
title('Profile 1')
ylabel('mab')
grid on
box on


s2=subplot(2,8,2)
hold on
for i=35:46
plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
plot(nanmean(ratio_aqd_freestream_cur(35:46,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
plot(nanmean(ratio_aqd_freestream_urms(35:46,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
set(gca,'xtick',[0.4:0.4:1],'yscale','log','yticklabel',[],'fontsize',10)  
ylim([0.04 1.5])
xlim([0 0.9])
title('Profile 2')
grid on
box on


s3=subplot(2,8,3)
hold on
for i=47:58
plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
plot(nanmean(ratio_aqd_freestream_cur(47:58,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
plot(nanmean(ratio_aqd_freestream_urms(47:58,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
set(gca,'xtick',[0.4:0.4:1],'yscale','log','yticklabel',[],'fontsize',10)  
ylim([0.04 1.5])
xlim([0 0.9])
title('Profile 3')
grid on
box on


s4=subplot(2,8,4)
hold on
for i=61:72
plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
plot(nanmean(ratio_aqd_freestream_cur(61:72,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
plot(nanmean(ratio_aqd_freestream_urms(61:72,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
set(gca,'xtick',[0.4:0.4:1],'yscale','log','yticklabel',[],'fontsize',10)  
ylim([0.04 1.5])
xlim([0 0.9])
title('Profile 4')
grid on
box on


s5=subplot(2,8,5)
hold on
for i=75:86
plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
plot(nanmean(ratio_aqd_freestream_cur(75:86,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
plot(nanmean(ratio_aqd_freestream_urms(75:86,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
set(gca,'xtick',[0.4:0.4:1],'yscale','log','yticklabel',[],'fontsize',10)  
ylim([0.04 1.5])
xlim([0 0.9])
title('Profile 5')
grid on
box on


s6=subplot(2,8,6)
hold on
for i=87:98
plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
plot(nanmean(ratio_aqd_freestream_cur(87:98,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
plot(nanmean(ratio_aqd_freestream_urms(87:98,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
set(gca,'xtick',[0.4:0.4:1],'yscale','log','yticklabel',[],'fontsize',10)  
ylim([0.04 1.5])
xlim([0 0.9])
title('Profile 6')
grid on
box on


s7=subplot(2,8,7)
hold on
for i=99:110
plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
plot(nanmean(ratio_aqd_freestream_cur(99:110,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
plot(nanmean(ratio_aqd_freestream_urms(99:110,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
set(gca,'xtick',[0.4:0.4:1],'yscale','log','yticklabel',[],'fontsize',10)  
ylim([0.04 1.5])
xlim([0 0.9])
title('Profile 7')
grid on
box on


s8=subplot(2,8,8)
hold on
for i=111:122
p1=plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     p2=plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
p3=plot(nanmean(ratio_aqd_freestream_cur(111:122,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
p4=plot(nanmean(ratio_aqd_freestream_urms(111:122,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
set(gca,'xtick',[0.4:0.4:1],'yscale','log','yticklabel',[],'fontsize',10)
ylim([0.04 1.5])
xlim([0 0.9])
title('Profile 8')
grid on
box on
l1=legend([p1 p2 p3 p4],[{'$\overline{\hat{u}}/\overline{u}_{\infty}$','$\hat{u}^R_{w}/u^R_{w,\infty}$',...
           '$\overline{u}_{c}/\overline{u}_{c,\infty}$','$u^R_{w}/u^R_{w,\infty}$'}],'interpreter','latex')
sgtitle({'$\overline{u}$ and $u_{rms}$ profiles'},...
    'interpreter','latex','fontsize',15)
       %%
s9=subplot(2,8,9)
[counts,centers] = hist(violin_vectrino_freestream_cur,0:.02:0.8);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',map_can,'showMM',2,...
    'ylabel','Attenuation parameter (\alpha)');
xlabel('$\hat{\overline{u}}/u_{top}$','Interpreter','Latex')
ylim([0 1])
set(gca, 'YTick',0:.1:1,'xtick',[], 'XTickLabel',[],'fontsize',12)
xlim([0.5 1.5])
grid on
box on

s10=subplot(2,8,12)
[counts,centers] = hist(violin_vectrino_freestream_urms,0:.02:0.8);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',map_bare,'showMM',2,'ylabel',' ');
ylim([0 1])
xlabel('$\hat{u}_{rms}/u_{rms,top}$','Interpreter','Latex','fontsize',12)
set(gca, 'YTick',0:.1:1,'yticklabel',[],'xtick',[], 'XTickLabel',[],'fontsize',12)
xlim([0.5 1.5])
% caxis([0 maximo])
grid on
box on
text(0.75,1.1,'Attenuation histogram','fontsize',15)


s11=subplot(2,8,14)
[counts,centers] = hist(violin_vectrino_aqd_cur,0:.02:1);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',map_can,'showMM',2,'ylabel',' ');
xlabel('$\alpha_{c}$','Interpreter','Latex','fontsize',12)
ylim([0 1])
set(gca, 'YTick',0:.1:1,'yticklabel',[],'xtick',[], 'XTickLabel',[],'fontsize',12)
xlim([0.5 1.5])
grid on
box on
% text(0.55,1.07,'(D)','fontsize',12)
% text(0.55,.05,{['$$\alpha_c = $$ ' num2str(round(avg1{2}(1)*100)/100)]},'Interpreter', 'LaTeX','fontsize',12)

s12=subplot(2,8,16)
[counts,centers] = hist(violin_vectrino_aqd_urms,0:.02:1);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',map_can,'showMM',2,'ylabel',' ');
ylim([0 1])
set(gca, 'YTick',0:.1:1,'yticklabel',[])
xlim([0.5 1.5])
set(gca,'xtick',[], 'XTickLabel',[])
xlabel('$\alpha^R_{w}$','Interpreter','Latex','fontsize',12)
grid on
box on
% text(0.55,1.07,'(D)','fontsize',12)
% text(0.55,.05,{['$$\alpha_c = $$ ' num2str(round(avg1{2}(1)*100)/100)]},'Interpreter', 'LaTeX','fontsize',12)
s1.Position = [0.08 0.4 0.1 0.5];
s2.Position = [0.19 0.4 0.1 0.5];
s3.Position = [0.3 0.4 0.1 0.5];
s4.Position = [0.41 0.4 0.1 0.5];
s5.Position = [0.52 0.4 0.1 0.5];
s6.Position = [0.63 0.4 0.1 0.5];
s7.Position = [0.74 0.4 0.1 0.5];
s8.Position = [0.85 0.4 0.1 0.5];
s9.Position = [0.08 0.06 0.16 0.25];
s10.Position = [0.255 0.06 0.16 0.25];
s11.Position = [0.425 0.06 0.16 0.25];
s12.Position = [0.595 0.06 0.16 0.25];
l1.Position = [0.7863    0.2237    0.1667    0.0866];

%%
print(f6,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\Figure6_alpha.tiff','-dtiff','-r300')


%% LESS PROFILES

% BrickerMonismith = load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\functions\BrickerMonismith.mat');
% BrickerMonismith_vectrino = load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\functions\BrickerMonismith_vectrino.mat');
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Quadpod_instruments.mat');
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments')
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

for idx=1:length(AQD11898.meanDepth_15min)
[AQD11898.Wave.Lp(idx),AQD11898.Wave.kp(idx),AQD11898.Wave.sigmatm02(idx)]=disper(AQD11898.meanDepth_15min(idx),AQD11898.Wave.Tm02_ADP(idx));
end
clear idx
plot(AQD11898.Wave.time,   AQD11898.Wave.sigmatm02);
hold on
AQD11898.Wave.Arms=AQD11898.Wave.ubr_sherwood./AQD11898.Wave.sigmatm02;
plot(AQD11898.Wave.time,   AQD11898.Wave.Arms);
for ii=1:length(Vectrino);
    avg_vectrino=nanmean(Vectrino(ii).MeanVel(11:15));
    avg_vectrino_rms=nanmean(Vectrino(ii).ubr(11:15));
    
    [~,idx] = min(abs(AQD11898.time_ft-Vectrino(ii).InitialTime_Matlab));
    ratio_vectrino_freestream_cur (ii) = avg_vectrino/AQD11898.freestream(idx);
    ratio_vectrino_freestream_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,30:35));
    ratio_aqd_freestream_cur (ii,:) = AQD11898.Mean_Current(idx,:)/AQD11898.freestream(idx);
    ratio_aqd_freestream_urms (ii,:) = AQD11898.ubr(idx,:)/nanmean(AQD11898.ubr(idx,30:35));
    
    if Vectrino(ii).BottomDistance<00.15 && Vectrino(ii).InitialTime_Matlab>737494.5
        violin_vectrino_freestream_cur (ii) = avg_vectrino/AQD11898.freestream(idx);
        violin_vectrino_freestream_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,30:35));
          
        violin_vectrino_aqd_Arms (ii) = violin_vectrino_freestream_urms(ii).*AQD11898.Wave.Tm02_ADP(idx)./(2*pi);
        plot(AQD11898.time_ft(idx),violin_vectrino_aqd_Arms(ii),'.r')
    else
        violin_vectrino_freestream_cur (ii) = nan;
        violin_vectrino_freestream_urms (ii) = nan;
        violin_vectrino_aqd_Arms (ii) = nan;
    end
    if Vectrino(ii).BottomDistance<0.15 && Vectrino(ii).InitialTime_Matlab>737494.5
        violin_vectrino_aqd_cur (ii) = avg_vectrino/nanmean(AQD11898.Mean_Current(idx,1:3));
        violin_vectrino_aqd_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,1:3));
      
    else
        violin_vectrino_aqd_cur (ii) = nan;
        violin_vectrino_aqd_urms (ii) = nan;
    end
    
end

%%
top_canopy_col = [221/256 96/256 49/256];
bot_canopy_col = [33/256 104/256 105/256];
bottom_col = [2/255 48/255 71/255];
free_stream_col = [255/255 183/255 3/255];
map_can = colorcet('L5', 'N', 64,'reverse',1);
map_bare = colorcet('L18', 'N', 64,'reverse',1);

f6=figure(6)
set(f6,'Units','centimeters','Position',[0 -5 18 25])

s1=subplot(2,8,1)
hold on
for i=23:34
plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col, 'markeredgecolor',bottom_col);
end
plot(nanmean(ratio_aqd_freestream_cur(23:34,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
plot(nanmean(ratio_aqd_freestream_urms(23:34,:)),AQD11898.z,'s',... 
         'color',free_stream_col)

set(gca,'xtick',[0.25:0.25:1],'yscale','log','ytick',[0.05 0.1 0.5 1],'yticklabel',[{'5','10','50','100'}], 'fontsize',10)  
yylabel=get(gca,'ytick')
ylim([0.04 1.5])
xlim([0 1]) 
title('Profile 1')
ylabel('cmab')
grid on
box on



s5=subplot(2,8,5)
hold on
for i=75:86
plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
plot(nanmean(ratio_aqd_freestream_cur(75:86,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
plot(nanmean(ratio_aqd_freestream_urms(75:86,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
set(gca,'xtick',[0.25:0.25:1],'yscale','log','yticklabel',[],'fontsize',10)   
ylim([0.04 1.5])
xlim([0 1]) 
title('Profile 5')
xlabel('$$u \hspace{0.1cm} ratios$$','interpreter','latex','fontsize',16)
grid on
box on


s7=subplot(2,8,7)
hold on
for i=99:110
plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
plot(nanmean(ratio_aqd_freestream_cur(99:110,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
plot(nanmean(ratio_aqd_freestream_urms(99:110,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
set(gca,'xtick',[0.25:0.25:1],'yscale','log','yticklabel',[],'fontsize',10)   
ylim([0.04 1.5])
xlim([0 1]) 
title('Profile 7')
grid on
box on



s8=subplot(2,8,8)
hold on
for i=111:122
p1=plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
     p2=plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
p3=plot(nanmean(ratio_aqd_freestream_cur(111:122,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
p4=plot(nanmean(ratio_aqd_freestream_urms(111:122,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
set(gca,'xtick',[0.25:0.25:1],'yscale','log','yticklabel',[],'fontsize',10)  
ylim([0.04 1.5])
xlim([0 1]) 
title('Profile 8')
grid on
box on
l1=legend([p1 p2 p3 p4],[{'$\hat{\overline{u}}/\overline{u}_{\infty}$','$\hat{u}^R/u^R_{\infty}$',...
           '$\overline{u}_{c}/\overline{u}_{c,\infty}$','$u^R/u^R_{\infty}$'}],'interpreter','latex','fontsize',14)
h=sgtitle({'$\overline{u}$ and $u^R$ profiles'},...
    'interpreter','latex','fontsize',15)
       %%
s9=subplot(2,8,9)
[counts,centers] = hist(violin_vectrino_freestream_cur,0:0.05:0.8);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',map_can,'showMM',2,...
    'ylabel','Attenuation parameter (\alpha)');
xlabel('$\hat{\overline{u}}/u_{top}$','Interpreter','Latex')
ylim([0 1])
set(gca, 'YTick',0:.1:1,'xtick',[], 'XTickLabel',[],'fontsize',12)
xlim([0.5 1.5])
grid on
box on

s10=subplot(2,8,12)
[counts,centers] = hist(violin_vectrino_freestream_urms,0:0.05:0.8);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',map_bare,'showMM',2,'ylabel',' ');
ylim([0 1])
xlabel('$\hat{u}^R/u_{top}^R$','Interpreter','Latex','fontsize',12)
set(gca, 'YTick',0:.1:1,'yticklabel',[],'xtick',[], 'XTickLabel',[],'fontsize',12)
xlim([0.5 1.5])
% caxis([0 maximo])
grid on
box on
text(0.75,1.1,'Attenuation histogram','fontsize',15)


s11=subplot(2,8,14)
[counts,centers] = hist(violin_vectrino_aqd_cur,0:0.05:1);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',map_can,'showMM',2,'ylabel',' ');
xlabel('$\alpha_{c}$','Interpreter','Latex','fontsize',12)
ylim([0 1])
set(gca, 'YTick',0:.1:1,'yticklabel',[],'xtick',[], 'XTickLabel',[],'fontsize',12)
xlim([0.5 1.5])
grid on
box on
% text(0.55,1.07,'(D)','fontsize',12)
% text(0.55,.05,{['$$\alpha_c = $$ ' num2str(round(avg1{2}(1)*100)/100)]},'Interpreter', 'LaTeX','fontsize',12)

s12=subplot(2,8,16)
[counts,centers] = hist(violin_vectrino_aqd_urms,0:0.05:1);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',map_can,'showMM',2,'ylabel',' ');
ylim([0 1])
set(gca, 'YTick',0:.1:1,'yticklabel',[])
xlim([0.5 1.5])
set(gca,'xtick',[], 'XTickLabel',[])
xlabel('$\alpha^R$','Interpreter','Latex','fontsize',12)
grid on
box on
% text(0.55,1.07,'(D)','fontsize',12)
% text(0.55,.05,{['$$\alpha_c = $$ ' num2str(round(avg1{2}(1)*100)/100)]},'Interpreter', 'LaTeX','fontsize',12)
%%
plot (violin_vectrino_aqd_Arms,violin_vectrino_aqd_urms,'.','markersize',10)
hold on
plot(violin_vectrino_aqd_Arms,2.452.*violin_vectrino_aqd_Arms+0.09382,'r')
xlabel('$A^{rms}_{\infty}=u_{\infty}*T/2\pi$','Interpreter','Latex','fontsize',14)
ylabel('$\alpha$','Interpreter','Latex','fontsize',14)
grid on
text(0.056,0.63,'R^2=0.8','fontsize',14)
xlim([0 0.3])
ylim([0 0.8])
arms=AQD11898.Wave.Arms(182:502);
ratio=AQD11898.ubr(182:502,1)./AQD11898.ubr(182:502,35);
plot(AQD11898.Wave.Arms(182:502), AQD11898.ubr(182:502,1)./AQD11898.ubr(182:502,35),'.k')
% print(gcf,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\Figure6_alpha_Arms.tiff','-dtiff','-r300')

%%
s1.Position = [0.08 0.42 0.2 0.49];
s5.Position = [0.29 0.42 0.2 0.49];
s7.Position = [0.5 0.42 0.2 0.49];
s8.Position = [0.71 0.42 0.2 0.49];
s9.Position = [0.08 0.06 0.16 0.25];
s10.Position = [0.255 0.06 0.16 0.25];
s11.Position = [0.425 0.06 0.16 0.25];
s12.Position = [0.595 0.06 0.16 0.25];
l1.Position = [0.7863    0.206    0.1667    0.0866];

%%
print(f6,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\Figure6_alpha.tiff','-dtiff','-r300')


