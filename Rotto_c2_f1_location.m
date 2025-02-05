load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\bathy.mat')
map=load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\bathycolormap.mat');
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments')
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Quadpod_instruments.mat')
load('C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd.mat')

awac_8m_model.Veg_dense = load ('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\model\v0\awacforcingVeg\awac_8m.dat');
awac_8m_model.Veg_dense(:,1) = datenum(num2str(awac_8m_model.Veg_dense(:,1),'%f'),'yyyymmdd.HHMMSS');

load 'C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_RBRs_waves.mat'
for i=1:10
    ST = sprintf('ST%d',i);
    results = load (['C:\Users\NeryNeto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\model\awac_nofric\rbr',...
        num2str(i) '_spc.dat']);
    RBR_Model.NoFric.JP.(ST).time_wave = datenum(num2str(results(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.NoFric.JP.(ST).Hsig = results(:,2);
    RBR_Model.NoFric.JP.(ST).Hrms = results(:,2)./sqrt(2);
    RBR_Model.NoFric.JP.(ST).Tp = results(:,3);
    RBR_Model.NoFric.JP.(ST).MeanDir = results(:,4);
    RBR_Model.NoFric.JP.(ST).PeakDir = results(:,5);
    clear results
end
%% Load images 

% Deployment Sites
caminho = 'C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Plots\chapter2\OS\';
[dense_seagrass.x, dense_seagrass.y, dense_seagrass.z] = read_kml(fullfile(caminho,'dense_seagrass.kml'));
[dense_seagrass.x_utm, dense_seagrass.y_utm, dense_seagrass.zone] = deg2utm(dense_seagrass.y,dense_seagrass.x);
[intermediate_seagrass.x, intermediate_seagrass.y, intermediate_seagrass.z] = read_kml(fullfile(caminho,'intermediate_sg.kml'));
[intermediate_seagrass.x_utm, intermediate_seagrass.y_utm, intermediate_seagrass.zone] = deg2utm(intermediate_seagrass.y,intermediate_seagrass.x);
[intermediate_seagrass.x2, intermediate_seagrass.y2, intermediate_seagrass.z2] = read_kml(fullfile(caminho,'intermediate_sg2.kml'));
[intermediate_seagrass.x_utm2, intermediate_seagrass.y_utm2, intermediate_seagrass.zone] = deg2utm(intermediate_seagrass.y2,intermediate_seagrass.x2);
[sparse_seagrass.x, sparse_seagrass.y, sparse_seagrass.z] = read_kml(fullfile(caminho,'sparse_sg.kml'));
[sparse_seagrass.x_utm, sparse_seagrass.y_utm, sparse_seagrass.zone] = deg2utm(sparse_seagrass.y,sparse_seagrass.x);

sta_rbr=[363817	6458205; %RBR1  #124169 
363864	6458112; %RBR2  #124170 
363908	6458022; %RBR3  #124171 
363951	6457929; %RBR4  #77766  
363991	6457838; %RBR5  #78097  
364067	6457660; %RBR6  #77841  
364119	6457570; %RBR7  #78096  
364147	6457482; %RBR8  #77957  
364273	6457221; %RBR9  #124872 
364383	6456977]; %RBR10 #77765
sta_rbr=flipud(sta_rbr);

sta_awac_rbr = [364498 	6456752];%depth = -15
sta_awac_rbr_advection = [364497	6457759];%depth = -8.4
sta_aqua_off = [364223	6457342]; 

sta_sawhorse = [sta_rbr(10,1)+sind(25)*50 sta_rbr(10,2)-cosd(25)*50]; %depth = -3.2

path_symbols = 'C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Plots\chapter2\symbols';

[symbol_cur_profile.plrY1,symbol_cur_profile.map3,symbol_cur_profile.alphamap3] = imread(fullfile(path_symbols,'current_profile.png'));
[symbol_wave.plrY1,symbol_wave.map3,symbol_wave.alphamap3] = imread(fullfile(path_symbols,'Wave.png'));
[symbol_scaff.plrY1,symbol_scaff.map3,symbol_scaff.alphamap3] = imread(fullfile(path_symbols,'Scaffold.png'));
[symbol_sand.plrY1,symbol_sand.map3,symbol_sand.alphamap3] = imread(fullfile(path_symbols,'SAND.png'));

RBR.ST9.distance = 104.2; 
RBR.ST8.distance = 100.18;
RBR.ST7.distance = 102.46; 
RBR.ST6.distance = 99.4;
RBR.ST5.distance = 193.55;
RBR.ST4.distance = 103.94;
RBR.ST3.distance = 92.35;
RBR.ST2.distance = 289.82;
RBR.ST1.distance = 267.65;
%% Load profile
for num=1:122 %each profile
    Vec_h(num,1)=nanmean(Vectrino(num).BottomDistance);
end
heights = 0.05:0.01:max(Vec_h)+0.01;heights=heights';

free_stream_col = [255, 183, 3]/256;
top_canopy_col = [221, 96, 49]/256;
bot_canopy_col = [33, 104, 105]/256;
bottom_col = [2, 48, 71]/256;
%% Plot bathymetry figure
f1=figure(1);
set(f1,'Units','centimeters','OuterPosition',[0 -10 21.5635 24.8179]);
f1.Color=[1 1 1];
bathy.z(bathy.z>0)=0;
s1=axes('Position',[0.043577430972389 0.03 0.598379351740696 0.936192017259987]);
set(s1,'xtick',[],'ytick',[],'color',[233 233 235]./255,'XColor',[233 233 235]./255,'YColor',[233 233 235]./255)
t=annotation('textbox','color','w','Position',[0.0481512605042016 0.93 0.035 0.035], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(a)','HorizontalAlignment','center','FitBoxToText','off','fontsize',14)%
s2=axes('Position',[0.11 0.35 0.49 0.57],'xtick',[],'ytick',[],'color',[233 233 235]./255)
contourf(bathy.x,bathy.y,bathy.z,20)
colormap(map.bathycolormap)
hold on
c=colorbar;
caxis([-18 0])

set(gca,'xtick',363300:700:365400,'xticklabel',[363300:400:365400],'ytick',6456000:800:6458600,'yticklabel',...
    ({'6456000',     '6456500',     '6457000',     '6457500',     '6458000',     '6458500'}),'fontsize',14)
ytickangle(90)
set(gca,'fontsize',10)
title('Rottnest Island','color',[3 106 139]./255,'fontsize',16)
ylabel(c,'Depth (m)','fontsize',12,'color',[100 147 167]./255)
xlabel('Longitude(m)','fontsize',12,'color',[100 147 167]./255,'fontweight', 'bold' );
ylabel('Latitude (m)','fontsize',12,'color',[100 147 167]./255,'fontweight', 'bold' );

pdense=plot(dense_seagrass.x_utm,dense_seagrass.y_utm,'linewidth',2,'color',[139/256 0 0]);
pint=plot(intermediate_seagrass.x_utm,intermediate_seagrass.y_utm,'linewidth',2,'color',[220/256 20/256 60/256]);
psparse=plot(sparse_seagrass.x_utm,sparse_seagrass.y_utm,'linewidth',2,'color',[255/256 160/256 122/256]);
plot([363189.475 363689.475],-[-6456533 -6456533],...
    'k','linewidth',2);
text(363400,6456630,'1 km','color','k','BackgroundColor', 'w', 'fontsize',10);

xlim([363090 364850])
ylim([6456497 6458600])


%RBRs
for i=1:length(sta_rbr)
rbr1=plot(sta_rbr(i,1),sta_rbr(i,2),'.k','markersize',30);
hold on
text(sta_rbr(i,1)+70,sta_rbr(i,2),['P' num2str(i)],'BackgroundColor', 'w', 'fontsize',10);
end

%AWAC 8m
awac8=plot(sta_awac_rbr_advection(1),sta_awac_rbr_advection(2),'.y','markersize',30);
text(sta_awac_rbr_advection(1)+70,sta_awac_rbr_advection(2),'CW3','BackgroundColor', 'w', 'fontsize',10);

%Aquadopp 11 m
aqd11=plot(sta_aqua_off(1),sta_aqua_off(2),'.y','markersize',30);
text(sta_aqua_off(1)+70,sta_aqua_off(2),'CW2','BackgroundColor', 'w', 'fontsize',10);

%AWAC 15m
awac15=plot(sta_awac_rbr(1),sta_awac_rbr(2),'.y','markersize',30);
text(sta_awac_rbr(1)+70,sta_awac_rbr(2),'CW1','BackgroundColor', 'w', 'fontsize',10);

%Scaffolding
scaff=plot(sta_sawhorse(1),sta_sawhorse(2),'.r','markersize',30);
text(sta_sawhorse(1)-230,sta_sawhorse(2),'CW4','BackgroundColor', 'w', 'fontsize',10);

l1=legend([pdense, pint, psparse rbr1 awac8 scaff],[{'Dense','Intermediate','Sparse'},{'Pressure','Currents and Waves','Scaffolding'}],...
    'fontsize',10,'orientation','horizontal','location','southoutside','numcolumns',3)
l1.Position=[0.085 0.25 0.516506609481638 0.0384304214835553]
set(l1,'box','off')
%% RBR points
% s4=axes('Position',[0.043577430972389  0.0409924487594391 0.598379351740696 0.187432578209277]);
s4=axes('Position',[0.110144057623049 0.0784789644012945 0.450180072028812 0.13915857605178]);
set(s4,'xtick',[],'ytick',[])
t=annotation('textbox','color','w','Position',[0.0481512605042016 0.0364077669902913 0.035 0.035], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(b)','HorizontalAlignment','center','VerticalAlignment','middle','FitBoxToText','off','fontsize',14)%

for i=1:length(sta_rbr)
    ST = sprintf('ST%d',i);
    if i<10
        distances10 (i) = RBR.(ST).distance;
    end
    distance_cumulative_aux (i,2) = nanmean(RBR.(ST).meandepth);    

end
 
distances_with_awac=[251 distances10];
distance_cumulative_with_awac(:,1) = [0 cumsum(fliplr(distances_with_awac(1:end)))];
distance_cumulative_with_awac(end,2)=15.1;
distance_cumulative_with_awac(1:10,2)=flip(distance_cumulative_aux(:,2));

samplingRateIncrease = 3;
newXSamplePoints_awac = [linspace(min(distance_cumulative_with_awac(:,1)), max(distance_cumulative_with_awac(:,1)),...
    length(distance_cumulative_with_awac(:,1)) * samplingRateIncrease)];
smoothedY_awac = [spline(distance_cumulative_with_awac(:,1),...
    distance_cumulative_with_awac(:,2), newXSamplePoints_awac)];

distance_p10_sta_sawhorse = sqrt((sta_sawhorse(1)-sta_rbr(10,1)).^2+(sta_sawhorse(2)-sta_rbr(10,2)).^2);
distance_p10_aqdOff = sqrt((sta_aqua_off(1)-sta_rbr(10,1)).^2+(sta_aqua_off(2)-sta_rbr(10,2)).^2);

plot(newXSamplePoints_awac,(-smoothedY_awac),'-k')
hold on
plot(distance_cumulative_with_awac (:,1),(-distance_cumulative_with_awac (:,2)),'sk','MarkerFaceColor','k')
ylim([-16 0])
xlim([0 1700])
for i=1:length(sta_rbr)
    if i==2
        text(distance_cumulative_with_awac (i,1)-25,-distance_cumulative_with_awac (i,2)-1.5,['P' num2str(-i+11)],'fontsize',10)
    elseif i== 4
        text(distance_cumulative_with_awac (i,1)-70,-distance_cumulative_with_awac (i,2)-1.5,['P' num2str(-i+11)],'fontsize',10)
    elseif i==3
        text(distance_cumulative_with_awac (i,1)+20,-distance_cumulative_with_awac (i,2)+1,['P' num2str(-i+11)],'fontsize',10)
    elseif i==5
        text(distance_cumulative_with_awac (i,1),-distance_cumulative_with_awac (i,2)+1.5,['P' num2str(-i+11)],'fontsize',10)
    elseif i==7
        text(distance_cumulative_with_awac (i,1)+20,-distance_cumulative_with_awac (i,2)+.2,['P' num2str(-i+11)],'fontsize',10)
    elseif i==10
        text(distance_cumulative_with_awac (i,1)-100,-distance_cumulative_with_awac (i,2)-.5,['P' num2str(-i+11)],'fontsize',10)
    elseif i==1
        text(distance_cumulative_with_awac (i,1)+10,-distance_cumulative_with_awac (i,2)+1.5,['P' num2str(-i+11)],'fontsize',10)    
    else
        text(distance_cumulative_with_awac (i,1)+30,-distance_cumulative_with_awac (i,2)+.2,['P' num2str(-i+11)],'fontsize',10)
    end
end
grid on
xlabel('Distance to shallowest station (m)','fontsize',12)
ylabel('Depth (m)','fontsize',12)
set(gca,'ytick',[-14:2:0],'yticklabel',[14:-2:0])

patch('Faces',[1 2 3 4],'Vertices',[0 -16;1046 -16;1046 0;0 0],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[1046 -16;1249 -16;1249 0;1046 0],'FaceColor',[220/256 20/256 60/256],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[1249 -16;1400 -16;1400 0;1249 0],'FaceColor',[255/256 160/256 122/256],'FaceAlpha',.3)

plot (distance_p10_sta_sawhorse, -nanmean(AQD11898.meanDepth_15min)-1,'sk','MarkerFaceColor','r')
text (distance_p10_sta_sawhorse+20, double(-nanmean(AQD11898.meanDepth_15min)-.5),'CW4')
plot (distance_p10_aqdOff, -11.5,'sk','MarkerFaceColor','y')
text(distance_p10_aqdOff, -10.5,'CW2')
plot(distance_cumulative_with_awac (11,1),-distance_cumulative_with_awac (11,2),'sk','MarkerFaceColor','y')
text(distance_cumulative_with_awac (11,1)-40,-distance_cumulative_with_awac (11,2)+1,'CW1')

%% CW4 measurements
s3=axes('Position',[0.7250    0.65    0.25    0.32]);
set(s3,'xtick',[],'ytick',[])
t=annotation('textbox','color','w','Position',[0.648691476590636 0.93 0.035 0.035], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(c)','HorizontalAlignment','center','VerticalAlignment','middle','FitBoxToText','off','fontsize',14)%

AQD11898.PRESSURE_RAW(AQD11898.PRESSURE_RAW<2.5)=nan;
A = (max(AQD11898.PRESSURE_RAW)+min(AQD11898.PRESSURE_RAW))/2;
a= 0.33*(max(AQD11898.PRESSURE_RAW)-A);
p=area(0:0.01:1.2,a.*sin(10*(0:0.01:1.2)+pi/2)+A,'facealpha',0.3, 'Facecolor',[142, 202, 230]/256)
xlim([0 1.2])
ylim([0 nanmax(AQD11898.PRESSURE_RAW)+.3])

v = [1 AQD11898.z(30);1.15 AQD11898.z(30); 1.15 AQD11898.z(37);1 AQD11898.z(37)];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor','none','EdgeColor',free_stream_col,'linestyle','-.','linewidth',2) % yellow
hold on

v1 = [0.5 0.45;0.65 0.45; .65 0.55;0.5 0.55];
patch('Faces',f,'Vertices',v1,'FaceColor','none','EdgeColor','r','linestyle','-.','linewidth',2) %red

v2= [0.3 min(heights);.45 min(heights); .45 max(heights);0.3 max(heights)];
patch('Faces',f,'Vertices',v2,'FaceColor','none','EdgeColor','k','linestyle','-.','linewidth',2) %,'FaceAlpha',0.7)

v2= [0.8 min(heights);1 min(heights); 1 0.85;0.8 0.85];
patch('Faces',f,'Vertices',v2,'FaceColor','none','EdgeColor',[0 176 80]./256,'linestyle','-.','linewidth',2) %,'FaceAlpha',0.7)

% set(gca,'xticklabel',[])
box on
% vel = AQD11898.Mean_Current./AQD11898.freestream;
% e=errorbar(nanmean(vel(177:517,1:2:end)),AQD11898.z(1:2:end),nanstd(vel(177:517,1:2:end))*0,'horizontal','-s','MarkerSize',8,...
%     'MarkerEdgeColor',free_stream_col,'MarkerFaceColor',free_stream_col)
% e.Color = 'k'; hold on
% clear vel
% 
% vel_rms = AQD11898.ubr./nanmean(AQD11898.ubr(:,30:35),2);
% e=errorbar(nanmean(vel_rms(177:517,1:2:end)),AQD11898.z(1:2:end),nanstd(vel_rms(177:517,1:2:end))*0,'horizontal','-*','MarkerSize',8,...
%     'MarkerEdgeColor',free_stream_col,'MarkerFaceColor',free_stream_col,'linestyle','--')
% e.Color = 'k';
% clear vel vel_rms
% h=xlabel([{'$$\overline{u}/\overline{u}_{\infty} (-)$$'},...
%     {'$$u_{rms}/u_{rms,\infty} (-)$$'}])

% set(h,'Interpreter','latex','fontsize',14)
ylabel('cmab')
set(gca,'ytick',0:0.5:4,'yticklabel',0:50:400,'fontsize',10,'xtick',[])

symbol_sand.player3=imagesc([0 1.2],[-0.2 0.05],(symbol_sand.plrY1));
set(symbol_sand.player3,'AlphaData',symbol_sand.alphamap3);
%%
% figure(2)
s5=axes('Position',[0.725 0.08 0.24 0.24]);
s5.Position=[0.7250    0.32    0.2400    0.2023]
set(s5,'xtick',[],'ytick',[])
% s3=axes('Position',[0.725 0.4517 0.2412 0.5182]);

map2 = colorcet('R2', 'N', 64);
WindRose(RBR_Model.NoFric.JP.ST2.PeakDir-45,RBR_Model.NoFric.JP.ST2.Hsig,'lablegend','','legendvariable','H_s','centeredin0','True','titlestring',['Modelled wave rose', '(h \approx 5.5 m)'," "],...
    'speedround',1,'nspeeds',2,'nfreq',5,'legendtype',2,'freqlabelangle',45,'maxfrequency',100,'anglenorth',0,'angleeast',90,'cmap',map2, ...
    'labels',{'N','E','S','W'},'axes',gca);

d=annotation('textbox','color','w','Position',[0.648691476590636 0.4936 .035 0.035], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(d)','HorizontalAlignment','center','VerticalAlignment','middle','FitBoxToText','off','fontsize',14)%

e=annotation('textbox','color','w','Position',[0.648691476590636 0.0364077669902913 .035 0.035], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(e)','HorizontalAlignment','center','VerticalAlignment','middle','FitBoxToText','off','fontsize',14)%

% export_fig('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Plots\rose.tiff','-r450') 

%%
set(gcf, 'InvertHardcopy', 'off');  % Ensure the background color is preserved
print(gcf,'C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Plots\Figure1_location.tiff','-dtiff','-r300')
% saveas(gcf,'\\uniwa.uwa.edu.au\userhome\students2\22371812\My Documents\Documents\PhD\GardenIsland_exp1\article\plots\Figure8.png')
% print(gcf,'C:\Users\22371812\OneDrive - The University of Western Australia\GardenIslandpaper\Figure7.tiff',...
%     '-dtiff','-r300')
% print(gcf,'C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Plots\Figure1_location_aux','-dtiff','-r900')
% export_fig('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Plots\Figure1_location_aux.tiff','-r900') 
