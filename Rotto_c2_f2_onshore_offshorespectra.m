load('C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m\AWAC_6734.mat')
% load('C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Quadpod_instruments.mat')
wind = readmatrix('C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\wind_2019_hilarys.csv');
load('C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd.mat')
%%
caminho='C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m';
fileList = dir(fullfile(caminho,'*.was'));
pwrFreq = importdata(fullfile(fileList.folder,fileList.name));clear fileList

fileList = dir(fullfile(caminho,'*.wds'));
pwrFreqDir = importdata(fullfile(fileList.folder,fileList.name));clear fileList

fileList = dir(fullfile(caminho,'*.wdr'));
pwrDir = importdata(fullfile(fileList.folder,fileList.name));clear fileList

fileList = dir(fullfile(caminho,'*.sen'));
sen = importdata(fullfile(fileList.folder,fileList.name));

tempo=datenum(sen(1,3),sen(1,1),sen(1,2),sen(1,4),sen(1,5),sen(1,6)):1/48:datenum(sen(end,3),sen(end,1),sen(end,2),sen(end,4),sen(end,5),sen(end,6));
tempo=tempo(1:1782);

map1 = colorcet('L16', 'N', 128, 'reverse',0);
map2 = colorcet('L16', 'N', 128, 'reverse',1);
%% Organize data
eli_antes=0;eli_dep=129;
pwrFreq(pwrFreq == -9.000000) = NaN; 
pwrFreqDir(pwrFreqDir == -9.000000) = NaN;
waveData.Frequency = pwrFreq(1,:)';
waveData.pwrSpectrum           = pwrFreq(eli_antes+2:end-eli_dep,:); 
waveData.fullSpectrumDirection = 4:4:360;
nFreqs   = size(pwrDir,2);
% waveData.fullSpectrum          = pwrFreqDir((eli_antes*nFreqs)+1:end-eli_dep*nFreqs,:);
waveData.fullSpectrum          = pwrFreqDir((eli_antes*nFreqs)+1:end-eli_dep*nFreqs,:);
% waveData.fullSpectrum (end,:) = [];
waveData.Frequency=waveData.Frequency(1:nFreqs,:);
nDirs    = length(waveData.fullSpectrumDirection);
nSamples = length(waveData.fullSpectrum) / nFreqs;
aux=waveData.fullSpectrum;
% aux=aux(:,1:2304);b=reshape(aux, nDirs, nFreqs, nSamples);
waveData.fullSpectrum = permute(reshape(aux, nDirs, nFreqs, round(nSamples)), [ 2 1 3]);

%% Unnormalize data
% E(f,theta) = S(f) * d(f,theta)     ,   with S= Energy spectra (*.was)    and d=full directional spectra (*.wds)
E=zeros(nFreqs,nDirs,nSamples);
for t=1:nSamples
    for f=1:nFreqs;
        for theta=1:nDirs;
            E(f,theta,t) = waveData.pwrSpectrum(t,f) * waveData.fullSpectrum(f,theta,t);
        end
    end
end

%% Aquadopp data
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments.mat');

clear caminho
caminho='C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Aquadopps\AQD2985';
fileList = dir(fullfile(caminho,'*.was'));
pwrFreq_aqd = importdata(fullfile(fileList.folder,fileList.name));clear fileList

fileList = dir(fullfile(caminho,'*.wds'));
pwrFreqDir_aqd = importdata(fullfile(fileList.folder,fileList.name));clear fileList

fileList = dir(fullfile(caminho,'*.wdr'));
pwrDir_aqd = importdata(fullfile(fileList.folder,fileList.name));clear fileList

fileList = dir(fullfile(caminho,'*.sen'));
sen_aqd = importdata(fullfile(fileList.folder,fileList.name));

tempo_aqd=datenum(sen_aqd(1,3),sen_aqd(1,1),sen_aqd(1,2),sen_aqd(1,4),sen_aqd(1,5),sen_aqd(1,...
    6)):1/48:datenum(sen_aqd(end,3),sen_aqd(end,1),sen_aqd(end,2),sen_aqd(end,4),sen_aqd(end,5),sen_aqd(end,6));
%% Plots
f2=figure(2)
set(f2,'Units','centimeters','OuterPosition',[5 1 18 19])

s1=subplot(4,1,1) %Wind
winddir = 270-wind(:,9);
windvel = wind(:,11);
[u,v]=veldir2uv(windvel,winddir);
map = colorcet('R2', 'N', 64);
colormap(map);

% quiverc([wind(1000:3000,2); 7.375158824952062e+05],[zeros(numel(wind(1000:3000,2)),1); 0.5],[u(1000:3000); -10],[v(1000:3000); 0]);
quiverc(wind(1000:3000,2),zeros(numel(wind(1000:3000,2)),1),u(1000:3000),v(1000:3000));

box on
h1=colorbar ;
set(gcf, 'InvertHardCopy', 'on')

% ind_vel = nanmax(windvel);
set(gca, 'CLim', [0, 15])
xticks=737492:1:737519;

% yticks=0:2.5:15;
% x1=get(h1,'position');
% x1(1)=0.2;x1(2)=0.787;x1(3)=0.2;
% set(h1,'position',x1)
% set(h1, 'YTick',yticks,'YTickLabel',yticks,'Fontsize',12)
% title('Wind speed [m/s]','Fontsize',14);
xlim([7.374923809201388e+05 7.375184325896990e+05]);
ylim([-0.6 1])
grid on
set(gca,'ytick',[],'YTicklabel',[],'XTick',xticks,'xticklabel',[])
datetick('x','dd','keepticks','keeplimits')

ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
% text(7.375158824952062e+05-1.2,0.7,'10 m/s')
h1.Ticks=0:3:15;
colorTitleHandle = get(h1,'Title');
titleString = 'Wind speed [m/s]';
set(colorTitleHandle ,'String',titleString);
h1.FontSize=10;
t=annotation('textbox','color','w','Position',[0.009 0.9501 0.027 0.037], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(a)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle')%)%

%%
s2=subplot(4,1,2) %offshore spectra

map2 = colorcet('L16', 'N', 128, 'reverse',0);
% map=load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\functions\bathycolormap.mat');
pcolor(tempo(1:end),waveData.Frequency(1:end),log(waveData.pwrSpectrum(:,1:98)'));
shading interp; colormap(s2, map2);
% caxis([0 2]);
caxis([-5 2]);
xlim([7.374923809201388e+05 7.375184325896990e+05]);
ylim([0.02 0.5]);
yticks=0:.1:.5;
% ccol=colorbar;
% ccol.Limits;
% title(ccol,{'log (energy)', '[m^2/Hz]'},'Fontsize',10);
ylabel('Frequecy [Hz]'); 
set(gca,'fontsize',10,'XTick',xticks,'xticklabel',[],'ytick',yticks,'YTickLabel',0:0.1:0.5)
t2=annotation('textbox','color','w','Position',[0.009 0.7 0.027 0.037], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(b)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle');%)%
grid on
title('Offshore wave spectra (CW1)','fontsize',14) 
set(gca,'layer','top')
%%
% s3=subplot(4,1,3)
% p=area(AWAC_6734.wave.time,AWAC_6734.wave.Hm0, 'Facecolor',[142, 202, 230]/256)
% hold on
% p2=patch('Faces',[1 2 3 4],'Vertices',[Vectrino(23).Time(1) 0;Vectrino(34).Time(end) 0;Vectrino(34).Time(end) 2.5;Vectrino(23).Time(1) 2.5],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
% patch('Faces',[1 2 3 4],'Vertices',[Vectrino(35).Time(1) 0;Vectrino(58).Time(end) 0;Vectrino(58).Time(end) 2.5;Vectrino(35).Time(1) 2.5],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
% patch('Faces',[1 2 3 4],'Vertices',[Vectrino(59).Time(1) 0;Vectrino(74).Time(end) 0;Vectrino(74).Time(end) 2.5;Vectrino(59).Time(1) 2.5],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
% patch('Faces',[1 2 3 4],'Vertices',[Vectrino(75).Time(1) 0;Vectrino(122).Time(end) 0;Vectrino(122).Time(end) 2.5;Vectrino(75).Time(1) 2.5],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
% xlim([7.374923809201388e+05 7.375184325896990e+05]);
% ylim([0.02 2.5])
% title('Offshore wave height (CW1)','Fontsize',14);
% set(gca,'fontsize',10,'XTick',xticks,'xticklabel',[])
% ylabel('H_{s,CW1} [m]')
% grid on
% t3=annotation('textbox','color','w','Position',[0.009 0.48 0.027 0.037], 'BackgroundColor',[100 147 167]./255,...
%     'EdgeColor',[100 147 167]./255,'String','C','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle');
% legend([p,p2],[{'Wave height (m)','In-canopy sampling'}],'fontsize',10)

%%
% s4=subplot(4,1,4)
% p=area(AWAC_6734.wave.time,AWAC_6734.wave.Tp)
% hold on
% p2=area(AWAC_6734.wave.time,AWAC_6734.wave.Tm02, 'Facecolor',[142, 202, 230]/256)
% p3=patch('Faces',[1 2 3 4],'Vertices',[Vectrino(23).Time(1) 0;Vectrino(34).Time(end) 0;Vectrino(34).Time(end) 18;Vectrino(23).Time(1) 18],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
% patch('Faces',[1 2 3 4],'Vertices',[Vectrino(35).Time(1) 0;Vectrino(58).Time(end) 0;Vectrino(58).Time(end) 18;Vectrino(35).Time(1) 18],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
% patch('Faces',[1 2 3 4],'Vertices',[Vectrino(59).Time(1) 0;Vectrino(74).Time(end) 0;Vectrino(74).Time(end) 18;Vectrino(59).Time(1) 18],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
% patch('Faces',[1 2 3 4],'Vertices',[Vectrino(75).Time(1) 0;Vectrino(122).Time(end) 0;Vectrino(122).Time(end) 18;Vectrino(75).Time(1) 18],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
% xlim([7.374923809201388e+05 7.375184325896990e+05]);
% ylim([0.02 18])
% title('Offshore wave period (CW1)','Fontsize',14);
% set(gca,'fontsize',10,'XTick',xticks)
% datetick('x','dd','keepticks','keeplimits')
% ylabel('Wave period [s]')
% grid on
% t4=annotation('textbox','color','w','Position',[0.01 0.25 0.027 0.037], 'BackgroundColor',[100 147 167]./255,...
%     'EdgeColor',[100 147 167]./255,'String','D','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle')%
% l=legend([p,p2],[{'T_p','T_{m02}'}],'fontsize',10,'orientation','vertical')
% 
% ax = gca;
% labels = string(ax.XAxis.TickLabels); % extract
% labels(2:2:end) = nan; % remove every other one
% ax.XAxis.TickLabels = labels; % set
% xlabel('Date in March and April 2019')

s3=subplot(4,1,3);
pcolor(AQD11898.time_15min,AQD11898.Wave.f_ADP(1,:),log(AQD11898.Wave.Spp_ADP)');
hold on
pcolor(RBR.ST2.time_wave,RBR.ST9.f_ss(:,1)',log(RBR.ST9.Syy));

shading interp;
colormap(s3,map1);
% caxis([0 .05])
caxis([-5 2]);
ylim([0.02 0.5]);
xlim([AQD11898.time_15min(1) RBR.ST2.time_wave(end)]);
yticks=0:.1:.5;
ylabel('Frequecy [Hz]')
xticks=737492:1:737519;
set(gca,'YTick',yticks,'XTick',xticks,'xticklabel',[],'fontsize',10)
col=colorbar;
title(col,{'log (energy)', '[m^2/Hz]'},'Fontsize',10);
title('Onshore wave spectra (CW4)','Fontsize',14)
grid on
set(gca,'layer','top')
t3=annotation('textbox','color','w','Position',[0.009 0.48 0.027 0.037], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(c)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle');
%% 
s4=subplot(4,1,4);
awac_interp=interp1(AWAC_6734.wave.time,AWAC_6734.wave.Hm0,AQD11898.time_15min);
awac_interp_rbr=interp1(AWAC_6734.wave.time,AWAC_6734.wave.Hm0,RBR.ST2.time_wave);
ratio1=AQD11898.Wave.Hm0_ADP./awac_interp;ratio1(ratio1>0.7)=nan;
ratio2=RBR.ST8.Hm0./awac_interp_rbr;ratio2(ratio2>0.7)=nan;

AQD11898.Wave.Hm0_ADP(AQD11898.Wave.Hm0_ADP>1 | AQD11898.Wave.Hm0_ADP<0)=nan;
p1=area(AWAC_6734.wave.time,AWAC_6734.wave.Hm0)
hold on

area(AQD11898.time_15min,AQD11898.Wave.Hm0_ADP, 'Facecolor',[142, 202, 230]/256)
p2=area(RBR.ST2.time_wave,RBR.ST9.Hm0, 'Facecolor',[142, 202, 230]/256);

p3=plot(AQD11898.time_15min,ratio1,'r','linewidth',1.5);
plot(RBR.ST2.time_wave,ratio2,'r','linewidth',1.5);

xlim([AQD11898.time_15min(1) RBR.ST2.time_wave(end)]);
% xlim([AQD11898.time_ft(1) AQD11898.time_ft(end)])
ylim([0.02 2.5])
% yh=ylabel({'H_{s,CW4} [m]';' '},'color','k');
yh=ylabel('H_{s} [m]','color','k');
% set(yh,'Units','Normalized');
% ylabPos = get(yh,'Position');
% text(ylabPos(1)-.01,ylabPos(2),'H_{s,CW4}/H_{s,CW1} [ ]',  'Units','Normalized',  'HorizontalAlignment','center',...
%     'Vert','bottom', 'Rotation',90,'color','r','fontsize',10)

grid on
set(gca,'YTick',0:.5:2.5,'XTick',xticks,'fontsize',10)

datetick('x','dd','keepticks','keeplimits')
xlabel('Date in March and April 2019')

ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set

ana(1,1)=737496.3;ana(1,2)=737498.114722178;
ana(2,1)=737500.865973262;ana(2,2)=737502.278119836;
ana(3,1)=737506.125001883;ana(3,2)=737508.511042647;
ana(4,1)=737511.676198762;ana(4,2)=737513.307471529;

% p3=patch('Faces',[1 2 3 4],'Vertices',[ana(1,1) 0;ana(1,2) 0;ana(1,2) 1;ana(1,1) 1],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[ana(2,1) 0;ana(2,2) 0;ana(2,2) 2.5;ana(2,1) 2.5],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[ana(3,1) 0;ana(3,2) 0;ana(3,2) 2.5;ana(3,1) 2.5],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[ana(4,1) 0;ana(4,2) 0;ana(4,2) 2.5;ana(4,1) 2.5],'FaceColor',[139/256 0 0],'FaceAlpha',.3)

% title('Onshore wave height (CW4)','Fontsize',14);
t4=annotation('textbox','color','w','Position',[0.009 0.26 0.027 0.037], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(d)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle');%)%
l=legend([p1,p2,p3], {'CW1', 'CW4', '$\frac{CW4}{CW1}$'}, ...
    'location', 'eastoutside', 'orientation', 'vertical', 'Interpreter', 'latex','fontsize',10)
%%
s1.Position = [0.1    0.79    0.75    0.17];
s2.Position = [0.1   0.52    0.75    0.17];
s3.Position = [0.1    0.31   0.75    0.17];
s4.Position = [0.1    0.09    0.75    0.17];

l.Position = [0.8563    0.1685    0.1379    0.0912];
col.Position = [0.9059    0.31    0.032    0.38];
h1.Position = [0.9049    0.7907    0.0320    0.1702];
%%
print(f2,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\Figure2_off_on_spectra.tiff','-dtiff','-r300')