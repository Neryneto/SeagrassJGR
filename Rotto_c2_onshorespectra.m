load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m\AWAC_6734.mat')
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd.mat')
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Quadpod_instruments.mat')
%%
caminho='C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m';
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

load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments.mat');

clear caminho
caminho='C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Aquadopps\AQD2985';
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
% tempo=tempo(1:1782);
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


eli_antes=0;eli_dep=0;
pwrFreq_aqd(pwrFreq_aqd == -9.000000) = NaN; 
pwrFreqDir_aqd(pwrFreqDir_aqd == -9.000000) = NaN;
waveData_aqd.Frequency = pwrFreq_aqd(1,:)';
waveData_aqd.pwrSpectrum           = pwrFreq_aqd(eli_antes+2:end-eli_dep,:); 
waveData_aqd.fullSpectrumDirection = 4:4:360;
nFreqs   = size(pwrDir,2);
% waveData_aqd.fullSpectrum          = pwrFreqDir((eli_antes*nFreqs)+1:end-eli_dep*nFreqs,:);
waveData_aqd.fullSpectrum          = pwrFreqDir_aqd((eli_antes*nFreqs)+1:end-eli_dep*nFreqs,:);
% waveData_aqd.fullSpectrum (end,:) = [];
waveData_aqd.Frequency=waveData_aqd.Frequency(1:nFreqs,:);
nDirs    = length(waveData_aqd.fullSpectrumDirection);
nSamples_aqd = length(waveData_aqd.fullSpectrum) / nFreqs;
aux=waveData_aqd.fullSpectrum;
% aux=aux(:,1:2304);b=reshape(aux, nDirs, nFreqs, nSamples);
waveData_aqd.fullSpectrum = permute(reshape(aux, nDirs, nFreqs, round(nSamples_aqd)), [ 2 1 3]);
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

E=zeros(nFreqs,nDirs,nSamples_aqd);
for t=1:nSamples_aqd
    for f=1:nFreqs;
        for theta=1:nDirs;
            E(f,theta,t) = waveData_aqd.pwrSpectrum(t,f) * waveData_aqd.fullSpectrum(f,theta,t);
        end
    end
end
%% Prepare RBR measured data

Hsig_rbr_interp_aux=[];
for i=1:10
    
    STi = sprintf('ST%d',i);
    
    Hsig_rbr_interp_aux (:,i) = interp1(RBR.(STi).time_wave, RBR.(STi).Hsig_swell_out, AWAC_6734.wave.time)	;
    if i<10
        distances10 (i) = RBR.(STi).distance;
    end
    distance_cumulative_aux (i,2) = nanmean(RBR.(STi).meandepth);    
end
Hsig_rbr_interp=fliplr(Hsig_rbr_interp_aux);
Hsig_measured=[Hsig_rbr_interp AWAC_6734.wave.Hm0];
% distances10(10)=0;
distances_with_awac=[251 distances10];
distance_cumulative_with_awac(:,1) = [0 cumsum(fliplr(distances_with_awac(1:end)))];
distance_cumulative_with_awac(end,2)=15.1;
distance_cumulative_with_awac(1:10,2)=flip(distance_cumulative_aux(:,2));

samplingRateIncrease = 3;

newXSamplePoints_awac = [linspace(min(distance_cumulative_with_awac(:,1)), max(distance_cumulative_with_awac(:,1)),...
    length(distance_cumulative_with_awac(:,1)) * samplingRateIncrease)];
smoothedY_awac = [spline(distance_cumulative_with_awac(:,1),...
    distance_cumulative_with_awac(:,2), newXSamplePoints_awac)];


%% Average Hsig each station
f4=figure(4)
set(f4,'Units','centimeters','Position',[12.8852    0   18   15])

fc='#016FB9';

s1=subplot(3,1,1)

yyaxis right
a=area(newXSamplePoints_awac,(-smoothedY_awac),-15.5,'Facecolor',[.5 .5 .5]);
a.EdgeColor=[.5 .5 .5];
a.FaceAlpha=0.3
hold on
% xlim([0 max(newXSamplePoints)])
xlim([0 1650])
% xlim([0 max(newXSamplePoints_awac)])
ylim([-15.3 -0])
ylabel('Depth [m]')

yyaxis left
%not normalized
errorbar(distance_cumulative_with_awac(:,1),nanmean(Hsig_measured),nanstd(Hsig_measured),'o--','linewidth',1.4,...
   'markeredgecolor',fc,'markerfacecolor',fc)
ax=gca;
ax.YAxis(2).Color = [0.5 0.5 0.5];
ax.YAxis(1).Color = fc;

set(gca,'fontsize',10,'xtick',[0:200:1600])
ylabel('H_{s,i} [m]')
grid on
ylim([0 1.5])
pdense=plot([1000 1000],[0 1.5],'-','linewidth',2,'color',[139/256 0 0]);
text(1480,1.2,'CW1','fontsize',10)
for i=1:10
%     if i==1
%         text(distance_cumulative_with_awac(i)-20,.5,['P' num2str(11-i)],'fontsize',12)
%     else
    text(distance_cumulative_with_awac(i)-15,1.2,['P' num2str(11-i)],'fontsize',10)
    end
% end
xlabel('Distance offshore from P10 [m]')

%%
s2=subplot(3,1,2)
pcolor(AQD11898.time_15min,AQD11898.Wave.f_ADP(1,:),log(AQD11898.Wave.Spp_ADP)')
hold on
pcolor(RBR.ST2.time_wave,RBR.ST9.f_ss(:,1)',log(RBR.ST9.Syy))

shading interp;
colormap(s2,map1);
% caxis([0 .05])
caxis([-5 0]);
ylim([0.02 0.5]);
xlim([AQD11898.time_15min(1) RBR.ST2.time_wave(end)]);
yticks=0:.1:.5;
ylabel('Frequecy [Hz]')
xticks=737492:1:737519;
set(gca,'YTick',yticks,'XTick',xticks,'xticklabel',[],'fontsize',10)
col=colorbar;
title(col,{'log (energy)', '[m^2/Hz]'},'Fontsize',10);
title('Onshore energy spectra (CW4)','Fontsize',14)
grid on
set(gca,'layer','top')

%% 
s3=subplot(3,1,3)
color3='#FF9505';
awac_interp=interp1(AWAC_6734.wave.time,AWAC_6734.wave.Hm0,AQD11898.time_15min);
awac_interp_rbr=interp1(AWAC_6734.wave.time,AWAC_6734.wave.Hm0,RBR.ST2.time_wave);
ratio1=AQD11898.Wave.Hm0_ADP./awac_interp;ratio1(ratio1>0.7)=nan;
ratio2=RBR.ST8.Hm0./awac_interp_rbr;ratio2(ratio2>0.7)=nan;

AQD11898.Wave.Hm0_ADP(AQD11898.Wave.Hm0_ADP>1 | AQD11898.Wave.Hm0_ADP<0)=nan;
area(AQD11898.time_15min,AQD11898.Wave.Hm0_ADP, 'Facecolor',color3)
hold on
p=area(RBR.ST2.time_wave,RBR.ST9.Hm0, 'Facecolor',color3)
p2=plot(AQD11898.time_15min,ratio1,'r','linewidth',1.5)
plot(RBR.ST2.time_wave,ratio2,'r','linewidth',1.5)

xlim([AQD11898.time_15min(1) RBR.ST2.time_wave(end)]);
% xlim([AQD11898.time_ft(1) AQD11898.time_ft(end)])
ylim([0.02 1])
yh=ylabel({'H_{s,CW4} [m]';' '},'color','k')
set(yh,'Units','Normalized');
ylabPos = get(yh,'Position');
text(ylabPos(1)-.07,ylabPos(2),'H_{s,CW4}/H_{s,CW1} [ ]',  'Units','Normalized',  'HorizontalAlignment','center',...
    'Vert','top', 'Rotation',90,'color','r','fontsize',10)

grid on
set(gca,'YTick',0:.25:1,'XTick',xticks,'xticklabel',[],'fontsize',10)

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
patch('Faces',[1 2 3 4],'Vertices',[ana(2,1) 0;ana(2,2) 0;ana(2,2) 1;ana(2,1) 1],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[ana(3,1) 0;ana(3,2) 0;ana(3,2) 1;ana(3,1) 1],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[ana(4,1) 0;ana(4,2) 0;ana(4,2) 1;ana(4,1) 1],'FaceColor',[139/256 0 0],'FaceAlpha',.3)

title('Onshore wave height (CW4)','Fontsize',14);
%%
s1.Position = [0.1353    0.75    0.7366    0.23];
s2.Position = [0.1353    0.38    0.7366    0.23];
s3.Position = [0.1353    0.09    0.7366    0.23]

t=annotation('textbox','color','w','Position',[0.01 0.94 0.027 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','A','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle')%)%
t=annotation('textbox','color','w','Position',[0.01 0.61 0.027 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','B','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle')%)%
t=annotation('textbox','color','w','Position',[0.01 0.31 0.027 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','C','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle')%)%

%%
print(f4,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\Figure4_onshorespectra.tiff','-dtiff','-r300')