%% Load measured data
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m\AWAC_6734.mat')
% load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_8m\AWAC_cur_wav.mat')
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Aquadopps\AQD2985\AQD2985.mat')
% load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\Sawhorse_instruments.mat')
load ('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\AQP11266_dir.mat')
% meteo = xlsread('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\RDW47_YEARLY_PROCESSED\wind.xlsx', 'D1:G56');
% date = meteo(:,2); windSp = meteo(:,4); windDir=meteo(:,1); 
%% Load modeled data
for i=1:10
    ST = sprintf('ST%d',i);
    results = load (['C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\model\direction_withvegetation\rbr',...
        num2str(i) '_spc.dat']);
    RBR_Model.(ST).time_wave = datenum(num2str(results(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.(ST).Hsig = results(:,2);
    RBR_Model.(ST).Hrms = results(:,2)./sqrt(2);
    RBR_Model.(ST).Tp = results(:,3);
    RBR_Model.(ST).MeanDir = results(:,4);
    RBR_Model.(ST).PeakDir = results(:,5);
    clear results
end

results = load (['C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\model\direction_withvegetation\aqd_spc.dat']);
RBR_Model.aqd_2985.time_wave = datenum(num2str(results(:,1),'%f'),'yyyymmdd.HHMMSS');
RBR_Model.aqd_2985.Hsig = results(:,2);
RBR_Model.aqd_2985.Hrms = results(:,2)./sqrt(2);
RBR_Model.aqd_2985.Tp = results(:,3);
RBR_Model.aqd_2985.MeanDir = results(:,4);
RBR_Model.aqd_2985.PeakDir = results(:,5);
clear results
%% Prepare Aquadopp data
samples_fft = 1024;
% Read and create waves
AQP11266.Wave.time = AQP11266.TIME_RAW(1:samples_fft:end); AQP11266.Wave.time(end) = [];
AQP11266.time_ft = AQP11266.TIME_RAW(1:samples_fft/2:end); %AQP11266.time_ft(end) = [];

heights = AQP11266.z(2);

% Waves
lf=.04;                 % Hz - low frequency cutoff
maxfac=200;             % Maximum value of factor scaling pressure to waves
minspec=0.01;           % m^2/Hz - minimum spectral level for computing direction and spreading
Ndir=0;                 %deg - direction offset (includes compass error and misalignment of cable probe relative to case
                        % the offset for the Aquadopp Profiler is 0

parms=[lf maxfac minspec Ndir];
dt=1;
hp=1.09;
hv=0;

AQP11266.Wave.Su=[];AQP11266.Wave.Sp=[];AQP11266.Wave.Dir=[];AQP11266.Wave.Spread=[];AQP11266.Wave.F=[];AQP11266.Wave.dF=[];

ID.layout = zeros(3,3);
ID.datatypes = {'velx' 'vely' 'pres'};
ID.fs = 1;
SM.freqs = AQP11266.Wave.f_ADP(1,:);
SM.dirs = 0:4:360;
SM.S = zeros(129,91);
SM.xaxisdir = 0;
SM.funit={'Hz'};
SM.dunit={'naut'};

EP=[];

for i=1:floor((size((AQP11266.PRESSURE_RAW),1))/samples_fft)
        
        if sum (isnan(nanmean (AQP11266.Despike_vel_x(i*samples_fft - (samples_fft - 1):i*samples_fft,1:5),2))) < 200
         
        auxu = nanmean (AQP11266.Despike_vel_x(i*samples_fft - (samples_fft - 1):i*samples_fft,1:5),2)';
        if isnan(auxu(1))==1; auxu(1)=0; end; if isnan(auxu(end))==1; auxu(end)=0; end;
        X = ~isnan(auxu);
        Y = cumsum(X-diff([1 X])/2); vu = interp1(1:nnz(X),auxu(X),Y); clear auxu X Y
        
        auxv = nanmean (AQP11266.Despike_vel_y(i*samples_fft - (samples_fft - 1):i*samples_fft,1:5),2)';
        if isnan(auxv(1))==1; auxv(1)=0; end; if isnan(auxv(end))==1; auxv(end)=0; end;
        X = ~isnan(auxv);
        Y = cumsum(X-diff([1 X])/2); vv = interp1(1:nnz(X),auxv(X),Y); clear auxv X Y
        
        auxp = AQP11266.PRESSURE_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft,1)';
        if isnan(auxp)==1; auxp(1)=0; end; if isnan(auxp(end))==1; auxp(end)=0; end;
        X = ~isnan(auxp);
        Y = cumsum(X-diff([1 X])/2); vp = interp1(1:nnz(X),auxp(X),Y); clear auxv X Y
        
        [AQP11266.Wave.Su(i,:),AQP11266.Wave.Sp(i,:),AQP11266.Wave.Dir(i,:),AQP11266.Wave.Spread,AQP11266.Wave.F,AQP11266.Wave.dF] = wds(vu',-vv',vp',dt,100,hp,hv,parms);
        ID.data = [vu' vv' vp'];
%         ID.depth = AQP11266.meanDepth_15min(i);
        ID.depth = 2;
        [SMout(i),EP]=dirspec(ID,SM,EP, 'MESSAGE',0);
        aux=SMout(i).S;aux(1:10,:)=nan;aux(50:end,:)=nan;
        maximum=max(nanmax(aux));
        [l, c]=find(SMout(i).S==maximum);
        AQP11266.Wave.DirTp(i)=SMout(i).dirs(c(1));
        clear l
        else
        AQP11266.Wave.Su(i,1:41)=nan;
        AQP11266.Wave.Sp(i,1:41)=nan;
        AQP11266.Wave.Dir(i,1:41)=nan;
        AQP11266.Wave.Spread(i,1:41)=nan;
        AQP11266.Wave.F(i,1:41)=nan;
        AQP11266.Wave.dF(i,1:41)=nan;
        AQP11266.Wave.DirTp(i)=nan;
%        SMout(i)=[];
%        EP(i)=[];
%        Dp(i)=[];

        end
end
   
% AQP11266.Wave.Dir = wrapTo360 (AQP11266.Wave.Dir);

% for i=1:1517
%     % for i=1:floor((size((AQP11266.PRESSURE_RAW),1))/samples_fft)
%     if isempty (SMout(i).S)~=1;
%        
%     else
%         AQP11266.Wave.DirTp(i)=nan;
%     end
% end
%% Prepare offshore Aquadopp for directional spectra

local_aqd='C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Aquadopps\AQD2985';
was_aqd = load(fullfile(local_aqd,'ADP0101.was'));
wds_aqd = load(fullfile(local_aqd,'ADP0101.wds'));

nFreqs=98;nDirs=90;nSamples=635;
E_offshore=zeros(nFreqs,nDirs,nSamples);
w.was=was_aqd(2:end,:);
w.wds = permute(reshape(wds_aqd', nDirs, nFreqs, nSamples), [ 2 1 3]);

for tt=1:nSamples
    for f=1:nFreqs;
        for theta=1:nDirs;
            E_offshore(f,theta,tt) = (w.was(tt,f) * w.wds(f,theta,tt))/4; %divide por 4 para distribuir a energia na direção de cada bin
        end      
    end
end
%% Plot Mean Dir TS Aquadopp and AWACs modeled and measured
f7 = figure (7);
f7.Position = [545.0000  156.2000  636.0000  558.4000];
map=load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\bathycolormap.mat');
xticks=floor(AWAC_6734.wave.time(1)):1:AWAC_6734.wave.time(end);

s1=subplot (2,2,1)
p1=plot(AWAC_6734.wave.time,AWAC_6734.wave.MeanDir,'color','#016FB9','linewidth',1.5);
hold on
p2=plot(AQD2985.wav.time,AQD2985.wav.DirTp,'color','#FF9505','linewidth',1.2);
p3=plot(RBR_Model.aqd_2985.time_wave, RBR_Model.aqd_2985.MeanDir,'color','#66101F','linestyle','--','linewidth',1.2);
xlim([AQD2985.wav.time(1)+2 AQD2985.wav.time(end)]);
ylim([50 240])
grid on
ylabel('[°]')
set(gca,'XTick',xticks,'xticklabel',[], 'ytick',[0:45:360])
t=annotation('textbox','color','w','Position',[0.03 0.93 0.04 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(a)','FitBoxToText','off','fontsize',10,...
    'horizontalalignment','center','verticalalignment','middle');%)%;

% s2=subplot (3,2,3)
dir_swell=mod(AQP11266.Wave.Dir(:,15)-100,360);
dir_swell(dir_swell<=120 | dir_swell>=160)=nan;

for i=1:length(SMout)
    if isempty(SMout(i).S)~=1
    [row, col]=max(nansum(SMout(i).S));
    dir_swell_aqd(i)=SMout(i).dirs(col);
    total(i)=sum(nansum(SMout(i).S));
    end
end
    
p4=plot(AQP11266.Wave.time(1:1517),movmean(mod(dir_swell_aqd-250,360),5,'omitnan'),'color','#EC4E20','linewidth',1.5);

hold on
aux=RBR_Model.ST9.MeanDir;aux(aux>300)=nan;
p5=plot(RBR_Model.aqd_2985.time_wave, aux,'color','#66101F','linestyle','--','linewidth',1.2)
xlim([AQD2985.wav.time(1)+2 AQD2985.wav.time(end)])
ylim([40 250])
grid on
title (['Mean wave direction'],'fontsize',14)
ylabel('[°]')
xlabel('Date')
set(gca,'XTick',xticks,'ytick',[0:45:360],'fontsize',10)
datetick('x','dd','keepticks','keeplimits')
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xlabel('Date in March and April 2019')

ana(2,1)=737500.865973262;ana(2,2)=737502.278119836;
ana(3,1)=737506.125001883;ana(3,2)=737508.511042647;
ana(4,1)=737511.676198762;ana(4,2)=737513.307471529;

% p3=patch('Faces',[1 2 3 4],'Vertices',[ana(1,1) 0;ana(1,2) 0;ana(1,2) 2;ana(1,1) 2],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[ana(2,1) 0;ana(2,2) 0;ana(2,2) 250;ana(2,1) 250],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[ana(3,1) 0;ana(3,2) 0;ana(3,2) 250;ana(3,1) 250],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[ana(4,1) 0;ana(4,2) 0;ana(4,2) 250;ana(4,1) 250],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
%%
s3=subplot(2,2,3);
% E_offshore(25:end,:,:)=nan;
tick = pcolor_spectral(was_aqd(1,1:23), 0:4:356, nanmean(E_offshore(1:23,:,:),3)');
colormap(s3,flip(colorcet('L16', 'N', 128)));
col=colorbar;
xlim([-0.27 0.27]);ylim([-0.27 0.27]);
ylabel(col,'Energy [m^2/Hz/°]','Fontsize',10);
caxis([0 0.0068])
t=annotation('textbox','color','w','Position',[0.03 0.35 0.04 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(b)','FitBoxToText','off','fontsize',10,'horizontalalignment',...
    'center','verticalalignment','middle');%)%;

s4=subplot(2,2,4);
tick_onshore = pcolor_spectral(SMout(30).freqs(7:65), SMout(30).dirs-250, (SMout(100).S(7:65,:)'));
colormap(s4,flip(map.bathycolormap))
caxis([0.0000 0.0007]);
xlim([-0.27 0.27]);ylim([-0.27 0.27]);
col2=colorbar;
ylabel(col2,'Energy [m^2/Hz/°]','Fontsize',10);
t=annotation('textbox','color','w','Position',[0.5 0.35 0.04 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(c)','FitBoxToText','off','fontsize',10,'horizontalalignment','center',...
    'verticalalignment','middle');%)%;

%%
l=legend([p1 p2 p4 p3], {['CW1 (h \approx 15 m)'],['CW2 (h \approx 8 m)'],...
    ['CW4 (h \approx 2 m)'], ['Modelled']},'orientation','horizontal','fontsize',10);

s1.Position = [0.1300    0.5860    0.7939    0.3500];
% s2.Position = [0.1300    0.5706    0.7930    0.1500];
s3.Position =[0.0760    0.0112    0.34    0.39];
s4.Position =[ 0.5356    0.0112    0.34    0.39];
l.Position=[0.0966 0.4516 0.8474 0.0405];
t=annotation('textbox','color','k','Position',[0.682 0.385 0.04 0.04], 'BackgroundColor','none','EdgeColor','none','String','N','fontsize',10,'horizontalalignment','center',...
    'verticalalignment','top');%)%;
t=annotation('textbox','color','k','Position',[0.2208 0.385 0.04 0.04], 'BackgroundColor','none','EdgeColor','none','String','N','fontsize',10,'horizontalalignment','center',...
    'verticalalignment','top');%)%;
t=annotation('textbox','color','k','Position',[0.682 0.0 0.04 0.04], 'BackgroundColor','none','EdgeColor','none','String','S','fontsize',10,'horizontalalignment','center',...
    'verticalalignment','top');%)%;
t=annotation('textbox','color','k','Position',[0.2208 0.0 0.04 0.04], 'BackgroundColor','none','EdgeColor','none','String','S','fontsize',10,'horizontalalignment','center',...
    'verticalalignment','top');%)%;
t=annotation('textbox','color','k','Position',[0.06 0.189 0.04 0.04], 'BackgroundColor','none','EdgeColor','none','String','W','fontsize',10,'horizontalalignment','center',...
    'verticalalignment','top');%)%;
t=annotation('textbox','color','k','Position',[0.515 0.189 0.04 0.04], 'BackgroundColor','none','EdgeColor','none','String','W','fontsize',10,'horizontalalignment','center',...
    'verticalalignment','top');%)%;
t=annotation('textbox','color','k','Position',[0.39 0.189 0.04 0.04], 'BackgroundColor','none','EdgeColor','none','String','E','fontsize',10,'horizontalalignment','center',...
    'verticalalignment','top');%)%;
t=annotation('textbox','color','k','Position',[0.855 0.189 0.04 0.04], 'BackgroundColor','none','EdgeColor','none','String','E','fontsize',10,'horizontalalignment','center',...
    'verticalalignment','top');%)%;

%%
print(f7,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\Figure3_modeledDir.tiff','-dtiff','-r300')

