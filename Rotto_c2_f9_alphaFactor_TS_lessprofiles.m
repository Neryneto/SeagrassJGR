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
le_vectrino=[];

for ii=1:length(Vectrino);
    avg_vectrino = nanmean(Vectrino(ii).MeanVel(11:15));
    avg_vectrino_rms = nanmean(Vectrino(ii).ubr(11:15));
    
    [~,idx] = min(abs(AQD11898.time_ft-Vectrino(ii).InitialTime_Matlab));
    [~,idx_vector] = min(abs(ADV1670.time_ft-Vectrino(ii).InitialTime_Matlab));
    
    [~,idx_le1] = min(abs(AQD11898.time_ft-Vectrino(ii).InitialTime_Matlab));
    [~,idx_le] = min(abs(AQP11266.time_15min-Vectrino(ii).InitialTime_Matlab));
    
    le_vectrino(ii) = le(idx_le);
    le_vectrino1(ii) = le(idx_le1);
    
    ratio_vectrino_freestream_cur (ii) = avg_vectrino/AQD11898.freestream(idx);
    ratio_vectrino_freestream_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,30:35));

    ratio_aqd_freestream_cur (ii,:) = AQD11898.Mean_Current(idx,:)/nanmean(AQD11898.Mean_Current(idx,30:35));
    ratio_aqd_freestream_urms (ii,:) = AQD11898.ubr(idx,:)/nanmean(AQD11898.ubr(idx,30:35));

    
    if Vectrino(ii).BottomDistance<0.4 && Vectrino(ii).InitialTime_Matlab>737494.5
    % if Vectrino(ii).BottomDistance<0.15 && Vectrino(ii).InitialTime_Matlab>737494.5
              % violin_vectrino_vector_cur (ii) = avg_vectrino/ADV1670.Current_Mean(idx_vector);

        aux_frequencies = Vectrino(ii).fx;
        aux_frequencies = ADV1670.fx(1:15);
        Wavek =  wavek(aux_frequencies,ADV1670.press_15min(idx_vector)+.55);         
        phij (:,ii) =    (cosh (Wavek.*0.55))./(cosh (Wavek.*Vectrino(ii).BottomDistance));
        for i=1:length(aux_frequencies)-1
            avg_spec_bot (:,i) = nanmean(Vectrino(ii).Svel(Vectrino(ii).fx>aux_frequencies(i) & Vectrino(ii).fx<=aux_frequencies(i+1)));        
            spec_x (:,i) = nanmean(ADV1670.Suu(idx_vector,ADV1670.fx>aux_frequencies(i) & ADV1670.fx<=aux_frequencies(i+1)));
            spec_y (:,i) = nanmean(ADV1670.Svv(idx_vector,ADV1670.fx>aux_frequencies(i) & ADV1670.fx<=aux_frequencies(i+1)));
            avg_spec_top(:,i) = sqrt(2.*(spec_x(:,i)'+spec_y(:,i)').*0.0061);
        end

        alphaj (:,ii) = sqrt((phij(1:end-1,ii).^2).*avg_spec_bot'./avg_spec_top');

        % violin_vectrino_vector_urms (ii) = sqrt(nansum(alphaj (:,ii).^2.*(avg_spec_bot(:,i)'))./avg_spec_top(:,i)');
        violin_vectrino_vector_cur (ii) = avg_vectrino/ADV1670.Current_Mean(idx_vector);
        violin_vectrino_vector_urms (ii) = avg_vectrino_rms/ADV1670.ubr(idx_vector);
        
        violin_vectrino_aqd_cur (ii) = avg_vectrino/AQD11898.freestream(idx);
        violin_vectrino_aqd_urms (ii) = avg_vectrino_rms/nanmean(AQD11898.ubr(idx,30:35));
    else
        violin_vectrino_vector_cur (ii) = nan;
        violin_vectrino_vector_urms (ii) = nan;
        
        violin_vectrino_aqd_cur (ii) = nan;
        violin_vectrino_aqd_urms (ii) = nan;
        %         violin_vectrino_aqd_Arms (ii) = nan;
    end
end

%%
top_canopy_col = [221/256 96/256 49/256];
bot_canopy_col = [33/256 104/256 105/256];
bottom_col = [2/255 48/255 71/255];
free_stream_col = [255/255 183/255 3/255];
map_can = load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\bathycolormap.mat');

% map_bare = colorcet('L18', 'N', 64,'reverse',1);

f9=figure(9)
set(f9,'Units','centimeters','Position',[5 5 18 10])


s5=subplot(1,3,1)
hold on
for i=75:86
p1=plot(ratio_vectrino_freestream_cur(i),Vectrino(i).BottomDistance,'*',... 
         'markersize',7, 'color', bottom_col);
   p2=  plot(ratio_vectrino_freestream_urms(i),Vectrino(i).BottomDistance,'s',... 
         'markersize', 7, 'markerfacecolor', bottom_col, 'markeredgecolor',bottom_col);
end
p3=plot(nanmean(ratio_aqd_freestream_cur(75:86,:)),AQD11898.z,'*',... 
         'color',free_stream_col)
p4=plot(nanmean(ratio_aqd_freestream_urms(75:86,:)),AQD11898.z,'s',... 
         'color',free_stream_col)
p5=plot([0 1],[0.43 0.43],'--k')
p6=plot([0 1],[nanmean(le_vectrino(75:86)) nanmean(le_vectrino(75:86))],'-k')
set(gca,'xtick',[0.25:0.25:1],'yscale','log','ytick',[0.05 0.1 0.5 1],'yticklabel',[{'5','10','50','100'}], 'fontsize',10)  
ylim([0.04 1.5])
xlim([0 1])
xlabel('[ ]')
% title('Velocity profile')
grid on
box on
% text(0.55,2.05,'Velocity profiles','fontsize',16)
% text(0.55,-0.06,'Normalised velocity [ ]','fontsize',12)
ylabel('Elevation above bed [cm]','fontsize',12)

l1=legend([p1 p2 p3 p4 p5 p6],[{'$\hat{\overline{u}}/\overline{u}_{C1}$','$\hat{u}_{rms}/u_{C1,rms}$',...
           '$\overline{u}_{c}/\overline{u}_{C1}$','$u_{rms}/u_{C1,rms}$','$l_v$','$l_{v,e}$'}],'interpreter','latex','fontsize',12);
tA=annotation('textbox','color','w','Position',[0.009 0.26 0.027 0.037], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(a)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle');%)%
% h=sgtitle('Velocity profiles','fontsize',15)
       %%
s9=subplot(1,3,2)
violin_vectrino_vector_cur(violin_vectrino_vector_cur>0.7)=nan;
[counts,centers] = hist(violin_vectrino_vector_cur,0:0.05:0.8);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',flip(map_can.bathycolormap),'showMM',6,...
    'ylabel','Attenuation parameter (\alpha [  ])');
xlabel('$\alpha_{c} [ \ ]$','Interpreter','Latex','fontsize',12)
ylim([0 1])
set(gca, 'YTick',0:.1:1,'xtick',[], 'XTickLabel',[],'fontsize',10)
xlim([0.5 1.5])
grid on
box on

s10=subplot(1,3,3)
violin_vectrino_vector_urms(violin_vectrino_vector_urms>0.95)=nan;
% [counts,centers] = hist(sqrt(violin_vectrino_vector_urms),0:0.05:1);counts=counts*100/sum(counts);
[counts,centers] = hist(violin_vectrino_vector_urms,0:0.05:1);counts=counts*100/sum(counts);
[avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
    'divFactor',centers,'globalnorm',2,...
    'colormap',flip(map_can.bathycolormap),'showMM',6,'ylabel',' ');
ylim([0 1])
set(gca, 'YTick',0:.1:1,'yticklabel',[],'xtick',[], 'XTickLabel',[],'fontsize',12)
xlim([0.5 1.5])
xlabel ('$\alpha_{w,R} [ \ ]$','Interpreter','Latex','fontsize',12)
% caxis([0 maximo])
grid on
box on
% text(0.75,1.2,'Attenuation histograms','fontsize',15)


% s11=subplot(2,8,14)
% [counts,centers] = hist(violin_vectrino_aqd_cur,0:0.05:1);counts=counts*100/sum(counts);
% [avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
%     'divFactor',centers,'globalnorm',2,...
%     'colormap',flip(map_can.bathycolormap),'showMM',2,'ylabel',' ');
% xlabel ('$\hat{\overline{u}}/\overline{u}_{\infty} [ \ ]$','Interpreter','Latex','fontsize',12)
% ylim([0 1])
% set(gca, 'YTick',0:.1:1,'yticklabel',[],'xtick',[], 'XTickLabel',[],'fontsize',12)
% xlim([0.5 1.5])
% grid on
% box on
% % text(0.55,1.07,'(D)','fontsize',12)
% % text(0.55,.05,{['$$\alpha_c = $$ ' num2str(round(avg1{2}(1)*100)/100)]},'Interpreter', 'LaTeX','fontsize',12)
% %%
% s12=subplot(2,8,16)
% [counts,centers] = hist(violin_vectrino_aqd_urms,0:0.05:1);counts=counts*100/sum(counts);
% [avg1]=distributionPlot({[centers',counts']},'histOpt',99,...
%     'divFactor',centers,'globalnorm',2,...
%     'colormap',flip(map_can.bathycolormap),'showMM',2,'ylabel',' ');
% ylim([0 1])
% set(gca, 'YTick',0:.1:1,'yticklabel',[])
% xlim([0.5 1.5])
% set(gca,'xtick',[], 'XTickLabel',[])
% xlabel('$\hat{u}_{rms}/u_{rms,\infty} [ \ ]$','Interpreter','Latex','fontsize',14)
% grid on
% box on
ccol = colorbar ('southoutside');ccol.FontSize=10;
hL = ylabel(ccol,'Frequency [%]','fontsize',10);     
tB=annotation('textbox','color','w','Position',[0.009 0.2 0.027 0.037], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(b)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle');%)%

tC=annotation('textbox','color','w','Position',[0.009 0.26 0.027 0.037], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(c)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle');%)%
%%
% s1.Position = [0.08 0.45 0.2 0.49];
s5.Position = [0.0732 0.16 0.19 0.76];
% s7.Position = [0.5 0.45 0.2 0.49];
% s8.Position = [0.71 0.45 0.2 0.49];
s9.Position = [0.56    0.2523    0.21    0.6668];
s10.Position = [0.78    0.2523    0.21    0.6668];
% s11.Position = [0.5 0.04 0.2 0.25];
% s12.Position = [0.71 0.04 0.2 0.25];
l1.Position = [0.27    0.595    0.1971    0.3243];
tA.Position = [0.0732    0.9356    0.0270    0.0370];
tB.Position = [0.562    0.9356    0.0270    0.0370];
tC.Position = [0.78    0.9356    0.0270    0.0370];
ccol.Position = [0.6    0.1200    0.39    0.0553];

%%
print(f9,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\Figure9_alpha_TS.tiff','-dtiff','-r300')

