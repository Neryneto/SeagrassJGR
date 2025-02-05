% load 'C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_RBRs_waves.mat'
load 'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd.mat'
%%
free_stream_col = '#AAC0AA';
top_canopy_col = 'k';
%%
[cxy2_5,fc2_5] = mscohere(detrend(RBR.ST9.data.depth),detrend(RBR.ST7.data.depth),hamming(2048),1024,RBR.ST10.f(:,1),2);
[cxy2_10,fc2_10] = mscohere(detrend(RBR.ST9.data.depth),detrend(RBR.ST1.data.depth),hamming(2048),1024,RBR.ST10.f(:,1),2);
%% Swell and Sea conditions / cross-correlation
f16=figure(16)
f16.Position = [0 0 845.6 707.2];
swell=1652;
sea=99;

for j=1:length(RBR.ST2.Syy)
    cumsum_st(:,j)=cumsum(RBR.ST2.Syy(:,j).*100)./nansum(RBR.ST9.Syy(:,j));
    if cumsum_st(206)<80
        flag_swell(j)=logical(1);
    else
        flag_swell(j)=logical(0);
    end
end
s1=subplot(2,2,1);
% plot(RBR.ST6.Hrms)
loglog(RBR.ST6.f(:,1),nanmean(RBR.ST1.Syy,2),'linewidth',1.7,'color',top_canopy_col)
hold on
loglog(RBR.ST6.f(:,1),nanmean(RBR.ST9.Syy,2),'linewidth',1.7,'color',free_stream_col)

xlim([0.001 1])
ylim([5*10^-5 0.5])
grid on
ylabel('E [m^2/Hz]')
xlabel('f [Hz]')
% title('Swell')
legend('P1','P9','location','northwest')
set(gca,'fontsize',12)
t=annotation('textbox','color','w','Position',[0.03 0.9 0.027 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(a)','FitBoxToText','off','fontsize',12,'horizontalalignment','center','verticalalignment','middle');%)%;
%%
s2=subplot(2,2,2)
% plot(RBR.ST6.Hrms)
% loglog(RBR.ST6.f(:,1),RBR.ST6.Syy(:,sea),'linewidth',1.7,'color',top_canopy_col)
% hold on
% loglog(RBR.ST6.f(:,1),RBR.ST10.Syy(:,sea),'linewidth',1.7,'color',free_stream_col)
% xlim([0 0.5])
% ylim([2*10^-5 0.2])
% grid on
% ylabel('E (m^2/Hz')
% xlabel('f (Hz)')
% title('Sea')
% set(gca,'fontsize',12)
t=annotation('textbox','color','w','Position',[0.47 0.9 0.027 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(b)','FitBoxToText','off','fontsize',12,'horizontalalignment','center','verticalalignment','middle');%)%;

% s3=subplot(2,2,3)
plot(1./fc2_5,cxy2_5,'linewidth',1.7,'color',top_canopy_col)
hold on
plot(1./fc2_10,cxy2_10,'linewidth',1.7,'color',free_stream_col)
xlim([0 25])
ylabel('Magnitude-squared coherence [ ]')
xlabel('Period [s]')
legend('P9-P7', 'P9-P1')
grid on
ylim([0 0.63])
set(gca,'fontsize',12)

%%
s3=subplot(2,2,3)

[l c]=find(RBR.ST6.time_wave>=737500 & RBR.ST6.time_wave<=737502);
[l2 c2]=find(RBR.ST6.time_wave>=737506 & RBR.ST6.time_wave<=737508.5);
[l3 c3]=find(RBR.ST6.time_wave>=737511 & RBR.ST6.time_wave<=737513);
% l_total=[c';c2'];
l_total=[c';c2'];

a1=loglog(RBR.ST6.f(:,l),nanmean(RBR.ST6.Syy(:,swell),2),'linewidth',1.7,'color',top_canopy_col);
hold on
a2=loglog(RBR.ST6.f(:,l3),nanmean(RBR.ST9.Syy(:,870),2),'linewidth',1.7,'color',free_stream_col);
a3=loglog(RBR.ST6.f(:,l3),nanmean(RBR.ST9.Syy(:,c3),2),'linewidth',1.7,'color','#FEB603');


% plot(RBR.ST6.f(:,1),cumsum(RBR.ST6.Syy(:,swell).*100)./nansum(RBR.ST6.Syy(:,swell)),'linewidth',1.7,'color',top_canopy_col)
% hold on
% plot(RBR.ST6.f(:,1),cumsum(RBR.ST6.Syy(:,sea).*100)./nansum(RBR.ST6.Syy(:,sea)),'linewidth',1.7,'color',free_stream_col)
set(gca,'XTick',[10^-5 10^-4 10^-3 10^-2 10^-1 10^-0]);
grid on
legend([a1(1),a2(1),a3(1)],{'Sea state 1','Sea state 2', 'Sea state 3'},'location','southwest')
ylabel('Cumulative Energy [%]')
xlabel('f [Hz]')
xlim([0.001 1])
ylim([5*10^-6 0.8])
set(gca,'fontsize',12)
% t=annotation('textbox','color','w','Position',[0.5 0.43 0.027 0.04], 'BackgroundColor',[100 147 167]./255,...
%     'EdgeColor',[100 147 167]./255,'String','D','FitBoxToText','off','fontsize',12,'horizontalalignment','center','verticalalignment','middle');%)%;
xlim([0 0.45])
t=annotation('textbox','color','w','Position',[0.2192 0.4093 0.0270 0.0400], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(c)','FitBoxToText','off','fontsize',12,'horizontalalignment','center','verticalalignment','middle');%)%;


% plot(RBR.ST6.f(:,1),cumsum(nanmean(RBR.ST6.Syy(:,flag),2).*100)./nansum(nanmean(RBR.ST6.Syy(:,flag),2)),'linewidth',1.7,'color',top_canopy_col)
% hold on
% plot(RBR.ST6.f(:,1),cumsum(nanmean(RBR.ST9.Syy(:,flag),2).*100)./nansum(nanmean(RBR.ST9.Syy(:,flag),2)),'linewidth',1.7,'color',free_stream_col)
% 
% [l c]=find(RBR.ST6.time_wave>=737500 & RBR.ST6.time_wave<=737502);
% [l2 c2]=find(RBR.ST6.time_wave>=737506 & RBR.ST6.time_wave<=737508.5);
% % [l3 c3]=find(RBR.ST6.time_wave>=737511 & RBR.ST6.time_wave<=737513);
% % l_total=[c';c2'];
% l_total=[c';c2'];
% 
% loglog(RBR.ST6.f(:,l),nanmean(RBR.ST6.Syy(:,l2),2),'linewidth',1.7,'color',top_canopy_col)
% hold on
% loglog(RBR.ST6.f(:,l3),nanmean(RBR.ST9.Syy(:,c3),2),'linewidth',1.7,'color',free_stream_col)

s3.Position=[0.3350    0.1109    0.3338    0.3412];
%%
print(f16,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\FigureS2_coherence_spectra.tiff','-dtiff','-r300')

