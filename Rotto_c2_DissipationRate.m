load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m\AWAC_6734.mat')
load 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd.mat'

%% Plot
xticks=RBR.ST1.time_wave(1):1:RBR.ST1.time_wave(end);
f9=figure(9)
% s1=subplot(2,2,1)
% 
% % p1=plot(RBR.ST1.time_wave,(RBR.ST2.E-RBR.ST3.E)./RBR.ST2.distance,'color','#016FB9','linewidth',1.5);
% p1=plot(RBR.ST1.time_wave,(RBR.ST2.delta)./(RBR.ST2.Hm0),'color','#016FB9','linewidth',1.5);
% hold on
% % p2=plot(RBR.ST1.time_wave,(RBR.ST8.E-RBR.ST9.E)./RBR.ST8.distance,'color','#FF9505','linewidth',1.5);
% p2=plot(RBR.ST1.time_wave,(RBR.ST8.delta)./(RBR.ST8.Hm0),'color','#FF9505','linewidth',1.5);
% xlim([RBR.ST1.time_wave(1)+2 RBR.ST1.time_wave(end)]);
% ylim([-2 22])
% grid on
% title (['Dissipation rate'])
% set(gca,'fontsize',14)
% ylabel('[W/m^2]')
% set(gca,'XTick',xticks,'ytick',[0:5:25],'fontsize',12)
% datetick('x','dd','keepticks','keeplimits')
% ax = gca;
% labels = string(ax.XAxis.TickLabels); % extract
% labels(2:2:end) = nan; % remove every other one
% ax.XAxis.TickLabels = labels; % set
% xlabel('Date in March and April 2019')
% legend('ST2-ST3','ST8-ST9','location','Northwest')
% t=annotation('textbox','color','w','Position',[0.03 0.9 0.027 0.07], 'BackgroundColor',[100 147 167]./255,...
%     'EdgeColor',[100 147 167]./255,'String','A','FitBoxToText','off','fontsize',14,'horizontalalignment','center','verticalalignment','top');%)%;
% 
% s2=subplot(2,2,3)
% s1=subplot(1,2,1)
% hrms_offshore=RBR.ST1.Hrms;
% dissipation_offshore = (RBR.ST2.E-RBR.ST3.E)./RBR.ST2.distance;
hrms_offshore=RBR.ST1.Hrms;
dissipation_offshore = RBR.ST7.delta;dissipation_offshore(RBR.ST7.ubr_LWT<0.1)=nan;

scatter(hrms_offshore,dissipation_offshore,5,'MarkerEdgeColor','#016FB9','MarkerFaceColor','#016FB9')
% hold on
% plot(hrms_offshore,1.127.*hrms_offshore.^2,'color','#66101F','linewidth',2)
% xlabel('H_{rms} at 15m deep')
xlabel('H_{rms,0}')
grid on; box on
ylabel ([{'\epsilon'}],'fontsize',16)
% ylim([0 4])
% t=annotation('textbox','color','w','Position',[0.03 0.35 0.027 0.07], 'BackgroundColor',[100 147 167]./255,...
%     'EdgeColor',[100 147 167]./255,'String','A','FitBoxToText','off','fontsize',14,'horizontalalignment','center','verticalalignment','top');%)%;
title (['P7 (-' num2str((round(nanmean(RBR.ST7.meandepth*10))/10)) ' m deep)'])
% set(gca,'XTick',0:0.5:2)

% s3=subplot(2,2,4)
% s2=subplot(1,2,2)
% % dissipation_onshore = (RBR.ST8.E-RBR.ST9.E)./RBR.ST8.distance;
% scatter(RBR.ST8.Hm0,RBR.ST8.delta,5,'MarkerEdgeColor','#FF9505','MarkerFaceColor','#FF9505')
% % xlabel('H_{rms} at 15m deep')
% xlabel('H_{rms}')
% grid on; box on
% ylabel ([{'Dissipation rate'},{'~ 3 m deep [W/m^2]'}])
% ylim([0 1])
% t=annotation('textbox','color','w','Position',[0.52 0.35 0.03 0.07], 'BackgroundColor',[100 147 167]./255,...
%     'EdgeColor',[100 147 167]./255,'String','B','FitBoxToText','off','fontsize',14,'horizontalalignment','center','verticalalignment','top');%)%;
% set(gca,'XTick',0:0.5:2)
% set(gca,'xtick', 1:8,'xticklabel', {'10','9','8','7','6','5','4','3'})
% xlabel('Station')
%%
% s1.Position=[ 0.1300    0.5571    0.7760    0.3679];
% s2.Position=[0.1264    0.1053    0.2600    0.2600];
% s3.Position=[0.6477    0.1053    0.2600    0.2600];

%%
print(gcf,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\Figure10_DissipationrateHrms_alt.tiff','-dtiff','-r300')
