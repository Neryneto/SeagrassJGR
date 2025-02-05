load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments')
load ('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\Sawhorse_instruments.mat')
load ('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd')
%%
lv = 0.43;
rhow=1026;
Cd=1.37;
Uc=nanmean(AQP11266.Mean_Current(:,18:30),2);
E=400000000;
t=0.3*10^-3;
% t=150*10^-6;
Fd=6*rhow*Cd*lv^3.*Uc.^2;
Fr=E*t^3;
Ca = Fd/Fr;

rhos=860;g=9.8;
% B=(12*(rhow-rhos)*g*lv^3)/(E*t^2);
Arms_infty = nanmean(AQP11266.ubr(:,1),2).*AQP11266.Wave.Tm02_ADP'./(2*pi);
L=lv./Arms_infty;

le=lv.*(Ca.*L).^(-1/4);
% le(le>0.4)=nan;
le2 = le;
le2(AQP11266.Wave.ubr_sherwood<0.17)=nan;

lambdaF_le=le*512*0.015;

%%
% E1=51*10^6;
% % L1=I./Arms_infty;
% 
% % Ca1=(12*rhow*Uc.^2*lv^3)/(E*t^3);
% % le2=lv.*((Ca1.*L).^(-1/4));
% 
% subplot(2,1,1)
% plot(AQP11266.time_15min,Ca.*L)
% hold on
% plot(AQP11266.time_15min,Ca1.*L)
% set(gca, 'YScale', 'log');
% legend('E=0.4 GPa','E=5.1 MPa')
% ylabel('CaL [ ]');xlabel('Date')
% set(gca,'XTick',xticks,'ytick',[0:.25:1],'fontsize',12)
% xticks=floor(AQP11266.time_15min(1))-1:1:AQP11266.time_15min(end);
% datetick('x','dd/mm','keepticks','keeplimits')
% grid on
% 
% subplot(2,1,2)
% plot(Ca.*L,le./lv,'.')
% hold on
% plot(Ca1.*L,le2./lv,'.')
% plot([0:2500],(0.94+0.06).*([0:2500].^-(0.25-0.02)),'--k')
% plot([0:2500],(0.94+0.06).*([0:2500].^-(0.25+0.02)),'--k')
% 
% set(gca, 'YScale', 'log','XScale', 'log');
% legend('E=0.4 GPa','E=5.1 MPa','l_e/l=(0.94±0.06).*(CaL.^-(0.25±0.02))','fontsize',8, ...
%     'orientation','vertical','Interpreter','latex')
% ylabel('l_e/l [ ]');xlabel('CaL')
% grid on
% set(gca,'fontsize',12)
% print(gcf,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\FigureCaL.tiff','-dtiff','-r300')
%%
f12=figure(12)
set(f12,'Units','centimeters','OuterPosition',[5 9 18 10])

yyaxis left
a3=area(AQP11266.time_15min, movmean(le*514*0.015,6,'omitnan'),'facecolor','#5297AD','EdgeColor','#5297AD');
ylim([0 3])
ylabel('$$\lambda_{f,e} [\;]$$', 'Interpreter','Latex','fontsize',14) 
set(gca,'ytick',0:.5:3,'fontsize',12)
ax=gca;ax.YAxis(1).Color = '#5297AD';
% set(gca, 'SortMethod', 'depth')
grid on
box on

yyaxis right
ub=AQP11266.ubr(:,1);ub(ub<0.1)=nan;
hold on
set(gca,'layer','top')

Arms_infty(Arms_infty<.02 | Arms_infty>0.5)=nan;
Arms_infty2=Arms_infty;
Arms_infty2(Arms_infty2<0.15)=nan;
a1=plot(AQP11266.time_15min,movmean(Arms_infty,4,'omitnan'),'-','color','k','LineWidth',1.4)
ax=gca;ax.YAxis(2).Color = 'k';
ax.GridColor = 'k'
xlim([AQP11266.time_15min(1) AQP11266.time_15min(end)])

set(gca,'xtick',xticks)
% area(AQP11266.time_15min, movmean(le,6,'omitnan'),'linewidth',2);
datetick('x','dd','keepticks','keeplimits')
xlabel ('Days in March 2019')
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set

xlim([AQP11266.time_15min(1) AQP11266.time_15min(end)])
ana(1,1)=737496.3;ana(1,2)=737498.114722178;
ana(2,1)=737500.865973262;ana(2,2)=737502.278119836;
ana(3,1)=737506.125001883;ana(3,2)=737508.511042647;
ana(4,1)=737511.676198762;ana(4,2)=737513.307471529;

% p3=patch('Faces',[1 2 3 4],'Vertices',[ana(1,1) 0;ana(1,2) 0;ana(1,2) 1;ana(1,1) 1],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[ana(2,1) 0;ana(2,2) 0;ana(2,2) 3;ana(2,1) 3],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[ana(3,1) 0;ana(3,2) 0;ana(3,2) 3;ana(3,1) 3],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
patch('Faces',[1 2 3 4],'Vertices',[ana(4,1) 0;ana(4,2) 0;ana(4,2) 3;ana(4,1) 3],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
ylim([0 0.25])
ylabel('$$A_{\infty,R}  [m]$$','Interpreter','Latex','fontsize',14);
%%
% f12=figure(12)
% set(f12,'Units','centimeters','OuterPosition',[5 9 18 10])
% 
% s1=subplot(2,1,1)
% yyaxis left
% p1=plot(AQP11266.time_15min, movmean(le/lv,5,'omitnan'),'.-','color','k');
% xlim([AQP11266.time_15min(1) AQP11266.time_15min(end)])
% % datetick('x','dd','keepticks','keeplimits')
% % xlabel ('Days in March 2019')
% ylabel('$$l_{e,v}/l_v [\;]$$','Interpreter','Latex','fontsize',14);
% set(gca,'xticklabel',[],'fontsize',12)
% grid on
% box on
% ylim([0 1])
% ax=gca;ax.YAxis(1).Color = 'k';
% set(gca,'XTick',xticks,'ytick',[0:.25:1],'fontsize',12)
% xticks=floor(AQP11266.time_15min(1))-1:1:AQP11266.time_15min(end);
% box on
% yyaxis right
% ub=AQP11266.ubr(:,1);ub(ub<0.1)=nan;
% hold on
% set(gca,'layer','top')
% 
% Arms_infty(Arms_infty<.02 | Arms_infty>0.5)=nan;
% Arms_infty2=Arms_infty;
% Arms_infty2(Arms_infty2<0.15)=nan;
% a1=area(AQP11266.time_15min,movmean(Arms_infty,4,'omitnan'),'facecolor','#FEB603')
% ax=gca;ax.YAxis(2).Color = '#FEB603';
% ylabel('$$A_{\infty,R}  [m]$$','Interpreter','Latex','fontsize',14);
% set(gca, 'SortMethod', 'depth')
% set(gca,'XTick',xticks,'ytick',[0:.1:.3],'fontsize',12)
% ylim([0 0.25])
% %
% s2=subplot(2,1,2)
% a3=area(AQP11266.time_15min, movmean(le*514*0.015,6,'omitnan'),'facecolor',[142, 202, 230]/256);
% hold on
% % a2=plot(AQP11266.time_15min,movmean(le2*514*0.015,6,'omitnan'),'.','color','#023047')
% set(gca,'xtick',xticks)
% % area(AQP11266.time_15min, movmean(le,6,'omitnan'),'linewidth',2);
% datetick('x','dd','keepticks','keeplimits')
% xlabel ('Days in March 2019')
% ax = gca;
% labels = string(ax.XAxis.TickLabels); % extract
% labels(2:2:end) = nan; % remove every other one
% ax.XAxis.TickLabels = labels; % set
% 
% ylabel('$$\lambda_{f,e} [\;]$$', 'Interpreter','Latex','fontsize',14) 
% set(gca,'ytick',0:.5:3,'fontsize',12)
% grid on
% box on
% xlim([AQP11266.time_15min(1) AQP11266.time_15min(end)])
% leg=legend([a1 p1 a3],{'$$A_{\infty,R}$$','$$l_{e,v}/l_v$$', ...
%     '$$\lambda_{f,e}$$'},'Interpreter','Latex','fontsize',14);
% % legend ([pbot ptop pk],{'Near-bed','Top of canopy','D_{10}'})
% ana(1,1)=737496.3;ana(1,2)=737498.114722178;
% ana(2,1)=737500.865973262;ana(2,2)=737502.278119836;
% ana(3,1)=737506.125001883;ana(3,2)=737508.511042647;
% ana(4,1)=737511.676198762;ana(4,2)=737513.307471529;
% 
% % p3=patch('Faces',[1 2 3 4],'Vertices',[ana(1,1) 0;ana(1,2) 0;ana(1,2) 1;ana(1,1) 1],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
% patch('Faces',[1 2 3 4],'Vertices',[ana(2,1) 0;ana(2,2) 0;ana(2,2) 3;ana(2,1) 3],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
% patch('Faces',[1 2 3 4],'Vertices',[ana(3,1) 0;ana(3,2) 0;ana(3,2) 3;ana(3,1) 3],'FaceColor',[139/256 0 0],'FaceAlpha',.3)
% patch('Faces',[1 2 3 4],'Vertices',[ana(4,1) 0;ana(4,2) 0;ana(4,2) 3;ana(4,1) 3],'FaceColor',[139/256 0 0],'FaceAlpha',.3)

%%
% leg.Position = [0.8363    0.268    0.1496    0.2718];
% s1.Position = [0.1300    0.5838    0.6993    0.34];
% s2.Position = [0.1300    0.2    0.6993    0.34];
% 
% t=annotation('textbox','color','w','Position',[0.03 0.88 0.03 0.04], 'BackgroundColor',[100 147 167]./255,...
%     'EdgeColor',[100 147 167]./255,'String','(a)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle');%)%;
% 
% t=annotation('textbox','color','w','Position',[0.03 0.5 0.03 0.04], 'BackgroundColor',[100 147 167]./255,...
%     'EdgeColor',[100 147 167]./255,'String','(b)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle');%)%;

%%
% save ('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_lev.mat','le','lambdaF_le')
print(f12,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\Figure12_lev','-dtiff','-r900')
% l_ve = lv/(Ca*L)^1/3;