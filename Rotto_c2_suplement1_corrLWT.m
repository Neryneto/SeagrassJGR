load 'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_RBRs_waves.mat'
load('C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments')
load('C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Quadpod_instruments.mat')
load ('C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\Sawhorse_instruments.mat')
%%
free_stream_col = [142, 202, 230]/256;
top_canopy_col = 'r';
bottom_col = '#023047';
% bottom_col = [2, 48, 71]/256;
%%
f2=figure (2)
set(f2,'Units','centimeters','OuterPosition',[5 1 18 15])
s1=subplot(2,1,1)
AQD11898.ubr(AQD11898.ubr<0.1)=nan;
a3=area(AQD11898.time_ft,nanmean(AQD11898.ubr(:,1),2));
a3.FaceColor = free_stream_col;
a3.LineWidth = .2;

hold on
AQD11898.Wave.ubr_sherwood(AQD11898.Wave.ubr_sherwood>0.3)=nan;
plot(AQD11898.Wave.time,AQD11898.Wave.ubr_sherwood,'color',top_canopy_col,'linewidth',1.5) %verde

for i=1:122
    plot(Vectrino(i).InitialTime_Matlab,nanmean(Vectrino(i).ubr(10:15)),'.','color',...
        bottom_col,'markersize',15)
    hold on
end

xlim([AQD11898.time_ft(1) AQD11898.time_ft(end)-.75]) 
xticks=floor(AQD11898.time_ft(1)):1:ceil(AQD11898.time_ft(end));
xlabel ('Days in March 2019')

set(gca,'XTick',xticks,'YTick',0:0.1:0.4,'fontsize',10);
datetick('x','dd','keepticks','keeplimits')

% plot(AQP11266.time_15min,AQP11266.Wave.ubr_sherwood) %verde
h = ylabel('$u_{R} [m/s]$','interpreter','latex','fontsize',14)
% set(h,'Interpreter','latex')

xlim([737494.55 737498.5])
ylim([0 .3]) 
grid on
l=legend('$u_{R, \infty}$','$u_{R,LWT}$', '$\hat{u}_{R}$','location', ...
    'northwest','interpreter','latex','fontsize',12);

t=annotation('textbox','color','w','Position',[0.03 0.9 0.027 0.04], ...
    'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(a)','FitBoxToText','off', ...
    'fontsize',10,'horizontalalignment','center','verticalalignment', ...
    'middle');%)%;
%%
s2=subplot(2,1,2)

scatter(AQD11898.Wave.ubr_sherwood,nanmean(AQD11898.ubr(:,1),2),5,'markerfacecolor',...
        bottom_col,'markeredgecolor',bottom_col)
hold on
plot([0 .5],[0 .5],'-r')

xlim([0 0.3]) 
ylim([0 0.3]) 
grid on
ylabel('$u_{R, \infty} [m/s]$','interpreter','latex','fontsize',14)
xlabel('$u_{R, LWT} [m/s]$','interpreter','latex','fontsize',14)
set(gca,'fontsize',10);
box on
t=annotation('textbox','color','w','Position',[0.03 0.45 0.027 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(b)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle');%)%;
sgtitle('$u_{R}$ velocities at CW4','interpreter','latex','fontsize',16)
%%
s2.Position = ([0.301    0.12    0.4015    0.33]);
l.Position = ([0.715    0.296    0.1887    0.1535]);
%%
print(gcf,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\FigureS1_LWT_Vectrino.tiff','-dtiff','-r300')
