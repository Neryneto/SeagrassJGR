load('C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m\AWAC_6734.mat')
load 'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd.mat'
%% Calculate using Falter methodology
for k=1:9 %for each instrument
    ST1 = sprintf('ST%d',k);
    ST2 = sprintf('ST%d',k+1);
    
    ratioF_Fr(:,k) = ((Model.(ST2).F_interp./Model.(ST1).F_interp)-1)./-((RBR.(ST2).Fr./RBR.(ST1).Fr)-1);
    ratioF_cg(:,k) = RBR.(ST2).Cg./RBR.(ST1).Cg;
%     ratioF_Fr(:,k) = -((RBR.(ST2).Fr-((Model.(ST2).F_interp./Model.(ST1).F_interp).*RBR.(ST1).Fr)))/RBR.(ST1).distance;

% test(:,k)=(RBR.(ST2).Fr-RBR.(ST1).Fr);
end
ratioF_Fr(ratioF_Fr<-1.86 | ratioF_Fr>3)=nan;
aux=ratioF_Fr(:);
a=hist(aux,30);
% histfit(aux,30)
%%
xticks=floor(RBR.ST1.time_wave(1):1:RBR.ST1.time_wave(end));

f6=figure(6)
set(f6,'Units','centimeters','Position',[12.8852    0   15   13])
% 
% s1=subplot(2,1,1)
% yyaxis left
% p1=plot(AWAC_6734.wave.time,AWAC_6734.wave.Hm0/1.41,'color','#016FB9','linewidth',1.5);
% ylabel('H_{s} [m]')

ana(2,1)=737500.865973262;ana(2,2)=737502.278119836;
ana(3,1)=737506.125001883;ana(3,2)=737508.511042647;
ana(4,1)=737511.676198762;ana(4,2)=737513.307471529;
% p3=patch('Faces',[1 2 3 4],'Vertices',[ana(1,1) 0;ana(1,2) 0;ana(1,2) 2;ana(1,1) 2],'FaceColor',[139/256 0 0],'FaceAlpha',.3)

% yyaxis right
% p2=plot(AWAC_6734.wave.time,AWAC_6734.wave.Tm02,'color','#FF9505','linewidth',1.5);
% ax=gca;ax.YAxis(2).Color = '#FF9505';
% ylabel('T_{m02} [s]')
% xlim([RBR.ST1.time_wave(1)+2 RBR.ST1.time_wave(end)]);
% grid on
% title('Offshore non-directional wave parameters','fontsize',14)
% set(gca,'XTick',xticks,'ytick',[0:2:8],'fontsize',10,'xticklabel',[])

%%
% ratioF_Fr(ratioF_Fr<0)=nan;
% map=load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\functions\bathycolormap.mat');
map = colorcet('D1A', 'N', 10,'reverse',1);
colormap(flip(map))
% colormap(flipud(parula))
pcolor(RBR.ST1.time_wave,[1:9],abs(flip(ratioF_Fr')))
shading('flat')
xlim([RBR.ST1.time_wave(1)+2 RBR.ST1.time_wave(end)]);
set(gca,'XTick',xticks,'ytick',[0:1:9])
datetick('x','dd','keepticks','keeplimits')
title('$$|\varphi_{i,R} / (F_{i+1}/F_{i} - 1)|$$','Interpreter','latex','fontsize',14)
% title('\epsilon_i','fontsize',16)
ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels(2:2:end) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set
xlabel('Date in March and April 2019')
c=colorbar;
ylabel(c,'Ratio [ ]')
ylabel('Between station P_{(i+1)} and P_{(i)}')
box on
caxis([0 1])
set(gca,'fontsize',10)
patch('Faces',[1 2 3 4],'Vertices',[ana(2,1) 0;ana(2,2) 0;ana(2,2) 9;ana(2,1) 9],'FaceColor','none','EdgeColor','y','LineWidth',2)
patch('Faces',[1 2 3 4],'Vertices',[ana(3,1) 0;ana(3,2) 0;ana(3,2) 9;ana(3,1) 9],'FaceColor','none','EdgeColor','y','LineWidth',2)
patch('Faces',[1 2 3 4],'Vertices',[ana(4,1) 0;ana(4,2) 0;ana(4,2) 9;ana(4,1) 9],'FaceColor','none','EdgeColor','y','LineWidth',2)


%%
% s1.Position = [0.17    0.7795    0.65    0.1623];
% s2.Position = [0.17    0.1187    0.65    0.6];
% t=annotation('textbox','color','w','Position',[0.02 0.9 0.035 0.038], 'BackgroundColor',[100 147 167]./255,...
%     'EdgeColor',[100 147 167]./255,'String','A','FitBoxToText','off','fontsize',10,...
%     'horizontalalignment','center','verticalalignment','middle');%)%;
% 
% t=annotation('textbox','color','w','Position',[0.02 0.58 0.035 0.035], 'BackgroundColor',[100 147 167]./255,...
%     'EdgeColor',[100 147 167]./255,'String','B','FitBoxToText','off','fontsize',10,...
%     'horizontalalignment','center','verticalalignment','middle');%)%;
% box on
% 
% t=annotation('textbox','color','k','Position',[0.05 0.72 0.035 0.035], 'BackgroundColor','w',...
%     'EdgeColor','w','String','(Onshore)','FitBoxToText','off','fontsize',10,...
%     'horizontalalignment','center','verticalalignment','middle');%)%;
% 
% t=annotation('textbox','color','k','Position',[0.05 0.08 0.035 0.035], 'BackgroundColor','w',...
%     'EdgeColor','w','String','(Offshore)','FitBoxToText','off','fontsize',10,...
%     'horizontalalignment','center','verticalalignment','middle');%)%;
%%
print(gcf,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\Figure6_outputFalter.tiff','-dtiff','-r300')