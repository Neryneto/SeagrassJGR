%% AFTER RYAN
bat_old=dlmread('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awac_forcing_spec_noFric_statio\bat_bickley_coarse.dat');
map=load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\functions\bathycolormap.mat');
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\direction_withvegetation\awaforcing_spc.mat')
%%

azymuth = 180-acosd((-6456977+6458205)./sqrt((363817-364383)^2+(6458205-6456977)^2));

sta_rbr=[363817	6458205; %RBR1  #
363864	6458112; %RBR2  #
363908	6458022; %RBR3  #
363951	6457929; %RBR4  #
363991	6457838; %RBR5  #78097
364067	6457660; %RBR6  #77841
364119	6457570; %RBR7  #78096
364147	6457482; %RBR8  #77957
364273	6457221; %RBR9  #124872
364383	6456977]; %RBR10 #77765];

caminho = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\OS';
[dense_seagrass.x, dense_seagrass.y, dense_seagrass.z] = read_kml(fullfile(caminho,'dense_seagrass.kml'));
[dense_seagrass.x_utm, dense_seagrass.y_utm, dense_seagrass.zone] = deg2utm(dense_seagrass.y,dense_seagrass.x);
[intermediate_seagrass.x, intermediate_seagrass.y, intermediate_seagrass.z] = read_kml(fullfile(caminho,'intermediate_sg.kml'));
[intermediate_seagrass.x_utm, intermediate_seagrass.y_utm, intermediate_seagrass.zone] = deg2utm(intermediate_seagrass.y,intermediate_seagrass.x);
[intermediate_seagrass.x2, intermediate_seagrass.y2, intermediate_seagrass.z2] = read_kml(fullfile(caminho,'intermediate_sg2.kml'));
[intermediate_seagrass.x_utm2, intermediate_seagrass.y_utm2, intermediate_seagrass.zone] = deg2utm(intermediate_seagrass.y2,intermediate_seagrass.x2);
[sparse_seagrass.x, sparse_seagrass.y, sparse_seagrass.z] = read_kml(fullfile(caminho,'sparse_sg.kml'));
[sparse_seagrass.x_utm, sparse_seagrass.y_utm, sparse_seagrass.zone] = deg2utm(sparse_seagrass.y,sparse_seagrass.x);

%%

% PLOT DEPTH
% contourf(Xp,Yp,Depth,[0:2:20])
% colormap(map.bathycolormap)
% c1=colorbar;
% colorTitleHandle = get(c1,'Title');titleString = 'Depth (m)';
% set(colorTitleHandle ,'String',titleString);
% hold on
% plot(sta_rbr(:,1),sta_rbr(:,2),'.k','markersize',10)
% plot(sta_rbr(2:5,1),sta_rbr(2:5,2),'.w','markersize',10)
% box off
% title('Depth used for the grid')
% print(gcf,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\appendix1_depth.tiff','-dtiff','-r300')
%% Figure 8 - PLOT DIRECTION
f8=figure(8)
[~,h]=contourf(Xp,Yp,Hsig_20190322_100400,'Linestyle','none')
colormap(flip(map.bathycolormap))
caxis([0.00 1.4])
c1=colorbar;
hold on

depth=Depth_20190308_140400;
depth(isnan(depth))=0;
[cc,hh]=contour(Xp,Yp,depth,'showtext','off','LineColor',[0 0 0])
set(hh,'ContourZLevel',10)
% contourf(Xp,Yp,Dir)
ylabel(c1,'H_s [m]', 'fontsize',12)
% plot(sta_rbr(:,1),sta_rbr(:,2),'.w','markersize',10)
plot(sta_rbr(:,1),sta_rbr(:,2),'.r','markersize',15)
plot(364498,6456752,'.k','markersize',15)

% uxtheta=wrapTo360(Dir_20190322_110400(1:10:end,1:10:end));uxmag=Propag_20190322_110400(1:10:end,1:10:end);
uxtheta=wrapTo360(Dir_20190401_110400(1:10:end,1:10:end));uxmag=Hsig_20190401_110400(1:10:end,1:10:end);
uxtheta(isnan(uxtheta))=0;uxmag(isnan(uxmag))=0;
for i=1:size(uxtheta,2)
[u(:,i),v(:,i)]=compass2cart(uxtheta(:,i),uxmag(:,i));
end
q=quiver(Xp(1:10:end,1:10:end),Yp(1:10:end,1:10:end),-u./sqrt(u.^2+v.^2),-v./sqrt(u.^2+v.^2),0.5,'k');
q.Color='k';
% caxis([100 200])
box on
title('Wave direction on 22/03/2019 11:04:00','fontsize',14)
set(gca,'xtick',363300:350:365400,'xticklabel',[363300:350:365400],'ytick',6456000:400:6458600,'yticklabel',...
    ({'6456000',     '6456400',     '6456800',     '6457200',     '6457600',     '6458000', '6458400' }),'fontsize',11)

for i=[1 5 8 10]
text(sta_rbr(i,1)+70,sta_rbr(i,2),['P' num2str(11-i)],'BackgroundColor', 'w', 'fontsize',10)
end
text(364498+70,6456752,'CW1', 'fontsize',10)

%%
print(gcf,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\Figure8_ModeledWaveDirection.tiff','-dtiff','-r300')