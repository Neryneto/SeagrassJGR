%% AFTER RYAN
map=load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\bathycolormap.mat');
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\model\direction_withvegetation\awaforcing_spc.mat')
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

caminho = 'C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Plots\chapter2\OS';
[dense_seagrass.x, dense_seagrass.y, dense_seagrass.z] = read_kml(fullfile(caminho,'dense_seagrass.kml'));
[dense_seagrass.x_utm, dense_seagrass.y_utm, dense_seagrass.zone] = deg2utm(dense_seagrass.y,dense_seagrass.x);
[intermediate_seagrass.x, intermediate_seagrass.y, intermediate_seagrass.z] = read_kml(fullfile(caminho,'intermediate_sg.kml'));
[intermediate_seagrass.x_utm, intermediate_seagrass.y_utm, intermediate_seagrass.zone] = deg2utm(intermediate_seagrass.y,intermediate_seagrass.x);
[intermediate_seagrass.x2, intermediate_seagrass.y2, intermediate_seagrass.z2] = read_kml(fullfile(caminho,'intermediate_sg2.kml'));
[intermediate_seagrass.x_utm2, intermediate_seagrass.y_utm2, intermediate_seagrass.zone] = deg2utm(intermediate_seagrass.y2,intermediate_seagrass.x2);
[sparse_seagrass.x, sparse_seagrass.y, sparse_seagrass.z] = read_kml(fullfile(caminho,'sparse_sg.kml'));
[sparse_seagrass.x_utm, sparse_seagrass.y_utm, sparse_seagrass.zone] = deg2utm(sparse_seagrass.y,sparse_seagrass.x);

%% Figure 8 - PLOT DIRECTION
f8=figure(8)
set(f8,'Units','centimeters','Position',[5    5    15    12])
% fprintf(f,'%s\t',datestr(RBR.ST3.time_wave(i),'yyyymmdd.HHMMSS'));
s=whos ('Hsig*');
dir=whos('Dir*');
C={};
Dir={};
Hs={};
u={};v={};
for i=2:length(s)
    C(1,i-1) = {s(i).name};
    C(2,i-1) = {eval(s(i).name)};
    
    Dir(1,i-1) = {dir(i).name};
    Hs(1,i-1) = {s(i).name};
    aux_dir = eval(dir(i).name);
    aux_hs = eval(s(i).name);
    aux_dir(isnan(aux_dir))=0;aux_hs(isnan(aux_hs))=0;
    for j=1:size(aux_dir,2)
        [u_aux(:,j),v_aux(:,j)]=compass2cart(aux_dir(:,j),aux_hs(:,j));
    end
    
    u(1,i) ={u_aux};
    v(1,i) ={v_aux};
    
end

umean={nanmean(cat(3,u{1,:}),3)};
vmean={nanmean(cat(3,v{1,:}),3)};
hs_mean={nanmean(cat(3,C{2,:}),3)};
[~,h]=contourf(Xp,Yp,hs_mean{1,1},'Linestyle','none')
colormap(flip(map.bathycolormap))
caxis([0.00 1.4])
c1=colorbar;
hold on

depth=Depth_20190308_140400;
depth(isnan(depth))=0;
[cc,hh]=contour(Xp,Yp,depth,'showtext','on','LineColor',[0 0 0])
set(hh,'ContourZLevel',10)
% contourf(Xp,Yp,Dir)
ylabel(c1,'H_s [m]', 'fontsize',10)
% plot(sta_rbr(:,1),sta_rbr(:,2),'.w','markersize',10)
plot(sta_rbr(:,1),sta_rbr(:,2),'.r','markersize',15)
plot(364498,6456752,'.k','markersize',15)


q=quiver(Xp(1:10:end,1:10:end),Yp(1:10:end,1:10:end),-umean{1,1}(1:10:end,1:10:end)./...
    sqrt(umean{1,1}(1:10:end,1:10:end).^2+vmean{1,1}(1:10:end,1:10:end).^2),...
    -vmean{1,1}(1:10:end,1:10:end)./sqrt(umean{1,1}(1:10:end,1:10:end).^2+vmean{1,1}(1:10:end,1:10:end).^2),0.5,'k');
q.Color='k';
% caxis([100 200])
box on
set(gca,'xtick',363300:350:365400,'xticklabel',[363300:350:365400],'ytick',6456000:400:6458600,'yticklabel',...
    ({'6456000',     '6456400',     '6456800',     '6457200',     '6457600',     '6458000', '6458400' }),'fontsize',11)

for i=[1 5 8 10]
text(sta_rbr(i,1)+70,sta_rbr(i,2),['P' num2str(11-i)],'BackgroundColor', 'w', 'fontsize',10)
end
text(364498+70,6456752,'CW1', 'fontsize',10)
grid on
xlabel('Easting [m]')
ylabel('Northing [m]')
set(gca, 'fontsize',10)

%%
print(gcf,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\Figure4_ModeledWaveDirection.tiff','-dtiff','-r300')