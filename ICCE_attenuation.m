%comparing fej estimates for ICCE
load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd')
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m\AWAC_6734.mat')
%%

Hsig_rbr_interp_aux=[];
for i=1:10
    
    STi = sprintf('ST%d',i);
    
    if i<10
        distances10 (i) = RBR.(STi).distance;
    end
    distance_cumulative_aux (i,2) = nanmean(RBR.(STi).meandepth);    
end
distances_with_awac=[251 distances10];
distance_cumulative_with_awac(:,1) = [0 cumsum((distances_with_awac(1:end)))];
distance_cumulative_with_awac(1,2)=15.1;
distance_cumulative_with_awac(2:11,2)=(distance_cumulative_aux(:,2));

samplingRateIncrease = 3;

newXSamplePoints_awac = [linspace(min(distance_cumulative_with_awac(:,1)), max(distance_cumulative_with_awac(:,1)),...
    length(distance_cumulative_with_awac(:,1)) * samplingRateIncrease)];
smoothedY_awac = [spline(distance_cumulative_with_awac(:,1),...
    distance_cumulative_with_awac(:,2), newXSamplePoints_awac)];
%%
rho=1026;g=9.8;

for j = 1:length(AWAC_6734.wave.Hm0)
    [AWAC_6734.Lj(j),AWAC_6734.kj(j),AWAC_6734.sigmaj(j)]=disper(AWAC_6734.PRESSURE(j),AWAC_6734.wave.Tz(j));
end
      
awac_depth=interp1(AWAC_6734.TIME,AWAC_6734.PRESSURE,AWAC_6734.wave.time);
AWAC_6734.E = 0.5.*rho.*g.*AWAC_6734.wave.Hm0.^2;
par = (2.*AWAC_6734.kj'.*awac_depth)./(sinh(2.*AWAC_6734.kj'.*awac_depth));
AWAC_6734.Cg = 0.5.*(1+par).*(AWAC_6734.sigmaj'./AWAC_6734.kj'); clear par omega
AWAC_6734.F = AWAC_6734.Cg.*AWAC_6734.E;
awac_F=interp1(AWAC_6734.wave.time,AWAC_6734.F,RBR.ST1.time_wave);
awacCg=interp1(AWAC_6734.wave.time,AWAC_6734.Cg,RBR.ST1.time_wave);

% s1(:,1) = awac_F-(distance_cumulative_with_awac(2,1).*rho.*1.37.*RBR.ST10.ubr_LWT.^3)./4;
% s2(:,1) = awac_F-(distance_cumulative_with_awac(2,1).*rho.*1.37.*3.14.*RBR.ST10.ubr_LWT.^3)./4;
% s3(:,1) = awac_F-(distance_cumulative_with_awac(2,1).*rho.*1.37.*3.14.*0.6^3.*RBR.ST10.ubr_LWT.^3)./4;
% s4(:,1) = awac_F-(distance_cumulative_with_awac(2,1).*rho.*1.37.*1.57.*0.6^3.*RBR.ST10.ubr_LWT.^3)./4;

s1(:,1) = awac_F-(distance_cumulative_with_awac(2,1).*rho.*(0.37).*RBR.ST10.ubr_LWT.^3)./4;
s2(:,1) = awac_F-(distance_cumulative_with_awac(2,1).*rho.*(0.51.*3.14*1).*RBR.ST10.ubr_LWT.^3)./4;
s3(:,1) = awac_F-(distance_cumulative_with_awac(2,1).*rho.*(1.37.*3.14*0.6^3).*RBR.ST10.ubr_LWT.^3)./4;
s4(:,1) = awac_F-(distance_cumulative_with_awac(2,1).*rho.*(1.37.*1.1*0.6^3).*RBR.ST10.ubr_LWT.^3)./4;

% h1(:,1) = sqrt(s1(:,1)*2./(4*rho.*g));
% h2(:,1) = sqrt(s2(:,1)*2./(4.*rho.*g));
% h3(:,1) = sqrt(s3(:,1)*2./(4.*rho.*g));
% h4(:,1) = sqrt(s4(:,1)*2./(4.*rho.*g));

for k=2:10 %for each instrument
    STi = sprintf('ST%d',k);
%     s1(:,k) = s1(:,k-1)-(distance_cumulative_with_awac(k+1,1).*rho.*1.*RBR.(STi).ubr_LWT'.^3)./4;
%     s2(:,k) = s2(:,k-1)-(distance_cumulative_with_awac(k+1,1).*rho.*1.*3.14.*RBR.(STi).ubr_LWT'.^3)./4;
%     s3(:,k) = s3(:,k-1)-(distance_cumulative_with_awac(k+1,1).*rho.*1.*3.14.*0.6^3.*RBR.(STi).ubr_LWT'.^3)./4;
%     s4(:,k) = s4(:,k-1)-(distance_cumulative_with_awac(k+1,1).*rho.*1.*1.57.*0.6^3.*RBR.(STi).ubr_LWT'.^3)./4;
    
    s1(:,k) = s1(:,k-1)-(distance_cumulative_with_awac(k+1,1).*rho.*0.37.*RBR.(STi).ubr_LWT'.^3)./4;
    s2(:,k) = s2(:,k-1)-(distance_cumulative_with_awac(k+1,1).*rho.*(1.5.*3.14*1).*RBR.(STi).ubr_LWT'.^3)./4;
    s3(:,k) = s3(:,k-1)-(distance_cumulative_with_awac(k+1,1).*rho.*(2.*3.14*0.7^3).*RBR.(STi).ubr_LWT'.^3)./4;
    s4(:,k) = s4(:,k-1)-(distance_cumulative_with_awac(k+1,1).*rho.*(2*2.5*0.7^3).*RBR.(STi).ubr_LWT'.^3)./4;

    s1(s1<0)=nan;s2(s2<0)=nan;s3(s3<0)=nan;s4(s4<0)=nan;
    
%     h1(:,k) = sqrt(s1(:,k)'*2./(RBR.(STi).Cg.*rho.*g)); 
%     h2(:,k) = sqrt(s2(:,k)'*2./(RBR.(STi).Cg.*rho.*g)); 
%     h3(:,k) = sqrt(s3(:,k)'*2./(RBR.(STi).Cg.*rho.*g)); 
%     h4(:,k) = sqrt(s4(:,k)'*2./(RBR.(STi).Cg.*rho.*g)); 

end

h1 = sqrt(nanmean(s1).*2./(4.*rho.*g));
h2 = sqrt(nanmean(s2).*2./(4.*rho.*g));
h3 = sqrt(nanmean(s3).*2./(4.*rho.*g));
h4 = sqrt(nanmean(s4).*2./(4.*rho.*g));

h1_std=std(sqrt((s1).*2./(rho.*g)))
h2_std=std(sqrt((s2).*2./(rho.*g)))
h3_std=std(sqrt((s3).*2./(rho.*g)))
h4_std=std(sqrt((s4).*2./(rho.*g)))
%%
f5=figure(5)
set(f5,'Units','centimeters','Position',[4.2971   0   30    17])
yyaxis right
a=area(newXSamplePoints_awac,flip(-smoothedY_awac),-15.5,'Facecolor',[.5 .5 .5]);
a.EdgeColor=[.5 .5 .5];
a.FaceAlpha=0.3
hold on
% xlim([0 max(newXSamplePoints)])
xlim([0 1600])
% xlim([0 max(newXSamplePoints_awac)])
ylim([-15.5 -0])
ylabel('Depth [m]')

yyaxis left
% errorbar(flip(abs(1495.8-distance_cumulative_with_awac(1:end-1,1))),flip(nanmean(s1)),flip(nanstd(s1)./sqrt(1800)),'o--k','linewidth',1.4,...
%    'markeredgecolor','k','markerfacecolor','k')
% errorbar(flip(abs(1495.8-distance_cumulative_with_awac(1:end-1,1))),flip(nanmean(h1)),flip(nanstd(h1)./sqrt(1800)),'o--k','linewidth',1.4,...
%    'markeredgecolor','k','markerfacecolor','k')

hold on
% errorbar(flip(abs(1495.8-distance_cumulative_with_awac(1:end-1,1))),flip(nanmean(s2)),flip(nanstd(s2)./sqrt(1800)),'o--r','linewidth',1.4,...
%    'markeredgecolor','r','markerfacecolor','r')
errorbar(flip(abs(1495.8-distance_cumulative_with_awac(1:end-1,1))),flip((h2)),flip((h2_std)./sqrt(1800)),'o--r','linewidth',1.4,...
   'markeredgecolor','r','markerfacecolor','r')

% errorbar(flip(abs(1495.8-distance_cumulative_with_awac(1:end-1,1))),flip(nanmean(s3)),flip(nanstd(s3)./sqrt(1800)),'o--','color',[111 187 214]/256,'linewidth',1.4,...
%    'markeredgecolor',[111 187 214]/256,'markerfacecolor',[111 187 214]/256)
errorbar(flip(abs(1495.8-distance_cumulative_with_awac(1:end-1,1))),flip((h3)),flip((h3_std)./sqrt(1800)),'o--','color',[111 187 214]/256,'linewidth',1.4,...
   'markeredgecolor',[111 187 214]/256,'markerfacecolor',[111 187 214]/256)

% errorbar(flip(abs(1495.8-distance_cumulative_with_awac(1:end-1,1))),flip(nanmean(s4)),flip(nanstd(s4)./sqrt(1800)),'o--b','linewidth',1.4,...
%    'markeredgecolor','b','markerfacecolor','b')
errorbar(flip(abs(1495.8-distance_cumulative_with_awac(1:end-1,1))),flip((h4)),flip((h4_std)./sqrt(1800)),'o--b','linewidth',1.4,...
   'markeredgecolor','b','markerfacecolor','b')
ylim([0 1.1])
ax=gca;
ax.YAxis(2).Color = [0.5 0.5 0.5];
xlabel('Distance offshore from P10 [m]')
ylabel('H_{m0} [m]')
grid on
set(gca,'fontsize',14)
%%
print(f5,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\disipation_ICCE.tiff','-dtiff','-r300')