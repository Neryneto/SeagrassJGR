load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd')
%% Create flags
% checking if the flow is turbulent. Re<1000 should be excluded
b=0.014;u=RBR.ST9.ubr_LWT;
v=9.35.*10^(-7);
Re = b.*u./v;
flag_Re=(Re>1000);

lambdaP = (pi*(d^2)/4)/(S+d)^2;
Ld1=2*hc_e*(1-lambdaP)/(Cd1*lambdaF);
KC=AQP11266.freestream.*AQP11266.Wave.Tm02_ADP'./Ld;
% for j=1:length(Re) 
% cumsum_st=cumsum(RBR.ST2.Syy(:,j).*100)./nansum(RBR.(ST).Syy(:,j));
% if cumsum_st(206)<80
% flag_swell(j)=1;
% else
% flag_swell(j)=0;
% end
% end

% flag_swell=RBR.ST9.Hrms>0.2;
flag_swell=RBR.ST9.ubr_LWT>0.1;

flag=and(flag_Re,flag_swell);
%%
f11=figure(11)
% s1=subplot(1,2,1);
% a=(RBR.ST7.Cd_lowe_Att+RBR.ST8.Cd_lowe_Att)/2;
a=(RBR.ST7.Cd_lowe_NoAtt+RBR.ST8.Cd_lowe_NoAtt)/2;

b=Re;
a(a<0)=nan;
plot (b,a,'.','color','#e9c46a');
% plot (b,a,'.r');
hold on
plot ([500:4500], 0.1+(925./[500:4500]).^3.16, '--','color','#2a9d8f','linewidth', 1.5) % bradley houser
plot ([500:4500], (925./[500:4500]).^3.52, '--','color','#e76f51','linewidth', 1.5) % bradley houser

% plot ([500:4500], 0.08+(2200./[500:4500]).^2.2, '--','color','#e76f51','linewidth', 1.5) %mendez 1999
% plot ([500:4500], 0.6152+(1019./[500:4500]).^1.243, '-','color','k','linewidth', 1.5) %mine
ylim([0.01 2])
% legend('$$\tilde{C}_{D}$$', 'Bradley and Houser (2009)', 'Mendez et al. (1999)', 'Interpreter','latex')
legend('$$\tilde{C}_{D}$$', 'Bradley and Houser (2009) 1', 'Bradley and Houser (2009) 2', 'Interpreter','latex')
xlabel('$$Re_{v}$$', 'Interpreter','latex','fontsize',16)
ylabel('$$\tilde{C}_{D}$$','Interpreter','latex','fontsize',16)
% s2=subplot(1,2,2);
% histfit((RBR.ST8.fer_lowe(flag)+RBR.ST7.fer_lowe(flag))/2,50,'normal') %1.093
% title('f_{er}')
% xlim([-1 2])
% ylabel('Counts')
% xlabel('f_{er}')
grid on
%%
%%
print(gcf,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\Figure18_Disc.tiff','-dtiff','-r300')
