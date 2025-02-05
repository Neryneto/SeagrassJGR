load ('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd')
%% Create flags
% checking if the flow is turbulent. Re<1000 should be excluded
b=0.014;u=RBR.ST9.ubr_LWT;
v=9.35.*10^(-7);
Re = b.*u./v;
flag_Re=(Re>1500);

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
% xticks=RBR.ST9.time_wave(1):1:RBR.ST9.time_wave(end);
% 
f7=figure (7)
set(f7,'Units','centimeters','Position',[5    5    12    10])
fer8=RBR.ST8.fer_lowe;
fer7=RBR.ST7.fer_lowe;
% % h2=histfit(fer8,50,'normal') %1.093
% hold on
% h2=histfit(fer8(flag),50,'normal') %1.093
% title('ST8')
% % xlim([-2 2])
% fitdist(fer8(flag)','normal') %1.093
% 
% s4=subplot(2,2,4)
hold on
fer78=((fer7(flag)'+fer8(flag)')./2);
fer78(fer78<-.3 | fer78>1.5)=nan;
h2=histogram(fer78,20,'Normalization','probability'); %1.093

% h2=histfit(fer78,20,'normal') %1.093
pd=fitdist(fer78,'normal') %1.093

x = -1:0.086:2;
y = pdf(pd, x);
y = y / sum( pdf(pd, x) );
plot(x, y, 'r','linewidth',1.5) % Plot distribution as a black dashed line

xlim([-.5 1.6])
ylabel('[%]')
set(gca,'fontsize',10)
% title('$$f_{e}^R$$ between P7 - P9','Interpreter','latex','fontsize',14)

xlabel('$$f_{e,R}$$ [ ]','Interpreter','latex','fontsize',12);
box on
%%
% f11=figure(11)
% % s1=subplot(1,2,1);
% a=(RBR.ST8.fer_lowe+RBR.ST7.fer_lowe)/2;
% b=Re;
% plot (Re,(RBR.ST8.fer_lowe+RBR.ST7.fer_lowe)/2,'.','color','#e9c46a');
% hold on
% plot ([500:4500], 0.1+(925./[500:4500]).^3.16, '--','color','#2a9d8f','linewidth', 1) % bradley houser
% plot ([500:4500], 0.08+(2200./[500:4500]).^2.2, '--','color','#e76f51','linewidth', 1) %mendez 1999
% plot ([500:4500], 0.4651+(997.2./[500:4500]).^2.077, '-','color','k','linewidth', 1.5) %mendez 1999
% ylim([0.01 5])
% legend('f_{er}', 'Bradley and Houser (2009)', 'Mendez et al. (1999)')
% ylabel('f_{er}')
% xlabel('Re_v')

% s2=subplot(1,2,2);
% histfit((RBR.ST8.fer_lowe(flag)+RBR.ST7.fer_lowe(flag))/2,50,'normal') %1.093
% title('f_{er}')
% xlim([-1 2])
% ylabel('Counts')
% xlabel('f_{er}')
%%
%%
print(gcf,'C:\Users\NeryNeto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\article\chapter2\Figures\Figure7_fer.tiff','-dtiff','-r300')
