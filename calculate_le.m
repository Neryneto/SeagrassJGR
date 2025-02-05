load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments')
load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\Sawhorse_instruments.mat')
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
B=(12*(rhow-rhos)*g*lv^3)/(E*t^2);
Arms_infty = nanmean(AQP11266.ubr(:,1),2).*AQP11266.Wave.Tm02_ADP'./(2*pi);
L=lv./Arms_infty;
le=lv.*(Ca.*L).^(-1/4);
% le(le>0.4)=nan;
%%
f1=figure
set(f1,'Units','pixels','OuterPosition',[40 302 1105 416])

s1=subplot(2,1,1)
% yyaxis left
ub=AQP11266.ubr(:,1);ub(ub<0.1)=nan;
hold on
% yyaxis right
Arms_infty(Arms_infty<.02 | Arms_infty>0.5)=nan;
a1=area(AQP11266.time_15min,movmean(Arms_infty,4,'omitnan'))
p1=plot(AQP11266.time_15min, movmean(le/lv,1,'omitnan'),'.-');
xlim([AQP11266.time_15min(1) AQP11266.time_15min(end)])
% datetick('x','dd','keepticks','keeplimits')
% xlabel ('Days in March 2019')
ylabel('[ ]');
set(gca,'xticklabel',[],'fontsize',12)
grid on
box on
ylim([0 1])

s2=subplot(2,1,2)
a2=area(AQP11266.time_15min, movmean(le*514*0.015,6,'omitnan'),'facecolor','#216869');
% area(AQP11266.time_15min, movmean(le,6,'omitnan'),'linewidth',2);
datetick('x','dd','keepticks','keeplimits')
xlabel ('Days in March 2019')
ylabel('\lambda_{f,e} [ ]', 'fontsize',12) 
set(gca,'ytick',0:.5:3,'fontsize',12)
grid on
box on
xlim([AQP11266.time_15min(1) AQP11266.time_15min(end)])
leg=legend([a1 p1 a2],{'A_\infty^{rms} [m/s]','l_{e,v}/l_v [ ]', '\lambda_{f,e} [ ]'},'location','northwest');

% legend ([pbot ptop pk],{'Near-bed','Top of canopy','D_{10}'})
%%
s1.Position =[0.07    0.56    0.78    0.35];
s2.Position =[0.07    0.1517    0.78    0.35];
leg.Position = [0.8557    0.6693    0.1061    0.2404];
%%
save ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_lev.mat','le')
print(gcf,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\Figure14_lev','-dtiff','-r900')
% l_ve = lv/(Ca*L)^1/3;