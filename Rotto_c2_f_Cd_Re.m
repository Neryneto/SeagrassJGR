load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd')
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_lev.mat')
load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\Sawhorse_instruments.mat')
%% Create flags
d=0.014;u=RBR.ST9.ubr_LWT;
v=9.35.*10^(-7);
Re = d.*u./v;
flag_Re=(Re>1000);
rho=1026;
flag_swell=RBR.ST9.ubr_LWT>0.1;
flag=and(flag_Re,flag_swell);

alphar=0.69;
% 0.8289.*violin_vectrino_aqd_Arms+0.04999
lambdaF_lv=3.2;
lambdaF_le_rbr = interp1(AQP11266.time_15min,lambdaF_le,RBR.ST2.time_wave);



Ld1=2.*le.*(1-lambdaP)./(Cd2.*lambdaF);
KC=AQP11266.freestream.*AQP11266.Wave.Tp_ADP'./d;
Kc_rbr = interp1(AQP11266.time_15min,KC,RBR.ST2.time_wave);

%%
% alphar=0.6;
for k=1:9 %for each instrument
    STi = sprintf('ST%d',k);
    STf = sprintf('ST%d',k+1);
    RBR.(STi).fej_lowe =[];RBR.(STi).fej_lowe1=[];RBR.(STi).fer_lowe=[];RBR.(STi).fer_lowe1=[];
    RBR.(STi).Cd_lowe_NoAtt = []; RBR.(STi).Cd_lowe_Att=[]; RBR.(STi).Cd_lowe_Noflex=[];
    RBR.(STi).slope = nanmean((RBR.(STf).meandepth-RBR.(STi).meandepth)./RBR.(STi).distance);

    h = (RBR.(STf).meandepth + RBR.(STi).meandepth)./2;
    
    ubj_avg = (RBR.(STf).ubj_LWT + RBR.(STi).ubj_LWT)./2;
    ubr_avg = (RBR.(STf).ubr_LWT + RBR.(STi).ubr_LWT)./2;
    a = (RBR.(STf).arms + RBR.(STi).arms)./2;

    RBR.(STi).deltaj = -(RBR.(STf).Fj - ((rbr_awac_kr(:,k)' +1).*RBR.(STi).Fj))/RBR.(STi).distance;
    RBR.(STi).delta = -(RBR.(STf).Fr - ((rbr_awac_kr(:,k)' +1).*RBR.(STi).Fr))/RBR.(STi).distance;

    RBR.(STi).Cd_lowe_NoAtt = (4.*RBR.(STi).delta)./(lambdaF_le_rbr.*rho.*1.*ubr_avg.^3);
    RBR.(STi).Cd_lowe_Att = (4.*RBR.(STi).delta)./(lambdaF_le_rbr.*rho.*alphar.^3.*ubr_avg.^3);
    RBR.(STi).Cd_lowe_Noflex = (4.*RBR.(STi).delta)./(lambdaF_lv.*rho.*alphar.^3.*ubr_avg.^3);
    
    RBR.(STi).fer_lowe_NoAtt =  RBR.(STi).Cd_lowe_NoAtt.*lambdaF_le_rbr;
    RBR.(STi).fer_lowe_Att = RBR.(STi).Cd_lowe_Att.*lambdaF_le_rbr.*alphar.^3;
    RBR.(STi).fer_lowe_Noflex = RBR.(STi).Cd_lowe_NoAtt.*lambdaF_lv.*alphar.^3;
end

%%
f11=figure(11)
% s1=subplot(1,2,1);
% a=(RBR.ST7.Cd_lowe_Att+RBR.ST8.Cd_lowe_Att)/2;
a=(RBR.ST7.Cd_lowe_NoAtt+RBR.ST8.Cd_lowe_NoAtt)/2;
a2=(RBR.ST7.Cd_lowe_Noflex+RBR.ST8.Cd_lowe_Noflex)/2;
Kc_rbr_flag = Kc_rbr(flag);a2_flag = a2(flag);
b=Re;
a(a<0)=nan;
a2(a2<0)=nan;

plot (b,a2,'.','color','#e9c46a');hold on
plot (b,a2,'.');
plot (Kc_rbr_flag,a2_flag,'.');


% plot (b(flag),a(flag),'.r');
plot ([500:4500], 0.1+(925./[500:4500]).^3.16, '--','color','#2a9d8f','linewidth', 1.5) % bradley houser
plot ([500:4500], (925./[500:4500]).^3.52, '--','color','#e76f51','linewidth', 1.5) % bradley houser

% plot ([500:4500], 0.08+(2200./[500:4500]).^2.2, '--','color','#e76f51','linewidth', 1.5) %mendez 1999
% plot ([500:4500], 0.6152+(1019./[500:4500]).^1.243, '-','color','k','linewidth', 1.5) %mine
ylim([0.01 2])
legend('$$\tilde{C}_{D}$$', 'Bradley and Houser (2009) 1', 'Bradley and Houser (2009) 2', 'Interpreter','latex')
xlabel('$$Re_{v}$$', 'Interpreter','latex','fontsize',16)
ylabel('$$\tilde{C}_{D}$$','Interpreter','latex','fontsize',16)
grid on
%%
%%
print(gcf,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\Figure18_Disc.tiff','-dtiff','-r300')
