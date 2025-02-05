load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd')
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_lev.mat')
load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\Sawhorse_instruments.mat')
%%
% b=0.014;u=RBR.ST9.ubr_LWT;
% v=9.35.*10^(-7);
% Re = b.*u./v;
% flag_Re=(Re>1500);
% rho=1026;
% flag_swell=RBR.ST9.ubr_LWT>0.1;
% flag=and(flag_Re,flag_swell);
for i=1:length(AQP11266.meanDepth_15min)
    Wavek =  wavek(AQP11266.freq_urms,AQP11266.meanDepth_15min(i)+1.09);  
    phij  =    (cosh (Wavek.*AQP11266.z(18)))./(cosh (Wavek.*AQP11266.z(32)));
    alphaj (:,i) = sqrt((phij.^2).*(sqrt(2*(AQP11266.Suu_ADP(i*257-256:i*257,32)+AQP11266.Svv_ADP(i*257-256:i*257,32))*(AQP11266.freq_urms(2)-AQP11266.freq_urms(1))))./...
        (sqrt(2*(AQP11266.Suu_ADP(i*257-256:i*257,18)+AQP11266.Svv_ADP(i*257-256:i*257,18))*(AQP11266.freq_urms(2)-AQP11266.freq_urms(1)))));

    alphar (i) = sqrt(sqrt(2*(nansum(phij.^2.*(AQP11266.Suu_ADP(i*257-256:i*257,32)+AQP11266.Svv_ADP(i*257-256:i*257,32))*(AQP11266.freq_urms(2)-AQP11266.freq_urms(1)))))./...
        (sqrt(2*(nansum(AQP11266.Suu_ADP(i*257-256:i*257,18)+AQP11266.Svv_ADP(i*257-256:i*257,18))*(AQP11266.freq_urms(2)-AQP11266.freq_urms(1))))));
end
% alphar(alphar>1 |alphar<0.4)=nan;
aux_frequencies = 1./[20:-.5:.5];
curve1=61.2.*aux_frequencies(1:34).^0.005573-59.88; 
curve2=61.2.*aux_frequencies(1:34).^0.005573-59.62; 

loglog(AQP11266.freq_urms,nanmean(alphaj,2))
hold on
loglog(aux_frequencies(1:34),curve1)
loglog(aux_frequencies(1:34),curve2)
loglog(RBR.ST2.f(:,1),61.2.*RBR.ST2.f(:,1).^0.005573-59.88,'b')
%%
flag_time1 = AQP11266.time_15min > 737500.865973262 & AQP11266.time_15min < 737502.278119836;
flag_time2 = AQP11266.time_15min > 737506.125001883 & AQP11266.time_15min < 737508.511042647;
flag_time3 = AQP11266.time_15min > 737511.676198762 & AQP11266.time_15min < 737513.307471529;
flag_time_aux=or(flag_time1,flag_time2);
flag_time_total=or(flag_time_aux,flag_time3);
Arms_infty = nanmean(AQP11266.ubr(:,18),2).*AQP11266.Wave.Tm01_ADP'/(2*pi);
arms_flag=Arms_infty(flag_time_total);alphar_flag=alphar(flag_time_total);
alphar_interp=interp1(AQP11266.time_15min,alphar,RBR.ST2.time_wave);
% for i=1:
alphaj_interp=interp1(AQP11266.time_15min,alphar,RBR.ST2.time_wave);
plot(AQP11266.time_15min,alphar)
hold on
plot(RBR.ST2.time_wave,alphar_interp)
% scatter(Arms_infty(flag_time_total),alphar(flag_time_total),'.')
[~,min_cutoff] = min(abs(RBR.ST10.f_out(:,1)-0.05));
[~,max_cutoff] = min(abs(RBR.ST10.f_out(:,1)-0.2));

% 0.8289.*violin_vectrino_aqd_Arms+0.04999
lambdaF_lv=3.14;
lambdaF_le_rbr = interp1(AQP11266.time_15min,lambdaF_le,RBR.ST2.time_wave);
 aux_frequencies = RBR.ST2.f (min_cutoff:max_cutoff,1);
% alphaj=61.2.*aux_frequencies.^0.005573-59.88; 

%%
% alphar=0.6;
rho=1026;g=9.8;
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

%         alphar = (2.452*(RBR.ST8.ubr_LWT.*RBR.ST8.Tm01/(2*pi))+0.09382)*2;

    RBR.(STi).s1 = (4.*RBR.(STi).delta)./(lambdaF_lv.*rho.*1.*ubr_avg.^3);
    RBR.(STi).s2 = (4.*RBR.(STi).delta)./(lambdaF_le_rbr.*rho.*1.*ubr_avg.^3);
        alphar=0.6;

        
%     RBR.(STi).s3 = (4.*RBR.(STi).delta)./(lambdaF_lv.*rho.*alphar_interp.^3.*ubr_avg.^3);
%     RBR.(STi).s4 = (4.*RBR.(STi).delta)./(lambdaF_le_rbr.*rho.*alphar_interp.^3.*ubr_avg.^3);
        
    RBR.(STi).s3 = (4.*RBR.(STi).delta)./...
        (lambdaF_lv.*rho.*alphar.*alphaj.^2.*ubr_avg.*ubj_avg.^2);
    RBR.(STi).s4 = (4.*RBR.(STi).delta)./...
        (lambdaF_le_rbr.*rho.*alphar.*alphaj.^2.*ubr_avg.*ubj_avg.^2);

%     RBR.(STi).fer_lowe_NoAtt =  RBR.(STi).Cd_lowe_NoAtt.*lambdaF_le_rbr;
%     RBR.(STi).fer_lowe_Att = RBR.(STi).Cd_lowe_Att.*lambdaF_le_rbr.*alphar.^3;
%     RBR.(STi).fer_lowe_Noflex = RBR.(STi).Cd_lowe_NoAtt.*lambdaF_lv.*alphar.^3;

%     RBR.(STi).Cd_lowe_NoAtt = (4.*RBR.(STi).delta)./(lambdaF_le_rbr.*rho.*1.*ubr_avg.^3);
%     RBR.(STi).Cd_lowe_Att = (4.*RBR.(STi).delta)./(lambdaF_le_rbr.*rho.*alphar.^3.*ubr_avg.^3);
%     RBR.(STi).Cd_lowe_Noflex = (4.*RBR.(STi).delta)./(lambdaF_lv.*rho.*alphar.^3.*ubr_avg.^3);
%     
%     RBR.(STi).fer_lowe_NoAtt =  RBR.(STi).Cd_lowe_NoAtt.*lambdaF_le_rbr;
%     RBR.(STi).fer_lowe_Att = RBR.(STi).Cd_lowe_Att.*lambdaF_le_rbr.*alphar.^3;
%     RBR.(STi).fer_lowe_Noflex = RBR.(STi).Cd_lowe_NoAtt.*lambdaF_lv.*alphar.^3;
end
%%
edges=[-0.05:0.07:1.95];
figure(8)
s1=subplot(2,2,1)
hold on
sce1_alpha1 = (RBR.ST7.s1(flag_time_total)+RBR.ST8.s1(flag_time_total))./2;
sce1_alpha1(sce1_alpha1 <-0.5 | sce1_alpha1> 0.4)=nan;
h1=histogram(sce1_alpha1,edges,'facecolor','#D0DDD7','edgecolor','#D0DDD7');
fitdist(sce1_alpha1','normal')
hold on
% title('Scenario 1, \mu = 0.29 ± 0.15')
title('Scenario 1, \mu = 0.11 ± 0.05')

box on
ylabel('Counts')
xlabel('C_D')
grid on
xlim([-0.55 2.3])
ylim([0 210])

s2=subplot(2,2,2)
hold on
sce2_alphainterp = (RBR.ST7.s2(flag_time_total)+RBR.ST8.s2(flag_time_total))'./2;
sce2_alphainterp(sce2_alphainterp>0.7)=nan;
% sce2_alphainterp(sce2_alphainterp>2.36 | sce2_alphainterp<-.06)=nan;
% histfit(sce2_alphainterp,25)
histogram(sce2_alphainterp,edges,'facecolor','#AAC0AA', 'edgecolor','#AAC0AA');
fitdist(sce2_alphainterp,'normal')
% legend('\mu = 1.51')
% xlim([-5 5])
% title ('Scenario 2, \mu = 0.98 ± 0.50')
title ('Scenario 2, \mu = 0.26 ± 0.1')

box on
ylabel('Counts')
xlabel('C_D')
grid on
xlim([-0.55 2.3])
ylim([0 90])

s3=subplot(2,2,3)
hold on
sce3_fullupright = ((RBR.ST7.s3(:,flag_time_total)+RBR.ST8.s3(:,flag_time_total))'./2); %PAY ATTENTION!
sce3_fullupright(sce3_fullupright>1.35)=nan;
histfit(sce3_fullupright,25)
histogram(sce3_fullupright,edges, 'facecolor','#216869', 'edgecolor','#216869');
fitdist(sce3_fullupright,'normal')
% title ('Scenario 3, \mu = 0.40 ± 0.19')
title ('Scenario 3, \mu = 0.51 ± 0.19')
box on
ylabel('Counts')
xlabel('C_D')
grid on
xlim([-0.55 2.3])
% ylim([0 70])

s4=subplot(2,2,4)
hold on
sce4_fullupright = (RBR.ST7.s4(:,flag_time_total)+RBR.ST8.s4(:,flag_time_total))'./2;
sce4_fullupright(sce4_fullupright<0.03)=nan;
% aux=find(sce3_fullupright
% histfit(sce3_fullupright,25)
histogram(sce4_fullupright,edges, 'facecolor','#023047', 'edgecolor','#023047');
fitdist(sce4_fullupright,'normal')
% title ('Scenario 3, \mu = 0.40 ± 0.19')
title ('Scenario 4, \mu = 1.37 ± 0.4')
box on
ylabel('Counts')
xlabel('C_D')
grid on
xlim([-0.55 2.3])
% ylim([0 70])

% [pdca_s1] = fitdist(sce1_alpha1','rayleigh')
% pdf_s1 = pdf(pdca_s1,edges);
% plot(edges, pdf_s1, 'color','#AAC0AA','linewidth',2)
% 
% % plot([-3*(edges(2)-edges(1)) edges(1:20)], [0 pdf_s1(1:20)], 'color','#AAC0AA','linewidth',2)
% hold on
% 
% [pdca_s2] = fitdist(sce2_alphainterp,'rayleigh')
% pdf_s2 = pdf(pdca_s2,edges)
% plot(edges, pdf_s2, 'color','#023047','linewidth',2)
% 
% [pdca_s3] = fitdist(sce3_fullupright,'rayleigh')
% pdf_s3 = pdf(pdca_s3,edges)
% plot(edges, pdf_s3, 'color','#216869','linewidth',2)
% 
% title ('PDF')
% xlim([-0.55 2.3])
% ylim([0 3])
% grid on
% legend('S 1','S 2', 'S 3')
% xlabel('C_D')

%%
% s1.Position = [0.1300    0.5979    0.3347    0.3269];
% s2.Position = [0.5703    0.5979    0.3347    0.3269];
% s3.Position = [0.3550    0.1257    0.3347    0.3269];
%%
print(gcf,'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Plots\chapter2\F16_Scenarios','-dtiff','-r900')

