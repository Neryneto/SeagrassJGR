BrickerMonismith = load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\BrickerMonismith.mat');
BrickerMonismith_vectrino = load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\BrickerMonismith_vectrino.mat');
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Quadpod_instruments.mat');
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments')
load ('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd')
%% Calculate spectra

Vec_h=[];Svel = [];

heights = 0.05:0.01:max(Vec_h)+0.01;heights=heights';

for num=1:122 %each profile
    Vec_h(num,1) = nanmean(Vectrino(num).BottomDistance);
end
for i = 1:floor(length(ADV1670.PRES)/8192)
    ADV1670.press_15min (i) = nanmean(ADV1670.PRES(i*8192-8191:i*8192));
end
%% Calculate attenuation from Vector to Vectrino for each frequency
% bottom = [35 47 61 75];
% top = [46 58 72 86];
bottom = [35 36 37 47 48 49 61 62 63 75 76 77];
top = [46 46 46 58 58 58 72 72 72 86 86 86];

avg=struct;
for j=1:12;
    aux_bot=bottom(j);
    aux_top=top(j);
    time_vectrino = Vectrino(aux_bot).InitialTime_Matlab;
    time_vectrino = Vectrino(aux_bot).InitialTime_Matlab;
%     [l c] = min(abs(time_vectrino-AQD11898.TIME_RAW));
%     if j==1 | j==4 |j==7 |j==10
        time_vectrino_end = Vectrino(aux_top).Time(end);
    [ll c] = min(abs(time_vectrino_end-ADV1670.time_ft));
%     [ll cc] = min(abs(time_vectrino_end-AQD11898.TIME_RAW));
    
%     for i=1:floor((cc-c)/(60*1))
% %     avg_minute_vel(j,i) = nanmean(sqrt(ADV1670.vel_maj(c+481*(i-1):c+481*i).^2+ADV1670.vel_min(c+481*(i-1):c+481*i).^2));
%     avg_minute_vel(j,i) = nanmean(sqrt(AQD11898.vel_maj(c+61*(i-1):c+61*i).^2+AQD11898.vel_min(c+61*(i-1):c+61*i).^2));
%     end
%     aux=avg_minute_vel(j,:);
%   minmax(j)=nanmin(aux(aux>0))./nanmax(aux);  
%     end
    aux_frequencies = 1./[20:-.5:.5];
    Wavek =  wavek(aux_frequencies,ADV1670.press_15min(c)+.55);
    
    p=BrickerMonismith_vectrino.S_uwave_uwave{1,aux_bot};p(p<0)=nan;BrickerMonismith_vectrino.S_uwave_uwave{1,aux_bot}=p;
    t=BrickerMonismith_vectrino.S_uwave_uwave{1,aux_top};t(t<0)=nan;BrickerMonismith_vectrino.S_uwave_uwave{1,aux_top}=t;
    
    for i=1:length(aux_frequencies)-1
        avg_spec_bot (:,i) = nanmean(BrickerMonismith_vectrino.S_uwave_uwave{1,aux_bot}(BrickerMonismith_vectrino.fm{aux_bot}>aux_frequencies(i) &...
            BrickerMonismith_vectrino.fm{aux_bot}<=aux_frequencies(i+1)));
        avg_spec_top (:,i) = nanmean(BrickerMonismith_vectrino.S_uwave_uwave{1,aux_top}(BrickerMonismith_vectrino.fm{aux_top}>aux_frequencies(i) &...
            BrickerMonismith_vectrino.fm{aux_top}<=aux_frequencies(i+1)));
    end
    
    phij (:,j) =    (cosh (Wavek.*Vectrino(aux_top).BottomDistance))./(cosh (Wavek.*Vectrino(aux_bot).BottomDistance));
    alphaj (:,j) = sqrt((phij(1:end-1,j).^2).*(avg_spec_bot')./avg_spec_top');
end
% mask = nanstd(alphaj,0,2)>0.35;alphaj(mask,:)=nan;
mask = nanstd(alphaj,2)>0.35;alphaj(mask,:)=nan;
curve1=61.2.*aux_frequencies(1:34).^0.005573-59.88; 
curve2=61.2.*aux_frequencies(1:34).^0.005573-59.62; 

%%
% checking if the flow is turbulent. Re<1000 should be excluded
b=0.014;u=RBR.ST9.ubr_LWT;
v=9.35.*10^(-7);
Re = b.*u./v;
flag_Re=(Re>1700);

flag_swell=RBR.ST9.ubr_LWT>0.1;

flag=and(flag_Re,flag_swell);

[~,min_cutoff] = min(abs(RBR.ST10.f_out(:,1)-0.05));
[~,max_cutoff] = min(abs(RBR.ST10.f_out(:,1)-0.2));

%%
f15=figure(15)
set(f15,'Outerposition',[495.4000   51.4000  747.2000  628.0000])
% aa=nanmean(RBR.ST8.fej_lowe(:,flag_swell),2); bb=nanmean(RBR.ST7.fej_lowe(:,flag_swell),2);
fe78_mean=(RBR.ST8.fej_lowe(:,flag_swell)+RBR.ST7.fej_lowe(:,flag_swell))./2;
fej78_error=nanstd(fe78_mean,2,1)./sqrt(size(fe78_mean,1)/sqrt(155)); %155 represents effective degrees of freedom

e=errorbar(RBR.ST9.f(min_cutoff:max_cutoff,1),nanmean(fe78_mean,2),fej78_error,'vertical','-sk','MarkerSize',8,...
    'MarkerEdgeColor','#e9c46a','MarkerFaceColor','#e9c46a','capsize',0);
e.Color = '#e9c46a';

hold on

a=aux_frequencies(1:end-1);b=nanmean(alphaj,2);

s2r=plot (aux_frequencies(1:34),repmat(0.11.*3.14.*1,length(aux_frequencies(1:34))),'color','#016FB9','linewidth',1.2); 
s2f=plot (aux_frequencies(1:34),repmat(0.29.*1.3.*1,length(aux_frequencies(1:34))),'color','#7FB069','linewidth',1.2) 

s3r=plot (aux_frequencies(1:34),0.59.*3.14.*1.*...
    (61.2.*aux_frequencies(1:34).^0.005573-59.75).^2,'color','k','linewidth',1.2) 

s3f=plot (aux_frequencies(1:34),1.4.*1.3.*0.6.*...
    (61.2.*aux_frequencies(1:34).^0.005573-59.75).^2,'color','r','linewidth',1.2) 

xlim([0.05 0.2])
ylim([-0.05 1.6])

hold on
xlabel('Frequency [Hz]')
ylabel('f_{e,j} [ ]')
set(gca,'fontsize',14)
grid on

legend([e s2r(1) s2f(1) s3r(1) s3f(1)],{['$$f_{e,j} \ P7-P8$$ and $$P8-P9$$'],...
    ['2-R ($$\alpha_{w,R}=1,\ \tilde{C}_{D}, \ \lambda_{f}$$)'],...
    ['2-F ($$\alpha_{w,R}=1,\ \tilde{C}_{D}, \ \lambda_{f,e}$$)'],...
    ['3-R ($$\alpha_{w,j},\ C_{D}, \ \lambda_{f}$$)'],...
    ['3-F ($$\alpha_{w,j},\ C_{D}, \ \lambda_{f,e}$$)']},...
    'Interpreter','latex','location','Northwest','Fontsize',12)

%%
print(gcf,'C:\Users\NeryNeto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\article\chapter2\Figures\Figure14_fej_compared.tiff','-dtiff','-r300')

