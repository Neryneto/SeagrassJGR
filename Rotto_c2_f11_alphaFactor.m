BrickerMonismith = load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\BrickerMonismith.mat');
BrickerMonismith_vectrino = load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\functions\BrickerMonismith_vectrino.mat');
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Quadpod\Quadpod_instruments.mat');
load('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments');
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

bottom = [35 36 37 47 48 49 61 62 63 75 76 77];
top = [46 46 46 58 58 58 72 72 72 86 86 86];
 aux_frequencies = 1./[20:-.5:.5];
avg=struct;
for j=1:length(bottom);
    aux_bot=bottom(j);
    aux_top=top(j);
    plot(j,Vectrino(aux_bot).BottomDistance,'.r')
    hold on
plot(j,Vectrino(aux_top).BottomDistance,'.b')
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
mask = nanstd(alphaj,2)>0.35;alphaj(mask,:)=nan;
f14=figure(14)
a=aux_frequencies(1:35);b=nanmean(alphaj(1:35,:),2);
e=errorbar(aux_frequencies(1:35),nanmean(alphaj(1:35,:),2),nanstd(alphaj(1:35,:),2)/sqrt(39),'vertical','-sk','MarkerSize',8,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
% b(b>1)=nan;
a=aux_frequencies(1:35);b=nanmean(alphaj(1:35,:),2);
hold on
p1=plot ([aux_frequencies(1:33)], ...
    [61.2.*aux_frequencies(1:33).^0.005573-59.75],'-r','linewidth',2) 
p123=plot([0.25 0.33],[1 1],'-r','linewidth',2)
% p1=plot (aux_frequencies(1:31),8.608.*aux_frequencies(1:31).^0.04522-7.082,'-r','linewidth',2) 

curve1=61.2.*aux_frequencies(1:33).^0.005573-59.88; 
curve2=61.2.*aux_frequencies(1:33).^0.005573-59.62; 

x2 = [aux_frequencies(1:33), fliplr(aux_frequencies(1:33))];
inBetween = [curve1, fliplr(curve2)];
p12=fill(x2, inBetween,'r');
p12.FaceAlpha=0.15;
p12.EdgeColor='none'

x3 = [aux_frequencies(33:35), fliplr(aux_frequencies(33:35))];
inBetween = [repmat(0.87,1,3), fliplr(repmat(1.13,1,3))];
p12=fill(x3, inBetween,'r');
p12.FaceAlpha=0.15;
p12.EdgeColor='none'
% p_min=plot (aux_frequencies(1:31),-148.1+149.2.*aux_frequencies(1:31).^0.001546,'-r','linewidth',2) 
freqs=aux_frequencies(1:end-1);
t1=sqrt(nanmin(alphaj,[],2));
t2=(nanmin(alphaj,[],2));

ylim([0 1.5])
xlim([0 0.35])
e.Color = 'k'; hold on
xlabel('Frequency [Hz]')
ylabel('$\alpha_{w,j}$','interpreter','latex')
set(gca,'fontsize',14)
grid on

legend([e p1],{'Measured $\alpha_{w,j}$','$R^2=0.82, p=0.001, \nu=16$'},'interpreter','latex','location','SouthEast','Fontsize',12)
set(gca, 'YTick',0:.25:1.5,'fontsize',10)
set(gca,'fontsize',12)
text (0.12,0.55,'$\alpha_{w,j}=61.2*f^{0.005573}-59.75$','interpreter','latex','Fontsize',14)
%%
% alphar=sqrt(phij.^2.*alphaj)
%%
print(f14,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\Figure11_alpha_w.tiff','-dtiff','-r300')

