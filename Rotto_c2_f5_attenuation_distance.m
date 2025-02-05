load 'C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd.mat'
load 'C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m\AWAC_cur_wav.mat'; 
AWAC=str; %1912 06 Mar 12:04 to 15th April 7:34
%% Load and prepare modeled RBR; synchronize RBR and AWAC
for i=1:10
    ST = sprintf('ST%d',i);
    results = load (['C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\model\awac_nofric\rbr',...
        num2str(i) '_spc.dat']);
    RBR_Model.nofricSpec.(ST).time_wave = datenum(num2str(results(:,1),'%f'),'yyyymmdd.HHMMSS'); %1779 rows  06 Mar 12:34 to 12th April 13:34
    RBR_Model.nofricSpec.(ST).Hsig = results(:,2);
    RBR_Model.nofricSpec.(ST).Tp = results(:,3);
    RBR_Model.nofricSpec.(ST).MeanDir = results(:,4);
    RBR_Model.nofricSpec.(ST).PeakDir = results(:,5);
end

[l c]=min(abs(RBR_Model.nofricSpec.ST1.time_wave(end)-AWAC.wav.time));
fn = fieldnames(AWAC.wav);
for k=1:numel(fn)
    if( isnumeric(AWAC.wav.(fn{k})) )
       AWAC.wav.(fn{k}) = AWAC.wav.(fn{k})(2:c) ;
    end
end

% fn = fieldnames(RBR_Model.nofricSpec.ST10);
% for i=1:10
%     ST = sprintf('ST%d',i);   
%     for k=2:numel(fn)
%         if ( isnumeric(RBR_Model.nofricSpec.(ST).(fn{k})) )
%             RBR_Model.nofricSpec.(ST).interp.(fn{k}) = interp1(RBR_Model.nofricSpec.ST2.time_wave, RBR_Model.nofricSpec.(ST).(fn{k}), RBR.(ST).time_wave );
%         end
%     end
% end

for k=1:10
   ST = sprintf('ST%d',k);
   ST2 = sprintf('ST%d',k+1);
   h_modeled_spec_aux (:,k) = RBR_Model.nofricSpec.(ST).Hsig./AWAC.wav.Hm0;  
   h_modeled_meanwaveheight (:,k) = nanmean(RBR_Model.nofricSpec.(ST).Hsig);
   h_modeled_stdwaveheight (:,k) = nanstd(RBR_Model.nofricSpec.(ST).Hsig);
   h_modeled_spec_notNorm_aux (:,k) = RBR_Model.nofricSpec.(ST).Hsig;
end
h_modeled_spec_notNorm=fliplr(h_modeled_spec_notNorm_aux);
h_modeled_spec=fliplr(h_modeled_spec_aux);
%% Prepare RBR measured data
Hsig_rbr_interp_aux=[];
for i=1:10
    
    STi = sprintf('ST%d',i);
    
    Hsig_rbr_interp_aux (:,i) = interp1(RBR.(STi).time_wave, RBR.(STi).Hsig_swell_out, RBR_Model.nofricSpec.ST2.time_wave)	;
    if i<10
        distances10 (i) = RBR.(STi).distance;
    end
    distance_cumulative_aux (i,2) = nanmean(RBR.(STi).meandepth);    
end
Hsig_rbr_interp=fliplr(Hsig_rbr_interp_aux);
Hsig_measured=[Hsig_rbr_interp AWAC.wav.Hm0];
% distances10(10)=0;
distances_with_awac=[251 distances10];
distance_cumulative_with_awac(:,1) = [0 cumsum(fliplr(distances_with_awac(1:end)))];
distance_cumulative_with_awac(end,2)=15.1;
distance_cumulative_with_awac(1:10,2)=flip(distance_cumulative_aux(:,2));

samplingRateIncrease = 3;
% newXSamplePoints = linspace(min(distance_cumulative(:,1)), max(distance_cumulative(:,1)), length(distance_cumulative(:,1)) * samplingRateIncrease);
% smoothedY = spline(distance_cumulative(:,1), distance_cumulative(:,2), newXSamplePoints);

newXSamplePoints_awac = [linspace(min(distance_cumulative_with_awac(:,1)), max(distance_cumulative_with_awac(:,1)),...
    length(distance_cumulative_with_awac(:,1)) * samplingRateIncrease)];
smoothedY_awac = [spline(distance_cumulative_with_awac(:,1),...
    distance_cumulative_with_awac(:,2), newXSamplePoints_awac)];
%%
free_stream_col = '#7FB069'; %green
top_canopy_col = 'k';
fc='#016FB9';
%%
f5=figure(5);
set(f5,'Units','centimeters','Position',[12.8852    0   18   13]);

% s1=subplot(3,1,1);

yyaxis right
a=area(newXSamplePoints_awac,(-smoothedY_awac),-15.5,'Facecolor',[.5 .5 .5]);
a.EdgeColor=[.5 .5 .5];
a.FaceAlpha=0.3;
hold on
% xlim([0 max(newXSamplePoints)])
xlim([0 1650]);
% xlim([0 max(newXSamplePoints_awac)])
ylim([-15.3 -0]);
ylabel('Depth [m]');

yyaxis left
%not normalized
errorbar(distance_cumulative_with_awac(:,1),nanmean(Hsig_measured),nanstd(Hsig_measured),'o--','linewidth',1.4,...
   'markeredgecolor',fc,'markerfacecolor',fc)
hold on
e=errorbar(distance_cumulative_with_awac(:,1)',flip([nanmean(Hsig_measured(:,end)) h_modeled_meanwaveheight]),...
    flip([nanstd(Hsig_measured(:,end)) h_modeled_stdwaveheight]./2),'o--','linewidth',1.4,...
   'markeredgecolor','k','markerfacecolor','k')
e.Color='k';
ax=gca;
ax.YAxis(2).Color = [0.5 0.5 0.5];
ax.YAxis(1).Color = fc;

set(gca,'fontsize',10,'xtick',[0:200:1600],'XTick',0:200:1600,'XTickLabel',[])
ylabel('H_{s,i} [m]')
grid on
ylim([0 1.5])
pdense=plot([1000 1000],[0 1.5],'-','linewidth',2,'color','#66101F');
text(1480,1.2,'CW1','fontsize',10)
for i=1:10
%     if i==1
%         text(distance_cumulative_with_awac(i)-20,.5,['P' num2str(11-i)],'fontsize',12)
%     else
    text(distance_cumulative_with_awac(i)-15,1.2,['P' num2str(11-i)],'fontsize',10)
    end
% end
% xlabel('Distance offshore from P10 [m]')
%%
s2=subplot(3,1,2)
for i=1:9
    e1=errorbar((distance_cumulative_with_awac(i+1,1)+distance_cumulative_with_awac(i,1))/2,nanmean(Hsig_measured(:,i)./Hsig_measured(:,i+1)),...
        nanstd(Hsig_measured(:,i)./Hsig_measured(:,i+1)),'o','linewidth',1.4,...
        'markeredgecolor',free_stream_col,'markerfacecolor',free_stream_col);
    e1.Color=free_stream_col;
    hold on
    
    e2=errorbar((distance_cumulative_with_awac(i+1,1)+distance_cumulative_with_awac(i,1))/2,nanmean(Hsig_measured(:,i)./Hsig_measured(:,i+1))./distances10(i),...
        nanstd(Hsig_measured(:,i)./Hsig_measured(:,i+1)),'o','linewidth',1.4,...
        'markeredgecolor','r','markerfacecolor','r');
    e2.Color='r';
    
    if i<10
        
        e2=errorbar((distance_cumulative_with_awac(i+1,1)+distance_cumulative_with_awac(i,1))/2,nanmean(h_modeled_spec_notNorm(2:end,i)./h_modeled_spec_notNorm(2:end,i+1)),...
            nanstd(h_modeled_spec_notNorm(2:end,i)./h_modeled_spec_notNorm(2:end,i+1)),'o',...
            'linewidth',1.4,'markeredgecolor',top_canopy_col,'markerfacecolor',top_canopy_col);
        e2.Color = top_canopy_col
    else
        
        e2=errorbar((distance_cumulative_with_awac(i+1,1)+distance_cumulative_with_awac(i,1))/2,nanmean(h_modeled_spec_notNorm(2:end,i)./Hsig_measured(2:end,i+1)),...
            nanstd(h_modeled_spec_notNorm(2:end,i)./Hsig_measured(2:end,i+1)),'o',...
            'linewidth',1.4,'markeredgecolor',top_canopy_col,'markerfacecolor',top_canopy_col);
        e2.Color = top_canopy_col;
    end
   if i==6 | i==4 | i==2
    text((distance_cumulative_with_awac(i+1,1)+distance_cumulative_with_awac(i,1))/2-55,...
        1.18,['P' num2str(10-i) '-P' num2str(11-i)],'fontsize',10);
   else
       text((distance_cumulative_with_awac(i+1,1)+distance_cumulative_with_awac(i,1))/2-55,...
        1.27,['P' num2str(10-i) '-P' num2str(11-i)],'fontsize',10);
   end
end
set(gca,'fontsize',10,'xtick',[0:200:1600],'xticklabel',[])
grid on
xlim([0 1650])
ylabel('H_{s,i+1}/H_{s,i} [ ]')
grid on
ylim([0.5 1.3])
pdense=plot([1000 1000],[0 1.5],'-','linewidth',2,'color','#66101F');
l=legend ([e1 e2], {['Measured time-averaged'],['Modeled without dissipation']},...
    'location','southwest','orientation','vertical','fontsize',10)
%%
s3=subplot(3,1,3)
for i=1:9
%Using equations
e3=errorbar((distance_cumulative_with_awac(i+1,1)+distance_cumulative_with_awac(i,1))/2,nanmean((h_modeled_spec_notNorm(2:end,i)./h_modeled_spec_notNorm(2:end,i+1))-1),...
            nanstd((h_modeled_spec_notNorm(2:end,i)./h_modeled_spec_notNorm(2:end,i+1))-1),'o','linewidth',1.4,...
   'markeredgecolor','#FEB603','markerfacecolor','#FEB603');
e3.Color='#FF9505';
hold on
end

xlabel('Distance offshore from P10 [m]')
ylabel('\phi_{i,R} [ ]')
grid on
ylim([-0.1 0.11])

pdense=plot([1000 1000],[-.2 1.3],'linewidth',2,'color','#66101F');
set(gca,'fontsize',10,'xtick',[0:200:1600],'ytick',[-.1:.05:.1])


xlim([0 1650])
%normalized by offshore AWAC
% errorbar(distance_cumulative_with_awac(:,1),flip([1 nanmean(Hrms_dense10./AWAC.wav.Hm0)]),...
%    flip([0 nanstd(Hrms_dense10./AWAC.wav.Hm0)]),'b','linewidth',1.4)
% errorbar(distance_cumulative_with_awac(:,1),flip([1 nanmean(h_modeled_spec(2:end,:))]),...
%     flip([0 nanstd(h_modeled_spec(2:end,:))]),'color',[39 45 45]./256,'linewidth',1.4)

%normalized by offshore RBR
% errorbar(distance_cumulative(:,1),flip(nanmean(Hrms_dense10./Hrms_dense10(:,1))),flip(nanstd(Hrms_dense10./AWAC.wav.Hm0)),'b','linewidth',1.4)
% errorbar(distance_cumulative(:,1),flip(nanmean(h_modeled_spec./Hrms_dense10(:,1))),flip(nanstd(h_modeled_spec)),'color',[39 45 45]./256,'linewidth',1.4)

%normalized by previous station
%%
s1.Position = [0.1    0.70    0.824    0.28];
s2.Position = [0.1    0.39    0.824    0.28];
s3.Position = [0.1    0.09    0.824    0.28];

t=annotation('textbox','color','w','Position',[0.01 0.94 0.027 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(a)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle')%)%
t=annotation('textbox','color','w','Position',[0.01 0.6323 0.027 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(b)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle')%)%
t=annotation('textbox','color','w','Position',[0.01 0.32 0.027 0.04], 'BackgroundColor',[100 147 167]./255,...
    'EdgeColor',[100 147 167]./255,'String','(c)','FitBoxToText','off','fontsize',10,'horizontalalignment','center','verticalalignment','middle')%)%

l.Position = [ 0.6045    0.3955    0.3159    0.083];
%%
% print(f5,'C:\Users\NeryNeto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\article\chapter2\Figures\Figure5_waveattenuation.tiff','-dtiff','-r300')
print(f5,'C:\Users\NeryNeto\OneDrive - Nortek AS\Documents\UWA backup\PhD\chapters\Rottnest\article\chapter2\Figures\Figure5_reveiwtest.tiff','-dtiff','-r300')

% 