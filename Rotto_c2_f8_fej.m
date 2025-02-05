load ('C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd')
%% Create flags
% checking if the flow is turbulent. Re<1000 should be excluded
b=0.014;u=RBR.ST9.ubr_LWT;
v=9.35.*10^(-7);
Re = b.*u./v;
flag_Re=(Re>1500);

for j=1:length(Re) 
% cumsum_st=cumsum(RBR.(ST).Syy(:,j).*100)./nansum(RBR.(ST).Syy(:,j));
% if cumsum_st(206)<80
% flag_swell(j)=1;
% else
% flag_swell(j)=0;
% end
end
flag_swell=RBR.ST9.ubr_LWT>0.1;
flag=and(flag_Re,flag_swell);

[~,min_cutoff] = min(abs(RBR.ST10.f_out(:,1)-0.05));
[~,max_cutoff] = min(abs(RBR.ST10.f_out(:,1)-0.2));
%%
figure (8)
hold on

fe78_mean=(RBR.ST8.fej_lowe(:,flag_swell)+RBR.ST7.fej_lowe(:,flag_swell))./2;
fej78_error=nanstd(fe78_mean,2,1)./sqrt(size(fe78_mean,1)/sqrt(155)); %155 represents effective degrees of freedom
% 
e=errorbar(RBR.ST9.f(min_cutoff:max_cutoff,1),nanmean(fe78_mean,2),fej78_error,'vertical','-sk','MarkerSize',8,...
    'MarkerEdgeColor','#e9c46a','MarkerFaceColor','#e9c46a','capsize',0);
e.Color = '#e9c46a';

% plot (RBR.ST9.f(min_cutoff:max_cutoff,1),nanmean(feh78_mean,2),'color','#2a9d8f','linewidth',1.5)
grid on
% title('$$f_{e,j}$$ between P7 - P9','Interpreter','latex','fontsize',14)
ylabel('$f_{e,j} [ \ ]$','Interpreter','latex','fontsize',14)
% plot (RBR.ST9.f(min_cutoff:max_cutoff,1),nanmean(RBR.ST7.fej_lowe(:,flag_swell),2),'color','#e9c46a','linewidth',1.5)
% title ('P-8')
grid on
xlabel('f [Hz]')
% ylabel('f_{e,j}')

xlim([0.04 0.21])
ylim([0 1.5])
box on
%%
print(gcf,'C:\Users\Neryneto\OneDrive - UWA\Documents\PhD\chapters\Rottnest\Plots\chapter2\Figure8_fej.tiff','-dtiff','-r300')
