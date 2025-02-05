load 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_RBRs_waves'
% load 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_RBRs'
load 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd'
load 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m\AWAC_6734.mat'
% cases=importdata ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\v2\tpar_awac.txt')
%%
rho=1023;
g=9.81;

% [a,n] = find(RBR.ST10.SB_cutoff(:,1));
% min_cutoff = min(a);max_cutoff = max(a);
[~,min_cutoff] = min(abs(RBR.ST10.f_out(:,1)-0.05));
[~,max_cutoff] = min(abs(RBR.ST10.f_out(:,1)-0.2));
% min_cutoff=65;max_cutoff=280;

deltaf = nanmean(diff(RBR.ST10.f_out(:,1)));
azymuth = 180-acosd((-6456977+6458205)./sqrt((363817-364383)^2+(6458205-6456977)^2));

RBR.ST9.lv = 0.43; RBR.ST9.av = 0.014; RBR.ST9.N = 542; RBR.ST9.modeled_theta = 147.5; RBR.ST9.distance = 104.2*cosd(abs(azymuth-RBR.ST9.modeled_theta)); 
RBR.ST8.lv = 0.45; RBR.ST8.av = 0.014; RBR.ST8.N = 542; RBR.ST8.modeled_theta = 142.5; RBR.ST8.distance = 100.18*cosd(azymuth-RBR.ST8.modeled_theta);
RBR.ST7.lv = 0.42; RBR.ST7.av = 0.014; RBR.ST7.N = 542; RBR.ST7.modeled_theta = 147.5; RBR.ST7.distance = 102.46*cosd(azymuth-RBR.ST7.modeled_theta); 
RBR.ST6.lv = 0.43; RBR.ST6.av = 0.014; RBR.ST6.N = 456; RBR.ST6.modeled_theta = 147.5; RBR.ST6.distance = 99.4*cosd(azymuth-RBR.ST6.modeled_theta);
RBR.ST5.lv = 0.45; RBR.ST5.av = 0.014; RBR.ST5.N = 456; RBR.ST5.modeled_theta = 152.5; RBR.ST5.distance = 193.55*cosd(azymuth-RBR.ST5.modeled_theta);
RBR.ST4.distance = 103.94; RBR.ST4.lv = 0.4; RBR.ST4.av = 0.014; RBR.ST4.N = 543;
RBR.ST3.distance = 92.35; RBR.ST3.lv = 0.3; RBR.ST3.av = 0.01; RBR.ST3.N = 543;
RBR.ST2.distance = 289.82; RBR.ST2.lv = 0.3; RBR.ST2.av = 0.01; RBR.ST2.N = 412;
RBR.ST1.distance = 267.65; RBR.ST1.lv = 0.07; RBR.ST1.av = 0.007;RBR.ST1.N = 100;
%% Between every station using the model with no friction using external AWAC as boundary condition

for k=1:10 %for each instrument
    ST = sprintf('ST%d',k);
    for i=1:length(RBR.(ST).meandepth) %for each time
        [RBR.(ST).L_m02(i,k),RBR.(ST).k_m02(i,k),RBR.(ST).sigma_m02(i,k)]=disper(RBR.(ST).meandepth(i),RBR.(ST).Tm01(i));
        [Lp(i,k),kp(i,k),sigmap(i,k)]=disper(RBR.(ST).meandepth(i),RBR.ST1.Tpswell(i));
        
for j = 1:size (RBR.(ST).S_out,1) %for each frequency
            [RBR.(ST).Lj(j,i),RBR.(ST).kj(j,i),RBR.(ST).sigmaj(j,i)]=disper(RBR.(ST).meandepth(i),1./RBR.(ST).f_out(j,1));
%             RBR.(ST).k_ry (j,i) = lindisp_explicit(RBR.(ST).f_out(j,1),RBR.(ST).meandepth(i));
        end
    end
end

for k=1:10 %for each instrument
    ST = sprintf('ST%d',k);
    %     clearvar RBR.(ST).fj
    RBR.(ST).Sj=[]; RBR.(ST).Hrms =[]; RBR.(ST).aj=[]; RBR.(ST).Ej=[];
    RBR.(ST).ubj_LWT=[];RBR.(ST).Cgj=[];    RBR.(ST).fej=[];RBR.(ST).fer=[];
    RBR.(ST).Fj=[];RBR.(ST).Hj=[];RBR.(ST).aj  =[];
    for i=1:length(RBR.(ST).meandepth) %for each time       
        
        RBR.(ST).Sj (:,i) = RBR.(ST).SB(min_cutoff:max_cutoff,i);

        RBR.(ST).Hrms (i) = sqrt(8).*sqrt(nansum(RBR.(ST).Sj(:,i).*deltaf));
        RBR.(ST).gama (i) = RBR.(ST).Hrms(i)./RBR.(ST).meandepth(i);
        RBR.(ST).Hj (:,i) = sqrt(8).*sqrt(RBR.(ST).Sj(:,i).*deltaf);
        RBR.(ST).arms (i) = sqrt(2.*nansum(RBR.(ST).Sj(:,i).*deltaf));
        RBR.(ST).aj (:,i) = sqrt(2.*RBR.(ST).Sj(:,i).*deltaf);
        
        omega=RBR.(ST).sigmaj(min_cutoff:max_cutoff,i);
        
        RBR.(ST).ubj_LWT(:,i) = (RBR.(ST).aj(:,i).*omega)./...
            (sinh(RBR.(ST).kj(min_cutoff:max_cutoff,i).*RBR.(ST).meandepth(i)));
        RBR.(ST).ubr_LWT(i) = sqrt(nansum(RBR.(ST).ubj_LWT(:,i).^2));
       
        RBR.(ST).Ej (:,i) = 0.5.*rho.*g.*(RBR.(ST).aj(:,i).^2);
        parj = (2.*RBR.(ST).kj(min_cutoff:max_cutoff,i).*RBR.(ST).meandepth(i))./(sinh(2.*RBR.(ST).kj(min_cutoff:max_cutoff,i).*RBR.(ST).meandepth(i)));
        RBR.(ST).Cgj (:,i)= 0.5.*(1+parj).*(omega./RBR.(ST).kj(min_cutoff:max_cutoff,i)); clear parj
        
        RBR.(ST).E (:,i) = 0.5.*rho.*g.*((RBR.(ST).Hrms_swell_out(i)/2).^2);
        par = (2.*RBR.(ST).k_m02(i,k).*RBR.(ST).meandepth(i))./(sinh(2.*RBR.(ST).k_m02(i,k).*RBR.(ST).meandepth(i)));
        RBR.(ST).Cg (:,i)= 0.5.*(1+par).*(RBR.(ST).sigma_m02(i,k)./RBR.(ST).k_m02(i,k)); clear par omega

            
%         if k<10 RBR.(ST).ks(:,i) = sqrt(RBR.(STf).Cg (:,i)./RBR.(STi).Cg (:,i)); end
    end
    RBR.(ST).omegaj = nansum(RBR.(ST).sigmaj(min_cutoff:max_cutoff,:).*...
        RBR.(ST).ubj_LWT.^2)./nansum(RBR.(ST).ubj_LWT.^2);

    RBR.(ST).Fj=RBR.(ST).Ej.*RBR.(ST).Cgj;
    
    RBR.(ST).Fr=RBR.(ST).E.*RBR.(ST).Cg;

end

%% Open model output
caminho = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awac_nofric';
caminho2 = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\v2';

rho=1026;
g=9.81;
for k=1:10 %for each instrument
    ST = sprintf('ST%d',k);
    Model.(ST).out = load(fullfile(caminho,['rbr' num2str(k) '_spc.dat']));
    depth = load(fullfile(caminho2,['rbr_points_' num2str(k) '.txt']));

    for i=1:length(Model.(ST).out)
        [Model.(ST).Lp(i),Model.(ST).kp(i),Model.(ST).sigmap(i)]=disper(depth(1,7),Model.(ST).out(i,4));
        Model.(ST).E (:,i) = 0.5.*rho.*g.*(((Model.(ST).out(i,2)./1.41)./2).^2);
        par = (2.*Model.(ST).kp(i).*depth(1,7))./(sinh(2.*Model.(ST).kp(i).*depth(1,7)));
        Model.(ST).Cg (:,i)= 0.5.*(1+par).*(Model.(ST).sigmap(i)./Model.(ST).kp(i)); clear omega par
    end
    
     Model.(ST).F =  Model.(ST).E .*Model.(ST).Cg;
     Model.(ST).F_interp = interp1(datenum(num2str(Model.(ST).out(:,1),'%f'),'yyyymmdd.HHMMSS'),Model.(ST).F, RBR.ST2.time_wave);

end

for k=1:9 %for each instrument
    ST1 = sprintf('ST%d',k);
    ST2 = sprintf('ST%d',k+1);
    
    kr(:,k) = (Model.(ST2).F./Model.(ST1).F)-1;
    if k>1 && k<5
        aux=kr(:,k);aux(aux>-0.01)=nan;
        kr(:,k)=aux;
    elseif k>=5 && k<9
        aux=kr(:,k);aux(aux>-0.05)=nan;
        kr(:,k)=aux;
    else
        aux=kr(:,k);aux(aux>0)=nan;
        kr(:,k)=aux;
    end
    rbr_awac_kr(:,k) = interp1(datenum(num2str(Model.(ST1).out(:,1),'%f'),'yyyymmdd.HHMMSS'),kr(:,k), RBR.ST2.time_wave);
end

% 
% for i=1:length(AWAC_6734.wave.Hm0)
%    [l c]=find(abs(AWAC_6734.wave.Hm0(i)-cases.data(:,2))<0.25);
%    [l2 c2]=find(abs(AWAC_6734.wave.Tp(i)-cases.data(:,3))<1);
%    [l3 c3]=find(abs(AWAC_6734.wave.DirTp(i)-cases.data(:,4))<10); 
%    
%    try
%    aux=intersect(intersect(l,l2),l3);
%    [a1 b]=find(min(abs(AWAC_6734.wave.DirTp(i)-cases.data(aux,4))));
%    a(i)=aux(b);
%    kr_aux (i) = kr(a(i),7);
%    catch
%    a(i)=nan;
%    kr_aux (i) = nan;
%    end 
% %    [~,rbr_awac_kr(i)]= min(abs(AWAC_6734.wave.time(i)-RBR.ST2.time_wave));
% end


% 
% rbr_awac_kr(rbr_awac_kr>-0.06)=-0.08;
%% Create flags
% checking if the flow is turbulent. Re<1000 should be excluded
b=0.014;u=RBR.ST8.ubr_LWT;
v=9.35.*10^(-7);
Re = b.*u./v;
% flag_aux=(Re>1000);
% 
% for j=1:length(Re) 
% cumsum_st=cumsum(RBR.ST9.Syy(:,j).*100)./nansum(RBR.ST9.Syy(:,j));
% if cumsum_st(206)<80
% flag_aux2(j)=1;
% else
% flag_aux2(j)=0;
% end
% end
% flag=and(flag_aux,flag_aux2);
d=0.014;
S=0.05;
hc=0.4;
lambdaP = (pi*d^2/4)/(S+d)^2;
% lambdaF_le = hc*d/((S+d)^2);
lambdaF_le = 0.15*d/(S+d)^2;

percentage_refraction = 0.11;

%from interpolation alpha factor
aux_frequencies = 1./[20:-.5:.5];
theoretical_alphaj=8.608.*aux_frequencies(1:end).^0.04522-7.082;
interp_alphaj=8.608.*RBR.ST8.f(min_cutoff:max_cutoff).^0.04522-7.082;
alphar=0.46;
for k=1:9 %for each instrument
    STi = sprintf('ST%d',k);
    STf = sprintf('ST%d',k+1);
    RBR.(STi).fej_lowe =[];RBR.(STi).fej_lowe1=[];RBR.(STi).fer_lowe=[];RBR.(STi).fer_lowe1=[];

    RBR.(STi).slope = nanmean((RBR.(STf).meandepth-RBR.(STi).meandepth)./RBR.(STi).distance);

    h = (RBR.(STf).meandepth + RBR.(STi).meandepth)./2;
    
    ubj_avg = (RBR.(STf).ubj_LWT + RBR.(STi).ubj_LWT)./2;
    ubr_avg = (RBR.(STf).ubr_LWT + RBR.(STi).ubr_LWT)./2;
    a = (RBR.(STf).arms + RBR.(STi).arms)./2;

    RBR.(STi).deltaj = -(RBR.(STf).Fj - ((rbr_awac_kr(:,k)' +1).*RBR.(STi).Fj))/RBR.(STi).distance;

    RBR.(STi).delta = -(RBR.(STf).Fr - ((rbr_awac_kr(:,k)' +1).*RBR.(STi).Fr))/RBR.(STi).distance;

    RBR.(STi).fe_jonsson = (3*pi.*RBR.(STi).delta)./(2*rho.*(ubr_avg.^3));  
    RBR.(STi).Cd_lowe = (4.*RBR.(STi).deltaj)./(lambdaF_le.*rho.*alphar.*ubr_avg.*((interp_alphaj'.*ubj_avg).^2));
%     RBR.(STi).Cd_lowe2 = (4.*RBR.(STi).delta)./(lambdaF_le.*rho.*0.8.*ubr_avg.^3);
%     RBR.(STi).Cd_lowe3 = RBR.(STi).delta./(lambdaF_le.*0.7.^3);
%     RBR.(STi).Cd_lowe4 = RBR.(STi).fe_jonsson./(lambdaF_le.*0.8.^3);

    RBR.(STi).fej_lowe = RBR.(STi).Cd_lowe.*lambdaF_le.*alphar.*(interp_alphaj'.^2);
%     RBR.(STi).fej_lowe = (4.*RBR.(STi).deltaj)./(rho.*ubr_avg.*(ubj_avg.^2));
    RBR.(STi).fer_lowe = nansum(RBR.(STi).fej_lowe.*ubj_avg.^2)./(nansum(ubj_avg.^2));
%     RBR.(STi).fer_lowe = nansum(RBR.(STi).fej_lowe1.*ubj_avg.^2)./(nansum(ubj_avg.^2));
end

%%
clearvars -except Model kr rbr_awac_kr RBR
save ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\RBR_with_Cd','-v7.3')