%% Read, filter, calculate Urms values and save Sawhorse instruments

% Sawhorse:
%     AQD 11266 (11 to 29/03)
%     ADV 6314 (11 to 29)
%     ADV 6333 (11 to 30)
%     OBS 2997 (20 to 28)

%% Sawhorse instruments - create NetCDF files
path = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\';
% load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\Sawhorse_instruments.mat')

%% Read NetCDF files
ADV6314 = ncreadall (fullfile(path,'\VEC6314\ADV6314_V0'));
ADV6333 = ncreadall (fullfile(path,'\VEC6333\ADV6333_V0'));

fn = fieldnames(ADV6333); %the instrument presented an error that made it keeping measuring after 7200 bursts
for k=1:numel(fn)
    if (size(ADV6333.(fn{k}),2) > 7200)
        ADV6333.(fn{k})(:,7201:end)=[];
    end
end

ADV6333.VEL_X (abs(ADV6333.VEL_Y)>.6) = NaN; ADV6333.VEL_Y (abs(ADV6333.VEL_Y)>.8) = NaN;
ADV6314.VEL_X (abs(ADV6314.VEL_X)>.6) = NaN; ADV6314.VEL_Y (abs(ADV6314.VEL_Y)>.8) = NaN;

for i = 1:size(ADV6333.VEL_X, 1)
    if i <= size(ADV6314.VEL_X,1)
        [ADV6314.Despike_vel_x(i,1:3600), ADV6314.Despike_vel_y(i,1:3600), ADV6314.Despike_vel_z(i,1:3600)] = ...
            func_despike_phasespace3d_3var (ADV6314.VEL_X(i,1:3600),...
            ADV6314.VEL_Y(i,1:3600), ADV6314.VEL_Z(i,1:3600), 0);
        
        [ADV6333.Despike_vel_x(i,1:3600), ADV6333.Despike_vel_y(i,1:3600), ADV6333.Despike_vel_z(i,1:3600)] = ...
            func_despike_phasespace3d_3var (ADV6333.VEL_X(i,1:3600),...
            ADV6333.VEL_Y(i,1:3600), ADV6333.VEL_Z(i,1:3600), 0);
        
        [ADV6314.Despike_vel_x(i,3601:7200), ADV6314.Despike_vel_y(i,3601:7200), ADV6314.Despike_vel_z(i,3601:7200)] = ...
            func_despike_phasespace3d_3var (ADV6314.VEL_X(i,3601:7200),...
            ADV6314.VEL_Y(i,3601:7200), ADV6314.VEL_Z(i,3601:7200), 0);
        
        [ADV6333.Despike_vel_x(i,3601:7200), ADV6333.Despike_vel_y(i,3601:7200), ADV6333.Despike_vel_z(i,3601:7200)] = ...
            func_despike_phasespace3d_3var (ADV6333.VEL_X(i,3601:7200),...
            ADV6333.VEL_Y(i,3601:7200), ADV6333.VEL_Z(i,3601:7200), 0);
    else
        [ADV6333.Despike_vel_x(i,1:3600), ADV6333.Despike_vel_y(i,1:3600), ADV6333.Despike_vel_z(i,1:3600)] = ...
            func_despike_phasespace3d_3var (ADV6333.VEL_X(i,1:3600),...
            ADV6333.VEL_Y(i,1:3600), ADV6333.VEL_Z(i,1:3600), 0);
        
        [ADV6333.Despike_vel_x(i,3601:7200), ADV6333.Despike_vel_y(i,3601:7200), ADV6333.Despike_vel_z(i,3601:7200)] = ...
            func_despike_phasespace3d_3var (ADV6333.VEL_X(i,3601:7200),...
            ADV6333.VEL_Y(i,3601:7200), ADV6333.VEL_Z(i,3601:7200), 0);
    end
    
end


[vel_maja, vel_mina , Current_Meana,...
    Current_thetaa, Wave_thetaa, ...
    Despike_cur_maja, Despike_cur_mina,...
    ADV6314.fx, ubra, ubr_specta,...
    Suua, Svva] = rotate_ADV_max_variance (ADV6314.VEL_X (:,1:3600),...
    ADV6314.VEL_Y (:,1:3600),...
    ADV6314.Despike_vel_x (:,1:3600), ADV6314.Despike_vel_y (:,1:3600),...
    ADV6314.CORR_X (:,1:3600), ADV6314.CORR_Y (:,1:3600),50,...
    4);

[vel_majb, vel_minb , Current_Meanb,...
    Current_thetab, Wave_thetab, ...
    Despike_cur_majb, Despike_cur_minb,...
    ADV6314.fx, ubrb, ubr_spectb,...
    Suub, Svvb] = rotate_ADV_max_variance (ADV6314.VEL_X (:,3601:7200),...
    ADV6314.VEL_Y (:,3601:7200),...
    ADV6314.Despike_vel_x (:,3601:7200), ADV6314.Despike_vel_y (:,3601:7200),...
    ADV6314.CORR_X (:,3601:7200), ADV6314.CORR_Y (:,3601:7200),50,...
    4);

[vel_majc, vel_minc , Current_Meanc,...
    Current_thetac, Wave_thetac, ...
    Despike_cur_majc, Despike_cur_minc,...
    ADV6333.fx, ubrc, ubr_spectc,...
    Suuc, Svvc] = rotate_ADV_max_variance (ADV6333.VEL_X (:,...
    1:3600), ADV6333.VEL_Y (:,1:3600),...
    ADV6333.Despike_vel_x (:,1:3600), ADV6333.Despike_vel_y (:,1:3600),...
    ADV6333.CORR_X (:,1:3600), ADV6333.CORR_Y (:,1:3600),50,...
    4);

[vel_majd, vel_mind , Current_Meand,...
    Current_thetad, Wave_thetad, ...
    Despike_cur_majd, Despike_cur_mind,...
    ADV6333.fx, ubrd, ubr_spectd,...
    Suud, Svvd] = rotate_ADV_max_variance (ADV6333.VEL_X (:,...
    3601:7200), ADV6333.VEL_Y (:,3601:7200),...
    ADV6333.Despike_vel_x (:,3601:7200), ADV6333.Despike_vel_y (:,3601:7200),...
    ADV6333.CORR_X (:,3601:7200), ADV6333.CORR_Y (:,3601:7200),50,...
    4);

%Calculate Reynold Stresses
for i=1:size(vel_majd,1)
    % Benilov
    
    [ADV6333_Benilov_upup_bara(i), ADV6333_Benilov_vpvp_bara(i), ADV6333_Benilov_wpwp_bara(i),...
        ADV6333_Benilov_upwp_bara(i), ADV6333_Benilov_uwave_wwave_bara(i)] = calculate_Benilov (...
        vel_majc(i,:), vel_minc(i,:),...
        ADV6333.Despike_vel_z(i,1:3600),...
        ADV6333.PRES(i,1:3600),...
        0.608,...
        4,...
        1,...
        1026);
    
    [ADV6333_Benilov_upup_barb(i), ADV6333_Benilov_vpvp_barb(i), ADV6333_Benilov_wpwp_barb(i),...
        ADV6333_Benilov_upwp_barb(i), ADV6333_Benilov_uwave_wwave_barb(i)] = calculate_Benilov (...
        vel_majd(i,:), vel_mind(i,:),...
        ADV6333.Despike_vel_z(i,3601:end),...
        ADV6333.PRES(i,3601:end),...
        0.608,...
        4,...
        1,...
        1026);
    
    %Bricker and Monismith
    [ADV6333_BM_upupa(i,:), ADV6333_BM_vpvpa(i,:), ADV6333_BM_wpwpa(i,:),...
        ADV6333_BM_upwpa(i,:) ADV6333_BM_uw_wavea(i,:)] = calculate_BrickerMonismith (vel_majc(i,:)',...
        vel_minc(i,:)', ADV6333.Despike_vel_z(i,1:3600)',...
        ADV6333.PRES(i,1:3600)',...
        4);
    
    [ADV6333_BM_upupb(i,:), ADV6333_BM_vpvpb(i,:), ADV6333_BM_wpwpb(i,:),...
        ADV6333_BM_upwpb(i,:) ADV6333_BM_uw_waveb(i,:)] = calculate_BrickerMonismith (vel_majd(i,:)',...
        vel_mind(i,:)', ADV6333.Despike_vel_z(i,3601:end)',...
        ADV6333.PRES(i,3601:end)',...
        4);
    
    %Young and Webster
    [ADV6333_YW_uua(i), ADV6333_YW_vva(i), ADV6333_YW_wwa(i), ADV6333_YW_uwa(i), ADV6333_YW_uw_wavea(i)] = ...
        calculate_youngWebster (vel_majc(i,:)',...
        vel_minc(i,:)', ADV6333.Despike_vel_z(i,1:3600)',...
        ADV6333.PRES(i,1:3600),...
        4);
    
    [ADV6333_YW_uub(i), ADV6333_YW_vvb(i), ADV6333_YW_wwb(i), ADV6333_YW_uwb(i), ADV6333_YW_uw_waveb(i)] = ...
        calculate_youngWebster (vel_majd(i,:)',...
        vel_mind(i,:)', ADV6333.Despike_vel_z(i,3601:end)',...
        ADV6333.PRES(i,3601:end),...
        4);
    
    if i<441
        
           [ADV6314_Benilov_upup_bara(i), ADV6314_Benilov_vpvp_bara(i), ADV6314_Benilov_wpwp_bara(i),...
               ADV6314_Benilov_upwp_bara(i), ADV6314_Benilov_uwave_wwave_bara(i)] = calculate_Benilov (...
               vel_maja(i,1:3600),vel_mina(i,1:3600),...
                ADV6314.Despike_vel_z(i,1:3600),...
                ADV6314.PRES(i,1:3600),...
                0.235,...
                4,...
                1,...
                1026);
        
        
           [ADV6314_Benilov_upup_barb(i), ADV6314_Benilov_vpvp_barb(i), ADV6314_Benilov_wpwp_barb(i),...
               ADV6314_Benilov_upwp_barb(i), ADV6314_Benilov_uwave_wwave_barb(i)] = calculate_Benilov (...
               vel_majb(i,1:3600),vel_minb(i,1:3600),...
                ADV6314.Despike_vel_z(i,1:3600),...
                ADV6314.PRES(i,1:3600),...
                0.235,...
                4,...
                1,...
                1026);
        
        [ADV6314_BM_upupa(i,:), ADV6314_BM_vpvpa(i,:), ADV6314_BM_wpwpa(i,:),...
            ADV6314_BM_upwpa(i,:) ADV6314_BM_uw_wavea(i,:)] = calculate_BrickerMonismith (vel_maja(i,:)',...
                vel_mina(i,:)', ADV6314.Despike_vel_z(i,1:3600)',...
                ADV6314.PRES(i,1:3600)',...
                4);
        
            [ADV6314_BM_upupb(i,:), ADV6314_BM_vpvpb(i,:), ADV6314_BM_wpwpb(i,:),...
                ADV6314_BM_upwpb(i,:) ADV6314_BM_uw_waveb(i,:)] = calculate_BrickerMonismith (vel_majb(i,:)',...
                vel_minb(i,:)', ADV6314.Despike_vel_z(i,3601:end)',...
                ADV6314.PRES(i,3601:end)',...
                4);
        
        [ADV6314_YW_uua(i), ADV6314_YW_vva(i), ADV6314_YW_wwa(i), ADV6314_YW_uwa(i), ADV6314_YW_uw_wavea(i)] = ...
            calculate_youngWebster (vel_maja(i,:)',...
            vel_mina (i,:)', ADV6314.Despike_vel_z(i,1:3600)',...
            ADV6314.PRES(i,1:3600),...
            4);
        
        [ADV6314_YW_uub(i), ADV6314_YW_vvb(i), ADV6314_YW_wwb(i), ADV6314_YW_uwb(i), ADV6314_YW_uw_waveb(i)] = ...
            calculate_youngWebster (vel_majb(i,:)',...
            vel_minb(i,:)', ADV6314.Despike_vel_z(i,3601:end)',...
            ADV6314.PRES(i,3601:end),...
            4);
        
        %Shaw and Trowbridge
%         [ADV6314_uw(i), ADV6314_vw(i),] = calculate_turbulentStress_ShawTrowbridge(...
%             [vel_maja(i,:)',vel_mina(i,:)',ADV6314.Despike_vel_z(i,1:3600)'],...
%             [vel_majc(i,:)',vel_minc(i,:)',ADV6333.Despike_vel_z(i,1:3600)'],...
%             3599)
    end
    
end

% Read and create waves
for i = 1:size(ADV6333.PRES,1)
    
    aux=nanmean(ADV6333.PRES(i,1:3600));
    [a1.Eta(i,:)]=PcorFFTFun(ADV6333.PRES(i,1:3600)',4,900,1024,...
        aux+1,1,1/25,1/2.5,1,'all','off','off');
    [Hm0.a1(i,1),Tm01.a1(i,1),Tm02.a1(i,1),Tp.a1(i,1),fp.ba1(i,1),f.a1(i,:),SA_1(i,:)]=WaveSpectraFun((a1.Eta(i,:))',4,900,...
        2^10,aux,0.608+0.214,.1,.8,1,-3,'off','off','off','off');
    
    aux2=nanmean(ADV6333.PRES(i,3601:end));
    [a2.Eta(i,:)]=PcorFFTFun(ADV6333.PRES(i,3601:end)',4,900,1024,...
        aux2+1,1,1/25,1/2.5,1,'all','off','off');
    [Hm0.a2(i,1),Tm01.a2(i,1),Tm02.a2(i,1),Tp.a2(i,1),fp.a2(i,1),f.a2(i,:),SA_2(i,:)]=WaveSpectraFun((a2.Eta(i,:))',4,900,...
        2^10,aux,0.608+0.214,.1,.8,1,-3,'off','off','off','off'); clear aux aux2
    
    if i<441
        aux=nanmean(ADV6314.PRES(i,1:3600));
        [b1.Eta(i,:)]=PcorFFTFun(ADV6314.PRES(i,1:3600)',4,900,1024,...
            aux+1,1,1/25,1/2.5,1,'all','off','off');
        [Hm0.b1(i,1),Tm01.b1(i,1),Tm02.b1(i,1),Tp.b1(i,1),fp.bb1(i,1),f.b1(i,:),SB_1(i,:)]=WaveSpectraFun((b1.Eta(i,:))',4,900,...
            2^10,aux,0.608+0.214,.1,.8,1,-3,'off','off','off','off');
        
        aux2=nanmean(ADV6314.PRES(i,3601:end));
        [b2.Eta(i,:)]=PcorFFTFun(ADV6314.PRES(i,3601:end)',4,900,1024,...
            aux2+1,1,1/25,1/2.5,1,'all','off','off');
        [Hm0.b2(i,1),Tm01.b2(i,1),Tm02.b2(i,1),Tp.b2(i,1),fp.b2(i,1),f.b2(i,:),SB_2(i,:)]=WaveSpectraFun((b2.Eta(i,:))',4,900,...
            2^10,aux,0.608+0.214,.1,.8,1,-3,'off','off','off','off'); clear aux aux2
    end
end


%organise currents
ADV6314.vel_maj = [vel_maja vel_majb]; ADV6314.vel_min = [vel_mina vel_minb];
ADV6333.vel_maj = [vel_majc vel_majd]; ADV6333.vel_min = [vel_minc vel_mind];
ADV6314.Despike_cur_maj = [Despike_cur_maja Despike_cur_majb];
ADV6314.Despike_cur_min = [Despike_cur_mina Despike_cur_minb];
ADV6333.Despike_cur_maj = [Despike_cur_majc Despike_cur_majd];
ADV6333.Despike_cur_min = [Despike_cur_minc Despike_cur_mind];

for i = 1:2:size(ADV6333.VEL_X, 1)*2
    
    if i <= size(ADV6314.VEL_X,1)*2
        ADV6314.Current_Mean (i,1) = Current_Meana ((i+1)/2);
        ADV6314.Current_Mean (i+1,1) = Current_Meanb ((i+1)/2);
        ADV6314.Current_theta (i,1) = wrapTo360( Current_thetaa ((i+1)/2));
        ADV6314.Current_theta (i+1,1) = wrapTo360(Current_thetab ((i+1)/2));
        ADV6314.Wave_theta (i) = wrapTo360(Wave_thetaa ((i+1)/2));
        ADV6314.Wave_theta (i+1) = wrapTo360(Wave_thetab ((i+1)/2));
        ADV6314.ubr (i) = ubra ((i+1)/2);
        ADV6314.ubr (i+1) = ubrb ((i+1)/2);
        ADV6314.ubr_spect (i) = ubr_specta ((i+1)/2);
        ADV6314.ubr_spect (i+1) = ubr_spectb ((i+1)/2);
        ADV6314.Suu (i,:) = Suua ((i+1)/2,:);
        ADV6314.Suu (i+1,:) = Suub ((i+1)/2,:);
        ADV6314.Svv (i,:) = Svva ((i+1)/2,:);
        ADV6314.Svv (i+1,:) = Svvb ((i+1)/2,:);
        ADV6314.time_ft (i) = ADV6314.TIME((i+1)/2,1);
        ADV6314.time_ft (i+1) = ADV6314.TIME ((i+1)/2,3601);
        
        ADV6314.Benilov.upup_bar (i) = ADV6314_Benilov_upup_bara ((i+1)/2);
        ADV6314.Benilov.upup_bar (i+1) = ADV6314_Benilov_upup_barb ((i+1)/2);              
        ADV6314.Benilov.vpvp_bar (i) = ADV6314_Benilov_vpvp_bara ((i+1)/2);
        ADV6314.Benilov.vpvp_bar (i+1) = ADV6314_Benilov_vpvp_barb ((i+1)/2);    
        ADV6314.Benilov.wpwp_bar (i) = ADV6314_Benilov_wpwp_bara ((i+1)/2);
        ADV6314.Benilov.wpwp_bar (i+1) = ADV6314_Benilov_wpwp_barb ((i+1)/2);
        ADV6314.Benilov.upwp_bar (i) = ADV6314_Benilov_upwp_bara ((i+1)/2);
        ADV6314.Benilov.upwp_bar (i+1) = ADV6314_Benilov_upwp_barb ((i+1)/2);  
        
        ADV6314.BrickerMonismith.upup_bar (i) = ADV6314_BM_upupa ((i+1)/2);
        ADV6314.BrickerMonismith.upup_bar (i+1) = ADV6314_BM_upupb ((i+1)/2);
        ADV6314.BrickerMonismith.vpvp_bar (i) = ADV6314_BM_vpvpa ((i+1)/2);
        ADV6314.BrickerMonismith.vpvp_bar (i+1) = ADV6314_BM_vpvpb ((i+1)/2);
        ADV6314.BrickerMonismith.wpwp_bar (i) = ADV6314_BM_wpwpa ((i+1)/2);
        ADV6314.BrickerMonismith.wpwp_bar (i+1) = ADV6314_BM_wpwpb ((i+1)/2);
        ADV6314.BrickerMonismith.upwp_bar (i) = ADV6314_BM_upwpa ((i+1)/2);
        ADV6314.BrickerMonismith.upwp_bar (i+1) = ADV6314_BM_upwpb ((i+1)/2);
        ADV6314.BrickerMonismith.uwave_wwave_bar (i) = ADV6314_BM_uw_wavea((i+1)/2);
        ADV6314.BrickerMonismith.uwave_wwave_bar (i+1) = ADV6314_BM_uw_waveb((i+1)/2);
        
        ADV6314.YoungWebster.upup_bar (i) = ADV6314_YW_uua ((i+1)/2);
        ADV6314.YoungWebster.upup_bar (i+1) = ADV6314_YW_uub ((i+1)/2);
        ADV6314.YoungWebster.vpvp_bar (i) = ADV6314_YW_vva ((i+1)/2);
        ADV6314.YoungWebster.vpvp_bar (i+1) = ADV6314_YW_vvb ((i+1)/2);
        ADV6314.YoungWebster.wpwp_bar (i) = ADV6314_YW_wwa ((i+1)/2);
        ADV6314.YoungWebster.wpwp_bar (i+1) = ADV6314_YW_wwb ((i+1)/2);
        ADV6314.YoungWebster.upwp_bar (i) = ADV6314_YW_uwa ((i+1)/2);
        ADV6314.YoungWebster.upwp_bar (i+1) = ADV6314_YW_uwb ((i+1)/2);
        ADV6314.YoungWebster.uwave_wwave_bar (i) = ADV6314_YW_uw_wavea((i+1)/2);
        ADV6314.YoungWebster.uwave_wwave_bar (i+1) = ADV6314_YW_uw_waveb((i+1)/2);
        
        ADV6333.Current_Mean (i,1) = Current_Meanc ((i+1)/2);
        ADV6333.Current_Mean (i+1,1) = Current_Meand ((i+1)/2);
        ADV6333.Current_theta (i,1) = wrapTo360(Current_thetac ((i+1)/2));
        ADV6333.Current_theta (i+1,1) = wrapTo360(Current_thetad ((i+1)/2));
        ADV6333.Wave_thetc (i) = wrapTo360(Wave_thetac ((i+1)/2));
        ADV6333.Wave_thetc (i+1) = wrapTo360(Wave_thetad ((i+1)/2));
        ADV6333.ubr (i) = ubrc ((i+1)/2);
        ADV6333.ubr (i+1) = ubrd ((i+1)/2);
        ADV6333.ubr_spect (i) = ubr_spectc ((i+1)/2);
        ADV6333.ubr_spect (i+1) = ubr_spectd ((i+1)/2);
        ADV6333.Suu (i,:) = Suuc ((i+1)/2,:);
        ADV6333.Suu (i+1,:) = Suud ((i+1)/2,:);
        ADV6333.Svv (i,:) = Svvc ((i+1)/2,:);
        ADV6333.Svv (i+1,:) = Svvd ((i+1)/2,:);
        ADV6333.time_ft (i) = ADV6333.TIME((i+1)/2,1);
        ADV6333.time_ft (i+1) = ADV6333.TIME ((i+1)/2,3601);
        
        ADV6333.Benilov.upup_bar (i) = ADV6333_Benilov_upup_bara ((i+1)/2);
        ADV6333.Benilov.upup_bar (i+1) = ADV6333_Benilov_upup_barb ((i+1)/2);
        ADV6333.Benilov.vpvp_bar (i) = ADV6333_Benilov_vpvp_bara ((i+1)/2);
        ADV6333.Benilov.vpvp_bar (i+1) = ADV6333_Benilov_vpvp_barb ((i+1)/2);
        ADV6333.Benilov.wpwp_bar (i) = ADV6333_Benilov_wpwp_bara ((i+1)/2);
        ADV6333.Benilov.wpwp_bar (i+1) = ADV6333_Benilov_wpwp_barb ((i+1)/2);
        ADV6333.Benilov.upwp_bar (i) = ADV6333_Benilov_upwp_bara ((i+1)/2);
        ADV6333.Benilov.upwp_bar (i+1) = ADV6333_Benilov_upwp_barb ((i+1)/2);
        ADV6333.Benilov.uwave_wwave_bar (i) = ADV6333_Benilov_uwave_wwave_bara((i+1)/2);
        ADV6333.Benilov.uwave_wwave_bar (i+1) = ADV6333_Benilov_uwave_wwave_barb((i+1)/2);
        
        ADV6333.BrickerMonismith.upup_bar (i) = ADV6333_BM_upupa ((i+1)/2);
        ADV6333.BrickerMonismith.upup_bar (i+1) = ADV6333_BM_upupb ((i+1)/2);
        ADV6333.BrickerMonismith.vpvp_bar (i) = ADV6333_BM_vpvpa ((i+1)/2);
        ADV6333.BrickerMonismith.vpvp_bar (i+1) = ADV6333_BM_vpvpb ((i+1)/2);
        ADV6333.BrickerMonismith.wpwp_bar (i) = ADV6333_BM_wpwpa ((i+1)/2);
        ADV6333.BrickerMonismith.wpwp_bar (i+1) = ADV6333_BM_wpwpb ((i+1)/2);
        ADV6333.BrickerMonismith.upwp_bar (i) = ADV6333_BM_upwpa ((i+1)/2);
        ADV6333.BrickerMonismith.upwp_bar (i+1) = ADV6333_BM_upwpb ((i+1)/2);
        ADV6333.BrickerMonismith.uwave_wwave_bar (i) = ADV6333_BM_uw_wavea((i+1)/2);
        ADV6333.BrickerMonismith.uwave_wwave_bar (i+1) = ADV6333_BM_uw_waveb((i+1)/2);
        
        ADV6333.YoungWebster.upup_bar (i) = ADV6333_YW_uua ((i+1)/2);
        ADV6333.YoungWebster.upup_bar (i+1) = ADV6333_YW_uub ((i+1)/2);
        ADV6333.YoungWebster.vpvp_bar (i) = ADV6333_YW_vva ((i+1)/2);
        ADV6333.YoungWebster.vpvp_bar (i+1) = ADV6333_YW_vvb ((i+1)/2);
        ADV6333.YoungWebster.wpwp_bar (i) = ADV6333_YW_wwa ((i+1)/2);
        ADV6333.YoungWebster.wpwp_bar (i+1) = ADV6333_YW_wwb ((i+1)/2); 
        ADV6333.YoungWebster.upwp_bar (i) = ADV6333_YW_uwa ((i+1)/2);
        ADV6333.YoungWebster.upwp_bar (i+1) = ADV6333_YW_uwb ((i+1)/2);
        ADV6333.YoungWebster.uwave_wwave_bar (i) = ADV6333_YW_uw_wavea((i+1)/2);
        ADV6333.YoungWebster.uwave_wwave_bar (i+1) = ADV6333_YW_uw_waveb((i+1)/2);
        
    else
        ADV6333.Current_Mean (i,1) = Current_Meanc ((i+1)/2);
        ADV6333.Current_Mean (i+1,1) = Current_Meand ((i+1)/2);
        ADV6333.Current_thetc (i,1) = wrapTo360(Current_thetac ((i+1)/2));
        ADV6333.Current_thetc (i+1,1) = wrapTo360(Current_thetad ((i+1)/2));
        ADV6333.Wave_thetc (i) = wrapTo360(Wave_thetac ((i+1)/2));
        ADV6333.Wave_thetc (i+1) = wrapTo360(Wave_thetad ((i+1)/2));
        ADV6333.ubr (i) = ubrc (1,(i+1)/2);
        ADV6333.ubr (i+1) = ubrd (1,(i+1)/2);
        ADV6333.ubr_spect (i) = ubr_spectc (1,(i+1)/2);
        ADV6333.ubr_spect (i+1) = ubr_spectd (1,(i+1)/2);
        ADV6333.Suu (i,:) = Suuc ((i+1)/2,:);
        ADV6333.Suu (i+1,:) = Suud ((i+1)/2,:);
        ADV6333.Svv (i,:) = Svvc ((i+1)/2,:);
        ADV6333.Svv (i+1,:) = Svvd ((i+1)/2,:);
        ADV6333.time_ft (i) = ADV6333.TIME((i+1)/2,1);
        ADV6333.time_ft (i+1) = ADV6333.TIME ((i+1)/2,3601);
        
        ADV6333.Benilov.upup_bar (i) = ADV6333_Benilov_upup_bara ((i+1)/2);
        ADV6333.Benilov.upup_bar (i+1) = ADV6333_Benilov_upup_barb ((i+1)/2);
        ADV6333.Benilov.vpvp_bar (i) = ADV6333_Benilov_vpvp_bara ((i+1)/2);
        ADV6333.Benilov.vpvp_bar (i+1) = ADV6333_Benilov_vpvp_barb ((i+1)/2);
        ADV6333.Benilov.wpwp_bar (i) = ADV6333_Benilov_wpwp_bara ((i+1)/2);
        ADV6333.Benilov.wpwp_bar (i+1) = ADV6333_Benilov_wpwp_barb ((i+1)/2);
        ADV6333.Benilov.upwp_bar (i) = ADV6333_Benilov_upwp_bara ((i+1)/2);
        ADV6333.Benilov.upwp_bar (i+1) = ADV6333_Benilov_upwp_barb ((i+1)/2);
        ADV6333.Benilov.uwave_wwave_bar (i) = ADV6333_Benilov_uwave_wwave_bara((i+1)/2);
        ADV6333.Benilov.uwave_wwave_bar (i+1) = ADV6333_Benilov_uwave_wwave_barb((i+1)/2);
        
        ADV6333.BrickerMonismith.upup_bar (i) = ADV6333_BM_upupa ((i+1)/2);
        ADV6333.BrickerMonismith.upup_bar (i+1) = ADV6333_BM_upupb ((i+1)/2);
        ADV6333.BrickerMonismith.vpvp_bar (i) = ADV6333_BM_vpvpa ((i+1)/2);
        ADV6333.BrickerMonismith.vpvp_bar (i+1) = ADV6333_BM_vpvpb ((i+1)/2);
        ADV6333.BrickerMonismith.wpwp_bar (i) = ADV6333_BM_wpwpa ((i+1)/2);
        ADV6333.BrickerMonismith.wpwp_bar (i+1) = ADV6333_BM_wpwpb ((i+1)/2);
        ADV6333.BrickerMonismith.upwp_bar (i) = ADV6333_BM_upwpa ((i+1)/2);
        ADV6333.BrickerMonismith.upwp_bar (i+1) = ADV6333_BM_upwpb ((i+1)/2);
        ADV6333.BrickerMonismith.uwave_wwave_bar (i) = ADV6333_BM_uw_wavea((i+1)/2);
        ADV6333.BrickerMonismith.uwave_wwave_bar (i+1) = ADV6333_BM_uw_waveb((i+1)/2);
        
        ADV6333.YoungWebster.upup_bar (i) = ADV6333_YW_uua ((i+1)/2);
        ADV6333.YoungWebster.upup_bar (i+1) = ADV6333_YW_uub ((i+1)/2);
        ADV6333.YoungWebster.vpvp_bar (i) = ADV6333_YW_vva ((i+1)/2);
        ADV6333.YoungWebster.vpvp_bar (i+1) = ADV6333_YW_vvb ((i+1)/2);
        ADV6333.YoungWebster.wpwp_bar (i) = ADV6333_YW_wwa ((i+1)/2);
        ADV6333.YoungWebster.wpwp_bar (i+1) = ADV6333_YW_wwb ((i+1)/2); 
        ADV6333.YoungWebster.upwp_bar (i) = ADV6333_YW_uwa ((i+1)/2);
        ADV6333.YoungWebster.upwp_bar (i+1) = ADV6333_YW_uwb ((i+1)/2);
        ADV6333.YoungWebster.uwave_wwave_bar (i) = ADV6333_YW_uw_wavea((i+1)/2);
        ADV6333.YoungWebster.uwave_wwave_bar (i+1) = ADV6333_YW_uw_waveb((i+1)/2);
    end
end
%organise waves
for i=1:2:916
    ADV6333.Wave.Hm0(i,2) = ADV6333.TIME((i+1)/2,1);
    ADV6333.Wave.Hm0(i+1,2) = ADV6333.TIME((i+1)/2,3601);
    
    ADV6333.Wave.SM(i,:) = SA_1((i+1)/2,:);
    ADV6333.Wave.SM(i+1,:) = SA_2((i+1)/2,:);
    
    ADV6333.Wave.Hm0(i,1) = Hm0.a1((i+1)/2,1);
    ADV6333.Wave.Hm0(i+1,1) = Hm0.a2((i+1)/2,1);
    
    ADV6333.Wave.Tm01(i,1) = Tm01.a1((i+1)/2,1);
    ADV6333.Wave.Tm01(i+1,1) = Tm01.a2((i+1)/2,1);
    
    ADV6333.Wave.Tm02(i,1) = Tm02.a1((i+1)/2,1);
    ADV6333.Wave.Tm02(i+1,1) = Tm02.a2((i+1)/2,1);
    
    ADV6333.Wave.Tp(i,1) = Tp.a1((i+1)/2,1);
    ADV6333.Wave.Tp(i+1,1) = Tp.a2((i+1)/2,1);
    
    if i<880
        ADV6314.Wave.Hm0(i,2) = ADV6314.TIME((i+1)/2,1);
        ADV6314.Wave.Hm0(i+1,2) = ADV6314.TIME((i+1)/2,3601);
        
        ADV6314.Wave.Hm0(i,1) = Hm0.b1((i+1)/2,1);
        ADV6314.Wave.Hm0(i+1,1) = Hm0.b2((i+1)/2,1);
        
        ADV6314.Wave.SM_ADV_middle(i,:) = SB_1((i+1)/2,:);
        ADV6314.Wave.SM_ADV_middle(i+1,:) = SB_2((i+1)/2,:);
        
        ADV6314.Wave.Tm01(i,1) = Tm01.b1((i+1)/2,1);
        ADV6314.Wave.Tm01(i+1,1) = Tm01.b2((i+1)/2,1);
        
        ADV6314.Wave.Tm02(i,1) = Tm02.b1((i+1)/2,1);
        ADV6314.Wave.Tm02(i+1,1) = Tm02.b2((i+1)/2,1);
        
        ADV6314.Wave.Tp(i,1) = Tp.b1((i+1)/2,1);
        ADV6314.Wave.Tp(i+1,1) = Tp.b2((i+1)/2,1);
    end
end

ADV6314.Benilov.TKE = 0.5.*(ADV6314.Benilov.upup_bar.^2+ADV6314.Benilov.vpvp_bar.^2+ADV6314.Benilov.wpwp_bar.^2);
ADV6314.BrickerMonismith.TKE = 0.5.*(ADV6314.BrickerMonismith.upup_bar.^2+ADV6314.BrickerMonismith.vpvp_bar.^2+ADV6314.BrickerMonismith.wpwp_bar.^2);
ADV6314.YoungWebster.TKE = 0.5.*(ADV6314.YoungWebster.upup_bar.^2+ADV6314.YoungWebster.vpvp_bar.^2+ADV6314.YoungWebster.wpwp_bar.^2);

ADV6333.Benilov.TKE = 0.5.*(ADV6333.Benilov.upup_bar.^2+ADV6333.Benilov.vpvp_bar.^2+ADV6333.Benilov.wpwp_bar.^2);
ADV6333.BrickerMonismith.TKE = 0.5.*(ADV6333.BrickerMonismith.upup_bar.^2+ADV6333.BrickerMonismith.vpvp_bar.^2+ADV6333.BrickerMonismith.wpwp_bar.^2);
ADV6333.YoungWebster.TKE = 0.5.*(ADV6333.YoungWebster.upup_bar.^2+ADV6333.YoungWebster.vpvp_bar.^2+ADV6333.YoungWebster.wpwp_bar.^2);

clearvars -except ADV6314 ADV6333 caminho_VEC path
%% Sawhorse instruments - Aquadopp
% Get_Nortek_ADPHR2NetCDF ([path 'ADP11266'],'AQDSAW02','AQP11266_V0','',0,0,.1,.2,0,60,58)
AQP11266 = ncreadall ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\ADP11266\AQD11266_V0')

for i = 1:size(AQP11266.VEL_XE_RAW, 2)
    [AQP11266.Despike_vel_x(:,i), AQP11266.Despike_vel_y(:,i), AQP11266.Despike_vel_z(:,i)] = func_despike_phasespace3d_3var_ADP (AQP11266.VEL_XE_RAW(:,i),...
        AQP11266.VEL_YN_RAW(:,i), AQP11266.VEL_ZU_RAW(:,i), 0);
end
blank = 0.096; probe_height=1.09; cell_size=0.03;
AQP11266.z = probe_height-blank-cell_size*[1:32];

[AQP11266.vel_maj, AQP11266.vel_min, AQP11266.Mean_Current...
    AQP11266.Current_theta, AQP11266.Wave_theta,...
    AQP11266.Despike_vel_maj, AQP11266.Despike_vel_min,...
    AQP11266.freq_urms, AQP11266.ubr, AQP11266.ubj...
    AQP11266.Suu_ADP, AQP11266.Svv_ADP] = rotate_ADP_max_variance (...
    AQP11266.VEL_XE_RAW, AQP11266.VEL_YN_RAW, AQP11266.Despike_vel_x, AQP11266.Despike_vel_y);

samples_fft = 1024;
% Read and create waves
AQP11266.Wave.time = AQP11266.TIME_RAW(1:samples_fft:end); %AQP11266.Wave.time(end) = [];
AQP11266.time_ft = AQP11266.TIME_RAW(1:samples_fft/2:end); %AQP11266.time_ft(end) = [];

AQP11266.freestream = nanmean (AQP11266.Mean_Current(:,1:5),2);
BURSTS = floor(size(AQP11266.Despike_vel_x,1)/1024);
  
AQP11266.Despike_vel = sqrt((AQP11266.Despike_vel_x.^2)+(AQP11266.Despike_vel_y.^2));
heights = AQP11266.z(1:5);

% Waves
lf=.04;                 % Hz - low frequency cutoff
maxfac=200;             % Maximum value of factor scaling pressure to waves
minspec=0.03;           % m^2/Hz - minimum spectral level for computing direction and spreading
Ndir=0;                 %deg - direction offset (includes compass error and misalignment of cable probe relative to case
                        % the offset for the Aquadopp Profiler is 0

parms=[lf maxfac minspec Ndir];
dt=1;
hp=0.14;
hv=1.09;

for i=1:round((size((AQP11266.PRESSURE_RAW),1))/samples_fft)
        AQP11266.meanDepth_15min(i)=nanmean(AQP11266.PRESSURE_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft));
        AQP11266.time_15min(i) = nanmean(AQP11266.TIME_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft));

        
        aux2 = PcorFFTFun(AQP11266.PRESSURE_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft),1,samples_fft,256,...
            AQP11266.meanDepth_15min(i)+1.09,1.09,1/25,1/2.5,1,'all','off','off');
        
        [AQP11266.Wave.Hm0_ADP(i),...
            AQP11266.Wave.Tm01_ADP(i), AQP11266.Wave.Tm02_ADP(i), AQP11266.Wave.Tp_ADP(i),...
            AQP11266.Wave.fp_ADP(i),AQP11266.Wave.f_ADP(i,:), AQP11266.Wave.Spp_ADP(i,:)] = WaveSpectraFun(aux2,...
            1,samples_fft,256,AQP11266.meanDepth_15min(i)+1.09,1.09,1/25,1/2.5,1,-3,'off','off','off','off');
        
        clear aux aux2
        
        if sum (isnan(nanmean (AQP11266.Despike_vel_x(i*samples_fft - (samples_fft - 1):i*samples_fft,1:5),2))) < 100
         
        auxu = nanmean (AQP11266.Despike_vel_x(i*samples_fft - (samples_fft - 1):i*samples_fft,1:5),2)';
        if isnan(auxu)==1; auxu(1)=0; end; if isnan(auxu(end))==1; auxu(end)=0; end;
        X = ~isnan(auxu);
        Y = cumsum(X-diff([1 X])/2); vu = interp1(1:nnz(X),auxu(X),Y); clear auxu X Y
        
        auxv = nanmean (AQP11266.Despike_vel_y(i*samples_fft - (samples_fft - 1):i*samples_fft,1:5),2)';
        if isnan(auxv)==1; auxv(1)=0; end; if isnan(auxv(end))==1; auxv(end)=0; end;
        X = ~isnan(auxv);
        Y = cumsum(X-diff([1 X])/2); vv = interp1(1:nnz(X),auxv(X),Y); clear auxv X Y
        
        auxp = AQP11266.PRESSURE_RAW(i*samples_fft - (samples_fft - 1):i*samples_fft,1)';
        if isnan(auxp)==1; auxp(1)=0; end; if isnan(auxp(end))==1; auxp(end)=0; end;
        X = ~isnan(auxp);
        Y = cumsum(X-diff([1 X])/2); vp = interp1(1:nnz(X),auxp(X),Y); clear auxv X Y
        
        [AQP11266.Wave.Su(i,:),AQP11266.Wave.Sp(i,:),AQP11266.Wave.Dir(i,:),AQP11266.Wave.Spread,AQP11266.Wave.F,AQP11266.Wave.dF] = wds(vu',vv',vp',dt,100,hp,hv,parms);
        else
        AQP11266.Wave.Su(i,1:41)=nan;
        AQP11266.Wave.Sp(i,1:41)=nan;
        AQP11266.Wave.Dir(i,1:41)=nan;
        AQP11266.Wave.Spread(i,1:41)=nan;
        AQP11266.Wave.F(i,1:41)=nan;
        AQP11266.Wave.dF(i,1:41)=nan;
        end
        [AQP11266.Wave.ubr_sherwood(i),AQP11266.Wave.Tbr_sherwood(i),AQP11266.Wave.Snn_sherwood(:,i)]=ubspecdat(AQP11266.meanDepth_15min(i)+1.09,AQP11266.Wave.Spp_ADP(i,:),AQP11266.Wave.f_ADP(1,:));
        

               
end
   
AQP11266.Wave.Dir = wrapTo360 (AQP11266.Wave.Dir);
%% OBS
path_sawhorse_obs = 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Sawhorse\FLNTU'
file = fileread([path_sawhorse_obs '\OBS_2997_SAW_01.raw']);

lines = regexp (file, '\n', 'split');

[mmddyy,rest] = strtok(lines,char(9));
%Creates vector time
formatIn='m/dd/yy';
for i = 2:length(lines)-1
    try
        OBSday(i-1,:)=datevec(mmddyy(i),formatIn);
    catch
        OBSday(i-1,1:6)=nan;
    end
end

[hhmm,rest2] = strtok(rest,char(9));
formatIn2='HH:MM:SS';
for i = 2:length(lines)-1
    try
        OBShour(i-1,:)=datevec(hhmm(i),formatIn2);
    catch
        OBShour(i-1,:)=nan;
    end
end

OBS_2997(:,1)=datenum([OBSday(:,1:3) OBShour(:,4:6)]);OBS_2997(end)=[];

fid = fopen([path_sawhorse_obs '\OBS_2997_SAW_01.raw'],'r');
var_dir=textscan(fid, '%s %s %d %d %d %d %d','HeaderLines',1);
aux = var_dir{1,6}; aux(19478)=[];aux(63919)=[];
OBS_2997(:,2) = aux (1:end-1);
OBS_2997(:,3) = 0.0242*(OBS_2997(:,2)-50);

clearvars -except ADV6314 ADV6333 path AQP11266 OBS_2997 Wave_sawhorse

save ([path '\Sawhorse_instruments'],'-v7.3')
