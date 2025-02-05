%% Forcing with AWAC
load 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_RBRs_waves_Lowe_JGR.mat'
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Aquadopps\AQD2985\AQD2985.mat')
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_8m\AWAC_cur_wav.mat')
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m\AWAC_6734.mat')
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Tripod\AQD11898\Tripod_instruments')

% load model outputs
awac_8m_model.NoFric = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcing\awac_8m.dat');
awac_8m_model.Veg_dense = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVeg\awac_8m.dat');
awac_8m_model.Veg_int = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegInterm\awac_8m.dat');
awac_8m_model.Veg_sparse3 = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegSparse3\awac_8m.dat');
awac_8m_model.Veg_sparse2 = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegSparse2\awac_8m.dat');
awac_8m_model.Veg_sparse1 = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegSparse1\awac_8m.dat');
awac_8m_model.Veg_fake_cd08 = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVegcd8\awac_8m.dat');
awac_8m_model.Veg_fake_cd05 = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVegcd5\awac_8m.dat');
awac_8m_model.Veg_fake_cd03 = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVeg\awac_8m.dat');
awac_8m_model.aqdawacforcing = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\aqdawac8forcing\awac_8mver.dat');
awac_8m_model.aqdforcing = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\aqdforcing\awac_8mver.dat');

awac_8m_model.NoFric(:,1) = datenum(num2str(awac_8m_model.NoFric(:,1),'%f'),'yyyymmdd.HHMMSS');
awac_8m_model.Veg_dense(:,1) = datenum(num2str(awac_8m_model.Veg_dense(:,1),'%f'),'yyyymmdd.HHMMSS');
awac_8m_model.Veg_int(:,1) = datenum(num2str(awac_8m_model.Veg_int(:,1),'%f'),'yyyymmdd.HHMMSS');
awac_8m_model.Veg_sparse3(:,1) = datenum(num2str(awac_8m_model.Veg_sparse3(:,1),'%f'),'yyyymmdd.HHMMSS');
awac_8m_model.Veg_sparse2(:,1) = datenum(num2str(awac_8m_model.Veg_sparse2(:,1),'%f'),'yyyymmdd.HHMMSS');
awac_8m_model.Veg_sparse1(:,1) = datenum(num2str(awac_8m_model.Veg_sparse1(:,1),'%f'),'yyyymmdd.HHMMSS');
awac_8m_model.Veg_fake_cd08(:,1) = datenum(num2str(awac_8m_model.Veg_fake_cd08(:,1),'%f'),'yyyymmdd.HHMMSS');
awac_8m_model.Veg_fake_cd05(:,1) = datenum(num2str(awac_8m_model.Veg_fake_cd05(:,1),'%f'),'yyyymmdd.HHMMSS');
awac_8m_model.Veg_fake_cd03(:,1) = datenum(num2str(awac_8m_model.Veg_fake_cd03(:,1),'%f'),'yyyymmdd.HHMMSS');
awac_8m_model.aqdawacforcing(:,1) = datenum(num2str(awac_8m_model.aqdawacforcing(:,1),'%f'),'yyyymmdd.HHMMSS');
awac_8m_model.aqdforcing(:,1) = datenum(num2str(awac_8m_model.aqdforcing(:,1),'%f'),'yyyymmdd.HHMMSS');

aqd_model.NoFric = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcing\aqd.dat');
aqd_model.Veg_dense = load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVeg\aqd.dat');
aqd_model.Veg_int = load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegInterm\aqd.dat');
aqd_model.Veg_sparse3 = load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegSparse3\aqd.dat');
aqd_model.Veg_sparse2 = load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegSparse2\aqd.dat');
aqd_model.Veg_sparse1 = load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegSparse1\aqd.dat');
aqd_model.Veg_fake_cd08 = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVegcd8\aqd.dat');
aqd_model.Veg_fake_cd05 = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVegcd5\aqd.dat');
aqd_model.Veg_fake_cd03 = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVeg\aqd1.dat');
aqd_model.Veg_fake_cd05 = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVegcd5\aqd.dat');
aqd_model.Veg_fake_cd03 = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVeg\aqd1.dat');
aqd_model.aqdawacforcing = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\aqdawac8forcing\aqdver.dat');
aqd_model.aqdforcing = load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\aqdforcing\aqdver.dat');

aqd_model.NoFric(:,1) = datenum(num2str(aqd_model.NoFric(:,1),'%f'),'yyyymmdd.HHMMSS');
aqd_model.Veg_dense(:,1) = datenum(num2str(aqd_model.Veg_dense(:,1),'%f'),'yyyymmdd.HHMMSS');
aqd_model.Veg_int(:,1) = datenum(num2str(aqd_model.Veg_int(:,1),'%f'),'yyyymmdd.HHMMSS');
aqd_model.Veg_sparse3(:,1) = datenum(num2str(aqd_model.Veg_sparse3(:,1),'%f'),'yyyymmdd.HHMMSS');
aqd_model.Veg_sparse2(:,1) = datenum(num2str(aqd_model.Veg_sparse2(:,1),'%f'),'yyyymmdd.HHMMSS');
aqd_model.Veg_sparse1(:,1) = datenum(num2str(aqd_model.Veg_sparse1(:,1),'%f'),'yyyymmdd.HHMMSS');
aqd_model.Veg_fake_cd08(:,1) = datenum(num2str(aqd_model.Veg_fake_cd08(:,1),'%f'),'yyyymmdd.HHMMSS');
aqd_model.Veg_fake_cd05(:,1) = datenum(num2str(aqd_model.Veg_fake_cd05(:,1),'%f'),'yyyymmdd.HHMMSS');
aqd_model.Veg_fake_cd03(:,1) = datenum(num2str(aqd_model.Veg_fake_cd03(:,1),'%f'),'yyyymmdd.HHMMSS');
aqd_model.aqdawacforcing(:,1) = datenum(num2str(aqd_model.aqdawacforcing(:,1),'%f'),'yyyymmdd.HHMMSS');
aqd_model.aqdforcing(:,1) = datenum(num2str(aqd_model.aqdforcing(:,1),'%f'),'yyyymmdd.HHMMSS');
%% Old model with bigger grid forcing with BoM buoy
% load 'C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_RBRs_waves_Lowe_JGR.mat'
% bom = readmatrix('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\RDW2019_Z.xlsx','Range','A2838:L5765');
% load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Aquadopps\AQD2985\AQD2985.mat')
% load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_8m\AWAC_cur_wav.mat')
% load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m\AWAC_6734.mat')
% % load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\bom_buoy.txt')
% bom(:,1)=bom(:,1)+693960;
% % load model outputs
% awac_15m_model=load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awac_15m.dat');
% awac_15m_model(:,1) = datenum(num2str(awac_15m_model(:,1),'%f'),'yyyymmdd.HHMMSS');
% 
% awac_8m_model=load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awac_8m.dat');
% awac_8m_model(:,1) = datenum(num2str(awac_8m_model(:,1),'%f'),'yyyymmdd.HHMMSS');
% 
% aqd_model=load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\aqd.dat');
% aqd_model(:,1) = datenum(num2str(aqd_model(:,1),'%f'),'yyyymmdd.HHMMSS');
%%
for i=1:10
    ST = sprintf('ST%d',i);
    results = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcing\rbr',...
        num2str(i) '.dat']);
    RBR_Model.nofric.(ST).time_wave = datenum(num2str(results(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.nofric.(ST).Hsig = results(:,2);
    RBR_Model.nofric.(ST).Hrms = results(:,2)./sqrt(2);
    RBR_Model.nofric.(ST).Tp = results(:,3);
    RBR_Model.nofric.(ST).MeanDir = results(:,4);
    RBR_Model.nofric.(ST).PeakDir = results(:,5);
    
    results = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcing_PM\rbr',...
        num2str(i) '_spc.dat']);
    RBR_Model.nofric_PM.(ST).time_wave = datenum(num2str(results(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.nofric_PM.(ST).Hsig = results(:,2);
    RBR_Model.nofric_PM.(ST).Hrms = results(:,2)./sqrt(2);
    RBR_Model.nofric_PM.(ST).Tp = results(:,3);
    RBR_Model.nofric_PM.(ST).MeanDir = results(:,4);
    RBR_Model.nofric_PM.(ST).PeakDir = results(:,5);
    
    results_vegetation = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVeg\rbr',...
        num2str(i) '.dat']);
    RBR_Model.vegetationDense.(ST).time_wave = datenum(num2str(results_vegetation(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.vegetationDense.(ST).Hsig = results_vegetation(:,2);
    RBR_Model.vegetationDense.(ST).Hrms = results_vegetation(:,2)./sqrt(2);
    RBR_Model.vegetationDense.(ST).Tp = results_vegetation(:,3);
    RBR_Model.vegetationDense.(ST).MeanDir = results_vegetation(:,4);
    RBR_Model.vegetationDense.(ST).PeakDir = results_vegetation(:,5);
    clear results_vegetation
    
        results_vegetation = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\aqdawac8forcing\rbr',...
        num2str(i) '.dat']);
    RBR_Model.aqdawacforcing.(ST).time_wave = datenum(num2str(results_vegetation(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.aqdawacforcing.(ST).Hsig = results_vegetation(:,2);
    RBR_Model.aqdawacforcing.(ST).Hrms = results_vegetation(:,2)./sqrt(2);
    RBR_Model.aqdawacforcing.(ST).Tp = results_vegetation(:,3);
    RBR_Model.aqdawacforcing.(ST).MeanDir = results_vegetation(:,4);
    RBR_Model.aqdawacforcing.(ST).PeakDir = results_vegetation(:,5);
    clear results_vegetation
    
        results_vegetation = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\aqdforcing\rbr',...
        num2str(i) '.dat']);
    RBR_Model.aqdforcing.(ST).time_wave = datenum(num2str(results_vegetation(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.aqdforcing.(ST).Hsig = results_vegetation(:,2);
    RBR_Model.aqdforcing.(ST).Hrms = results_vegetation(:,2)./sqrt(2);
    RBR_Model.aqdforcing.(ST).Tp = results_vegetation(:,3);
    RBR_Model.aqdforcing.(ST).MeanDir = results_vegetation(:,4);
    RBR_Model.aqdforcing.(ST).PeakDir = results_vegetation(:,5);
    clear results_vegetation
%     results_vegetation = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVegInterm\rbr',...
%         num2str(i) '.dat']);
%     RBR_Model.vegetationIntermFK.(ST).time_wave = datenum(num2str(results_vegetation(:,1),'%f'),'yyyymmdd.HHMMSS');
%     RBR_Model.vegetationIntermFK.(ST).Hsig = results_vegetation(:,2);
%     RBR_Model.vegetationIntermFK.(ST).Hrms = results_vegetation(:,2)./sqrt(2);
%     RBR_Model.vegetationIntermFK.(ST).Tp = results_vegetation(:,3);
%     RBR_Model.vegetationIntermFK.(ST).MeanDir = results_vegetation(:,4);
%     RBR_Model.vegetationIntermFK.(ST).PeakDir = results_vegetation(:,5);
%     clear results_vegetation
    
    results_vegetation = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegSparse3\rbr',...
        num2str(i) '.dat']);
    RBR_Model.vegetationSparse3.(ST).time_wave = datenum(num2str(results_vegetation(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.vegetationSparse3.(ST).Hsig = results_vegetation(:,2);
    RBR_Model.vegetationSparse3.(ST).Hrms = results_vegetation(:,2)./sqrt(2);
    RBR_Model.vegetationSparse3.(ST).Tp = results_vegetation(:,3);
    RBR_Model.vegetationSparse3.(ST).MeanDir = results_vegetation(:,4);
    RBR_Model.vegetationSparse3.(ST).PeakDir = results_vegetation(:,5);
    clear results_vegetation
    
%     results_vegetation = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegSparse2\rbr',...
%         num2str(i) '.dat']);
%     RBR_Model.vegetationSparse2.(ST).time_wave = datenum(num2str(results_vegetation(:,1),'%f'),'yyyymmdd.HHMMSS');
%     RBR_Model.vegetationSparse2.(ST).Hsig = results_vegetation(:,2);
%     RBR_Model.vegetationSparse2.(ST).Hrms = results_vegetation(:,2)./sqrt(2);
%     RBR_Model.vegetationSparse2.(ST).Tp = results_vegetation(:,3);
%     RBR_Model.vegetationSparse2.(ST).MeanDir = results_vegetation(:,4);
%     RBR_Model.vegetationSparse2.(ST).PeakDir = results_vegetation(:,5);
%     clear results_vegetation
    
    results_vegetation = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegSparse1\rbr',...
        num2str(i) '.dat']);
    RBR_Model.vegetationSparse1.(ST).time_wave = datenum(num2str(results_vegetation(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.vegetationSparse1.(ST).Hsig = results_vegetation(:,2);
    RBR_Model.vegetationSparse1.(ST).Hrms = results_vegetation(:,2)./sqrt(2);
    RBR_Model.vegetationSparse1.(ST).Tp = results_vegetation(:,3);
    RBR_Model.vegetationSparse1.(ST).MeanDir = results_vegetation(:,4);
    RBR_Model.vegetationSparse1.(ST).PeakDir = results_vegetation(:,5);
    clear results_vegetation
    
     results_kn = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVegSparseKn\rbr',...
        num2str(i) '.dat']);
    RBR_Model.vegetationSparsekn.(ST).time_wave = datenum(num2str(results_kn(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.vegetationSparsekn.(ST).Hsig = results_kn(:,2);
    RBR_Model.vegetationSparsekn.(ST).Hrms = results_kn(:,2)./sqrt(2);
    RBR_Model.vegetationSparsekn.(ST).Tp = results_kn(:,3);
    RBR_Model.vegetationSparsekn.(ST).MeanDir = results_kn(:,4);
    RBR_Model.vegetationSparsekn.(ST).PeakDir = results_kn(:,5);
    clear results_vegetation
       
     results_kn = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVegcd3\rbr',...
        num2str(i) '.dat']);
    RBR_Model.Veg_fake_cd03.(ST).time_wave = datenum(num2str(results_kn(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.Veg_fake_cd03.(ST).Hsig = results_kn(:,2);
    RBR_Model.Veg_fake_cd03.(ST).Hrms = results_kn(:,2)./sqrt(2);
    RBR_Model.Veg_fake_cd03.(ST).Tp = results_kn(:,3);
    RBR_Model.Veg_fake_cd03.(ST).MeanDir = results_kn(:,4);
    RBR_Model.Veg_fake_cd03.(ST).PeakDir = results_kn(:,5);
    clear results_vegetation

         results_kn = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVegcd5\rbr',...
        num2str(i) '.dat']);
    RBR_Model.Veg_fake_cd05.(ST).time_wave = datenum(num2str(results_kn(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.Veg_fake_cd05.(ST).Hsig = results_kn(:,2);
    RBR_Model.Veg_fake_cd05.(ST).Hrms = results_kn(:,2)./sqrt(2);
    RBR_Model.Veg_fake_cd05.(ST).Tp = results_kn(:,3);
    RBR_Model.Veg_fake_cd05.(ST).MeanDir = results_kn(:,4);
    RBR_Model.Veg_fake_cd05.(ST).PeakDir = results_kn(:,5);
    clear results_vegetation
    
     results_kn = load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVegcd8\rbr',...
        num2str(i) '.dat']);
    RBR_Model.Veg_fake_cd08.(ST).time_wave = datenum(num2str(results_kn(:,1),'%f'),'yyyymmdd.HHMMSS');
    RBR_Model.Veg_fake_cd08.(ST).Hsig = results_kn(:,2);
    RBR_Model.Veg_fake_cd08.(ST).Hrms = results_kn(:,2)./sqrt(2);
    RBR_Model.Veg_fake_cd08.(ST).Tp = results_kn(:,3);
    RBR_Model.Veg_fake_cd08.(ST).MeanDir = results_kn(:,4);
    RBR_Model.Veg_fake_cd08.(ST).PeakDir = results_kn(:,5);
    clear results_vegetation
 
end

%azymuth angle
azymuth = 180-acosd((-6456977+6458205)./sqrt((363817-364383)^2+(6458205-6456977)^2));
angleGrid = 180-azymuth;
%%
save ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_RBRs_model.mat','RBR_Model',...
    '-v7.3')

%% calculate \DeltaR (modeled)
rho=1023;
g=9.81;

for k=1:10 %for each instrument 
    ST = sprintf('ST%d',k);
    RBR_Model.nofric.(ST).depth_interp = interp1(RBR.(ST).time_wave,RBR.(ST).meandepth,RBR_Model.nofric.(ST).time_wave);
    RBR_Model.nofric.(ST).Hrms_interp = interp1(RBR_Model.nofric.(ST).time_wave,RBR_Model.nofric.(ST).Hrms,RBR.(ST).time_wave);
    
    for i=1:length(RBR_Model.nofric.(ST).depth_interp) %for each time
        if i>length(RBR.(ST).meandepth)
            [RBR_Model.nofric.(ST).Lp(i),RBR_Model.nofric.(ST).kp(i),RBR_Model.nofric.(ST).sigmap(i)]=disper(RBR_Model.nofric.(ST).depth_interp(i),...
                RBR_Model.nofric.(ST).Tp(i));
        else
            [RBR.(ST).Lr_p(i),RBR.(ST).kjr_p(i),RBR.(ST).sigmap(i)]=disper( ...
                RBR.(ST).meandepth(i), RBR.(ST).Tp(i));
            if isnan(RBR_Model.nofric.(ST).depth_interp(i))==0
                [RBR_Model.nofric.(ST).Lp(i),RBR_Model.nofric.(ST).kp(i),RBR_Model.nofric.(ST).sigmap(i)]=disper(RBR_Model.nofric.(ST).depth_interp(i),...
                    RBR_Model.nofric.(ST).Tp(i));
            else
                RBR_Model.nofric.(ST).Lp(i)=nan;
                RBR_Model.nofric.(ST).kp(i)=nan;
                RBR_Model.nofric.(ST).sigmap(i)=nan;
            end
        end
    end
    
    RBR_Model.nofric.(ST).E = (1/8).*rho.*g.*(RBR_Model.nofric.(ST).Hrms.^2);
    par = (2.*RBR_Model.nofric.(ST).kp'.*RBR_Model.nofric.(ST).depth_interp)./...
        (sinh(2.*RBR_Model.nofric.(ST).kp'.*RBR_Model.nofric.(ST).depth_interp));
    RBR_Model.nofric.(ST).Cg = 0.5.*(1+par).*(RBR_Model.nofric.(ST).sigmap./RBR_Model.nofric.(ST).kp)';
    
    RBR.(ST).Er = (1/8).*(rho.*g.*(RBR.(ST).Hrms.^2));
    parc = (2.*RBR.(ST).kjr_p.*RBR.(ST).meandepth)./(sinh(2.* ...
        RBR.(ST).kjr_p.*RBR.(ST).meandepth));
    RBR.(ST).Cr = 0.5.* (1+parc) .*(RBR.(ST).sigmap./(RBR.(ST).kjr_p));
    RBR.(ST).Fjg_p = RBR.(ST).Er.*RBR.(ST).Cr;
    
    RBR_Model.nofric.(ST).F = RBR_Model.nofric.(ST).E.*RBR_Model.nofric.(ST).Cg;
    RBR_Model.nofric.(ST).F_interp = interp1(RBR_Model.nofric.(ST).time_wave,RBR_Model.nofric.(ST).F,RBR.(ST).time_wave);

    RBR.(ST).Fr = RBR.(ST).Er.*RBR.(ST).Cr;
end

for k=1:9 %for each instrument
    ST = sprintf('ST%d',k);
    STf = sprintf('ST%d',k+1);
    
    RBR_Model.nofric.(ST).DeltaF = (RBR_Model.nofric.(STf).F_interp - RBR_Model.nofric.(ST).F_interp)./RBR.(ST).distance;
    RBR.(ST).DeltaF = (RBR.(STf).Fr - RBR.(ST).Fr)./RBR.(ST).distance;
    
    RBR_Model.nofric.(ST).DeltaF_norm = (RBR_Model.nofric.(STf).F_interp - RBR_Model.nofric.(ST).F_interp)./RBR_Model.nofric.(STf).Hrms_interp;
%     RBR.(ST).DeltaF_norm = (RBR.(STf).Fr - RBR.(ST).Fr)./RBR.(STf).Hrms;

%     RBR_Model.nofric.(ST).refractionIndex = RBR_Model.nofric.(ST).DeltaF./RBR.(ST).DeltaF;
%     RBR_Model.nofric.(ST).refractionIndex_norm = RBR_Model.nofric.(ST).DeltaF_norm./RBR.(ST).DeltaF_norm;

%     RBR_Model.nofric.(ST).refractionIndex(RBR_Model.nofric.(ST).refractionIndex>1 | RBR_Model.nofric.(ST).refractionIndex<0)=nan;
%     sum(isnan(RBR_Model.nofric.(ST).refractionIndex))
%     figure(k)
%     plot(RBR_Model.nofric.(ST).refractionIndex)
%         figure(k+10)

%     fitdist(RBR_Model.nofric.(ST).refractionIndex','Kernel','Kernel','epanechnikov')
end

%%
save ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_RBRs_model_with_CD.mat','RBR_Model',...
    '-v7.3')

% %% jonswap
% g=9.806; Beta=5/4;
% index=650;
% Hs=RBR.ST4.Hm0(:,index);
% Gamma=1;
% Tz=RBR.ST4.Tp(:,index);
% Tm= Tz*(0.327*exp(-0.315*Gamma)+1.17);
% alphabar = 5.058*(1-.287*log(Gamma))*(Hs/Tm^2)^2  ;
% Omega=2.*pi.*RBR.ST4.f(:,1);
% Omegam    = 2*pi/Tm;
% SigmaA=0.07;  %spectral width parameter
% SigmaB=0.09;  %spectral width parameter
% sigma = (Omega<=Omegam)*SigmaA+(Omega>Omegam)*SigmaB;
% A     = exp(-((Omega/Omegam-1)./(sigma*sqrt(2))).^2);  
% 
% S     = alphabar*g^2 .* Omega.^-5 .* exp(-(Beta*(Omega/Omegam).^-4)) .* Gamma.^A;      %spectra m^2.s
% S     = 0.0081*g^2 .* Omega.^-5 .* exp(-(0.74*((g/15)./Omega).^-4)) .* Gamma.^A;      %spectra m^2.s

