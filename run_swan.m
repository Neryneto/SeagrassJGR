
% for CD=[5:8]
for CD=1:8

    for h=1:4
    cd(fullfile(['C:\Users\22371812\Desktop\cdh\cd' num2str(CD) 'h' num2str(h)]))

    a=strcat({'swanrun aqdforcing'});

    [s w]=dos(cell2mat(a));
    end
    
end

for CD=1:8
    for h=1:4
        for rbr=1:8
            
            cdhrbr = sprintf('Cd_%d_h_%d_rbr_%d',CD,h,rbr);
%             results=load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\cdh\cd' num2str(CD) 'h' num2str(h) '\rbr' num2str(rbr) '.dat']);
            results=load (['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\cdh2\cd' num2str(CD) 'h' num2str(h) '\rbr' num2str(rbr) '_spc.dat']);

            RBR_Modeled.(cdhrbr).time_wave = datenum(num2str(results(:,1),'%f'),'yyyymmdd.HHMMSS');
            RBR_Modeled.(cdhrbr).Hsig = results(:,2);
            RBR_Modeled.(cdhrbr).Hrms = results(:,2)./sqrt(2);
            RBR_Modeled.(cdhrbr).Tp = results(:,3);
            RBR_Modeled.(cdhrbr).MeanDir = results(:,4);
            RBR_Modeled.(cdhrbr).PeakDir = results(:,5);  
            clear results
            
        end
    end
    
end

save ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\cdh2\RBR_modeled','RBR_Modeled','-v7.3')

%%
load ('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_RBRs_model.mat')
%%

for CD=1:8
    for h=1:4        
        cdh = sprintf('Cd_%d_h_%d_rbr_%d',CD,h,2);
        RBR_Modeled.(cdh).Hsig_interp = interp1(RBR_Modeled.(cdh).time_wave,RBR_Modeled.(cdh).Hsig,RBR.ST2.time_wave);
        RBR_Modeled.(cdh).Hsig_interp(RBR_Modeled.(cdh).Hsig_interp<0.131)=nan;
        Wilmott (h,CD) = 1 - ((nansum((RBR.ST2.Hm0 - RBR_Modeled.(cdh).Hsig_interp).^2))./...
            (nansum((abs(RBR_Modeled.(cdh).Hsig_interp-nanmean(RBR.ST2.Hm0))+abs(RBR.ST2.Hm0-nanmean(RBR.ST2.Hm0))).^2)));
        Murphy(h,CD) = 1 - ((nansum((RBR_Modeled.(cdh).Hsig_interp-RBR.ST2.Hm0).^2))./...
            (nansum((RBR.ST2.Hm0-nanmean(RBR.ST2.Hm0)).^2)));
    end
end

