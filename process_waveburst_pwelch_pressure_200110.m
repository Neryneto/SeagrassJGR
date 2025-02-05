function [Hsig_sea,Hsig_swell,Hrms_swell,Hsig_IG,Hrms_IG,Hsig_fIG,Hrms_fIG,Tp_swell,Ta_swell,S,f_welch,depth]=process_waveburst_pwelch_surface_200110(pressure_raw,dt,chunk,h,swell_freqs,IG_freqs,fIG_freqs)

% Spectral analysis for pressure sensor data

g=9.81;
fs=1/dt;  % Sample frequency in Hz

% Determine fundamental frequency
f0=1/(dt*chunk);

    z=nanmean(pressure_raw,1);
    depth=z+h;   % h=height of pressure sensor above the bottom;
    pressure=detrend(pressure_raw);  % Detrend data

    [Spp,f_welch]=pwelch(pressure,hann(chunk),[],chunk,'onesided',fs);

        for i=1:length(f_welch)
            k_depthav(i)=lindisp_explicit(f_welch(i),depth);
            K_kin_depthav(i)=cosh(k_depthav(i)*(depth-z))/cosh(k_depthav(i)*depth);
        end
%    kmaxL=pi./(z-h);
% KpminL=cosh(kmaxL*0.05)/cosh(kmaxL*h); % Minimum Limit for K_p calculated based on linear wave theory
% K_kin_depthav(K_kin_depthav < KpminL) = KpminL; % Check to avoid large amplification, Kp should be larger than minimum K_p calculated based on linear wave theory

        S=(Spp)./(K_kin_depthav'.^2);

        
        % Swell ************************************************************
        % Find point nearest low freq
        [junk,ind_low]=min(abs(f_welch-swell_freqs(1)));
        % Find point nearest high freq
        [junk,ind_high]=min(abs(f_welch-swell_freqs(2)));
        Var_tot1=trapz(f_welch(ind_high:end),S(ind_high:end));
        Var_tot=trapz(f_welch(ind_low:ind_high),S(ind_low:ind_high));
        mom1_tot=trapz(f_welch(ind_low:ind_high),S(ind_low:ind_high).*f_welch(ind_low:ind_high));
        Hrms_swell=sqrt(8*Var_tot);
        Hsig_swell=4*sqrt(Var_tot);
        Hsig_sea = 4*sqrt(Var_tot1);
        % Peak period
        [dum,max_I]=max(S);
        Tp_swell=1/f_welch(max_I);
        % Mean period
        Ta_swell=trapz(f_welch(ind_low:ind_high),S(ind_low:ind_high))/trapz(f_welch(ind_low:ind_high),f_welch(ind_low:ind_high).*S(ind_low:ind_high));     
        % ***************************************************************

        % IG ************************************************************
        % Find point nearest low freq
        [junk,ind_low]=min(abs(f_welch-IG_freqs(1)));
        % Find point nearest high freq
        [junk,ind_high]=min(abs(f_welch-IG_freqs(2)));
        Var_tot=trapz(f_welch(ind_low:ind_high),S(ind_low:ind_high));
        Hrms_IG=sqrt(8*Var_tot);
        Hsig_IG=4*sqrt(Var_tot);
        
         % fIG ************************************************************
        % Find point nearest low freq
        [junk,ind_low]=min(abs(f_welch-fIG_freqs(1)));
        % Find point nearest high freq
        [junk,ind_high]=min(abs(f_welch-fIG_freqs(2)));
        Var_tot=trapz(f_welch(ind_low:ind_high),S(ind_low:ind_high));
        Hrms_fIG=sqrt(8*Var_tot);
        Hsig_fIG=4*sqrt(Var_tot);
        
        S=S';  % Correct output format

end
   