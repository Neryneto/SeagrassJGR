function [vel_maj, vel_min, Current_Mean,...
    Current_theta, Wave_theta,...
    Despike_curr_maj, Despike_curr_min, fx, ubr,ubr_spect,...
    Suu,Svv] = rotate_continuous_ADV_max_variance(VELX, VELY, velx_despike,vely_despike,frequency)

%% ANALYSIS_REYNOLDS_STRESS_ADP
%  This function rotates the ADV data into the maximal direction, removes
%  bad data and then calculates the mean and RMS velocity of the signal.

% Determine the size of data (all ADVs are the same)
VELX=VELX';VELY=VELY';velx_despike=velx_despike';vely_despike=vely_despike';
time_length = 2.^nextpow2(frequency*60*15);
[BURSTS] = floor(size(VELX,1)/time_length);

for bi = 1:BURSTS(1)
    
    
    % Post an update in the command window
    if rem(bi,100) == 0;
        clc; display(sprintf('Burst %g of %g', bi, BURSTS));
    end
    
    % Extract data burst
    xtmp     = VELX(bi*time_length-(time_length-1):bi*time_length,1);
    ytmp     = VELY(bi*time_length-(time_length-1):bi*time_length,1);
    xtmp2    = velx_despike(bi*time_length-(time_length-1):bi*time_length,1);
    ytmp2    = vely_despike(bi*time_length-(time_length-1):bi*time_length,1);
    
    Current_Mean (bi,:) = sqrt((nanmean(xtmp2).^2) + (nanmean(ytmp2).^2));
    
    % Only continue if there is enough data
    
    if sum (~isnan(ytmp2))>length(ytmp2)/2
        %  ROTATE THE DATA INTO MAXIMAL DIRECTION FOR WAVES AND CURRENTS
        %  Determine the direction of the waves in the boundary layer
        
        w = complex(xtmp2,ytmp2);
        w (w~=isfinite(w)==0);
        
        % Covariance matrix
        cv = cov(real(w),imag(w),'omitrows');
        
        % Theta is direction of maximum variance
        Wave_theta(bi) = wrapTo360(180./pi)*(0.5*atan2(2*cv(1,2), (cv(1,1) - cv(2,2))));
        
        wr = w.*exp(-1i.*Wave_theta(bi)*pi/180);
        vel_maj (bi*time_length-(time_length-1):bi*time_length,1) = real(wr);
        vel_min (bi*time_length-(time_length-1):bi*time_length,1) = imag(wr);
        
        % Determine the current angle
        Current_theta (bi) = wrapTo360(atan2d(nanmean(ytmp2),nanmean(xtmp2)));
        Despike_curr_maj (bi*time_length-(time_length-1):bi*time_length,1) = cosd(Current_theta(bi)).*xtmp2 - ...
            sind(Current_theta(bi)).*ytmp2;
        Despike_curr_min (bi*time_length-(time_length-1):bi*time_length,1) = sind(Current_theta(bi)).*xtmp2 + ...
            cosd(Current_theta(bi)).*ytmp2;
        
        % Length of signal
        aux = nandetrend(vel_maj(bi*time_length-(time_length-1):bi*time_length,1)); 
        if isnan(vel_maj(bi*time_length-(time_length-1)))==1; aux(1)=0; end 
        if isnan(vel_maj(bi*time_length))==1; aux(end)=0; end;
        X = ~isnan(aux);
        Y = cumsum(X-diff([1; X])/2); eta = interp1(1:nnz(X),aux(X),Y);
        [Suu(bi,:),fx]=calculatespectrum(frequency,eta,512); clear aux X Y eta
        
        aux = nandetrend(vel_min(bi*time_length-(time_length-1):bi*time_length,1)); 
        if isnan(vel_min(bi*time_length-(time_length-1)))==1; aux(1)=0; end; 
        if isnan(vel_min((bi*time_length)))==1; aux(end)=0; end;
        X = ~isnan(aux);
        Y = cumsum(X-diff([1; X])/2); eta = interp1(1:nnz(X),aux(X),Y);
        [Svv(bi,:),fy]=calculatespectrum(frequency,eta,512); clear L aux X Y eta
        
        %calculate representative wave orbital velocities
        ubr_spect (bi) = sqrt(2*(nansum(Suu(bi,:)+Svv(bi,:))*(fx(2)-fx(1)))); 
        
        %calculate spectral wave orbital velocities (high frequency)
        ubr (bi) = sqrt(2*(nanvar(vel_maj (bi*time_length-(time_length-1):bi*time_length,1)) + ...
            nanvar(vel_min (bi*time_length-(time_length-1):bi*time_length,1))));
        
        clear xtmp xtmp2 ytmp ytmp2
        
    else
        Wave_Urms(bi) = NaN;
        Current_Mean(bi) = NaN;
        Wave_theta(bi) = NaN;       
        Current_theta(bi) = NaN;
        
        vel_maj(bi*time_length-(time_length-1):bi*time_length,:) = NaN;
        vel_min(bi*time_length-(time_length-1):bi*time_length,:) = NaN;
        Despike_curr_maj (bi*time_length-(time_length-1):bi*time_length,:) = NaN;
        Despike_curr_min (bi*time_length-(time_length-1):bi*time_length,:) = NaN;
        Suu(bi,1:257) = NaN;Svv(bi,1:257) = NaN;
        ubr_spect (bi) = nan;
        ubr (bi) = nan;
        clear xtmp xtmp2 ytmp ytmp2
    end
    
end


end



