function [upup, vpvp, wpwp, upwp, uw_wave, TKEturb] = calculate_BrickerMonismith_aqp (U,V,W,P,fs,Tp)

%Function adpted from Python

N = length (U);

%% Turbulent reynolds stresses
uw = zeros (N,1);
vw = zeros (N,1);
uv = zeros (N,1);
uu = zeros (N,1);
vv = zeros (N,1);
ww = zeros (N,1);

% Wave reynolds stresses
uw_wave = zeros(N,1);
vw_wave = zeros(N,1);
uv_wave = zeros(N,1);
uu_wave = zeros(N,1);
vv_wave = zeros(N,1);
ww_wave = zeros(N,1);

% Interpolating out nans
if isnan(P(1))==1; P(1)=0; end
if isnan(P(end))==1; P(end)=0; end
nanp = isnan(P); t = 1:numel(P);
P(nanp) = interp1(t(~nanp), P(~nanp), t(nanp));

if isnan(U(1))==1; U(1)=0; end
if isnan(U(end))==1; U(end)=0; end
nanx = isnan(U); t = 1:numel(U);
U(nanx) = interp1(t(~nanx), U(~nanx), t(nanx));

if isnan(V(1))==1; V(1)=0; end
if isnan(V(end))==1; V(end)=0; end
nany = isnan(V); t = 1:numel(V);
V(nany) = interp1(t(~nany), V(~nany), t(nany));

if isnan(W(1))==1; W(1)=0; end
if isnan(W(end))==1; W(end)=0; end
nanw = isnan(W); t = 1:numel(W);
W(nanw) = interp1(t(~nanw), W(~nanw), t(nanw));

% dbar = 1e4 * (nanmean(P) + doffp )/(rho*g) %Average water depth;

%Making sure average depth is positive
% try dbar<0
%     nfft = length(U);
% catch
%     disp('error');
% end
%
dt = 1/fs; % sampling interval (s)
nfft = length(U)-1;
jj = 0:nfft;
df = 1./(dt*(length(W)-1));
WIN = hanning(length(jj)-1);
nny = ceil((length(W)./2));
f = (0:df:df*(length(W)-1)); % frequency vector
fm = f(1:nny);

%% Computing full spectra
Amu = fft(U-nanmean(U))./nfft;
Amv = fft(V-nanmean(V))./nfft;
Amw = fft(W-nanmean(W))./nfft;
Amp = fft(P-nanmean(P))./nfft;
 
%% Compute power densities
Suu=2*abs(Amu(1:nny)).^2./df ;

Svv=2*abs(Amv(1:nny)).^2./df; 

Sww=2*abs(Amw(1:nny)).^2./df;

Spp=2*abs(Amp(1:nny)).^2./df;

%% Compute phase, output in radians 
Uph = atan2(imag(Amu),real(Amu));  
Vph = atan2(imag(Amv),real(Amv));
Wph = atan2(imag(Amw),real(Amw));

%% Computing the cross spectra
Suv = 2.*real(Amu(1:nny).*conj(Amv(1:nny)))/(df);
Suw = 2.*real(Amu(1:nny).*conj(Amw(1:nny)))/(df);
Svw = 2.*real(Amv(1:nny).*conj(Amw(1:nny)))/(df);
Spu = 2.*real(Amp(1:nny).*conj(Amu(1:nny)))/(df);
Spv = 2.*real(Amp(1:nny).*conj(Amv(1:nny)))/(df);
Spw= 2.*real(Amp(1:nny).*conj(Amw(1:nny)))/(df);

%% Filtering
% Search for the wave peak within a reasonable range of frequencies -- 
% Range of frequencies needs to be adjusted for each data set
fmax = 1/Tp;  % maximum frequency of wave peak
% Lower end of the wave range as a fraction of the peak frequency, 
%   must adjust for each burst according to peak width
wr_low = 0.3; 
%Upper end of the wave range as a fraction of the peak frequency, 
%   adjust for each burst
wr_high = 1.2; 
waverangemin = fmax-(fmax*wr_low);
waverangemax = fmax+(fmax*wr_high);
waverange = find(fm > waverangemin & fm < waverangemax); %indices
fWave = fm(waverange); %frequencies for waves

% Calculate range of freq for linear interpretation (set limits):
%   must adjust for each data set

%range for turb slope, adjust depending on where wave peak sits
interprange = find(fm >(fmax +(fmax*wr_high)) & fm < 0.6);  
anchor = find(fm < (fmax -(fmax*wr_low)) & fm > (fmax-2*(fmax*wr_low)));
interprange = cat(1,anchor',interprange');
fInter = fm(interprange');
Suu_inter = Suu(interprange);
Suu_wave = Suu(waverange);
% figure
% loglog(fm,Suu,fWave,Suu_wave,'r-',fInter,Suu_inter','g-x');  

% Linear interpolation of turbulent spectra beneath wave peak
F = log10(fInter)'; 
Sl = log10(Suu_inter);
% plot(F,Sl,'x');
[P(1:2),s] = polyfit(F,Sl,1);
    slope = P(1);
    interc = P(2);
logy = (slope*(log10(fWave)) + interc);   
y = 10.^(logy);

% figure
% loglog(fm,Suu,'k-',fWave, Suu_wave,'r-',fWave,y,'b--');
% title('Suu')
% ylabel('S_u_u (m^2 s ^-^2 Hz^-^1)');
% xlabel('f (Hz)');
% hb = gca;
% set(hb,'position',[.15 .30 .7 .40]);
% hb.FontSize = 12;

[g,h] = size(Suu_wave);
Suu_wavecomp = zeros(g,h);
y = y';
for b = 1:g   
   Suu_wavecomp(b) = (Suu_wave(b) - y(b));  
end
%wave fourier component
Amuu_wave = sqrt((Suu_wavecomp+0j).*(df));


%now for vv    
Svv_inter = Svv(interprange);
Svv_wave = Svv(waverange);

%Linear interpolation over turbulent spectra
Sl = log10(Svv_inter);
[P(1:2),s] = polyfit(F,Sl,1);
    slope = P(1);
    interc = P(2);
logy = (slope*(log10(fWave)) + interc);   
y = 10.^(logy);

% figure
% loglog(fm,Svv,'k-',fWave,Svv_wave,'r-',fWave,y,'b-')
% title('Svv')
% [g,h] = size(Svv_wave);
% Svv_wavecomp = zeros(g,h);
y = y';
for b=1:g
    Svv_wavecomp(b) = (Svv_wave(b) - y(b));
end
%wave fourier component
Amvv_wave = sqrt((Svv_wavecomp+0j).*(df));      

%now for ww    
Sww_inter = Sww(interprange);
Sww_wave = Sww(waverange);

%Linear interpolation over turbulent spectra
Sl = log10(Sww_inter);
[P(1:2),s] = polyfit(F,Sl,1);
    slope = P(1);
    interc = P(2);
logy = (slope*(log10(fWave)) + interc);   
y = 10.^(logy);

% figure
% loglog(fm,Sww,'k-',fWave,Sww_wave,'r-',fWave,y,'b-')
% title('Sww')
% ylabel('S_w_w (m^2 s ^-^2 Hz^-^1)');
% xlabel('f (Hz)');
% hb = gca;
% set(hb,'position',[.15 .30 .7 .40]);

[g,h] = size(Sww_wave);
Sww_wavecomp = zeros(g,h);
y = y';
for b=1:g   
   Sww_wavecomp(b) = (Sww_wave(b) - y(b));  
end
%wave fourier component
Amww_wave = sqrt((Sww_wavecomp+0j).*(df));

%now for cc
Spp_inter = Spp(interprange);
Spp_wave = Spp(waverange);

%Linear interpolation over turbulent spectra
Sl = log10(Spp_inter);
[P(1:2),s] = polyfit(F,Sl,1);
    slope = P(1);
    interc = P(2);
logy = (slope*(log10(fWave)) + interc);   
y = 10.^(logy);

% figure
% loglog(fm,Spp,'k-',fWave,Spp_wave,'r-',fWave,y,'b-')
% title('Spp')
% ylabel('S_c_c (({\mu}mol L^-^1)^2 Hz^-^1)');
% xlabel('f (Hz)');
% hb = gca;
% set(hb,'position',[.15 .30 .7 .40]);

[g,h] = size(Spp_wave);
Spp_wavecomp = zeros(g,h);
y = y';
for b = 1:g
    Spp_wavecomp(b) = (Spp_wave(b) - y(b)); 
end
% wave fourier component
Amcc_wave = sqrt((Spp_wavecomp + 0j).*(df));

%Wave Magnitudes
Um_wave = sqrt(real(Amuu_wave).^2 + imag(Amuu_wave).^2);
Vm_wave = sqrt(real(Amvv_wave).^2 + imag(Amvv_wave).^2);
Wm_wave = sqrt(real(Amww_wave).^2 + imag(Amww_wave).^2);
Cm_wave = sqrt(real(Amcc_wave).^2 + imag(Amcc_wave).^2);
   
%Wave reynolds stresses
uw_wave = nansum(Um_wave.*Wm_wave.*cos(Wph(waverange) - Uph(waverange)));
uv_wave =  nansum(Um_wave.*Vm_wave.*cos(Vph(waverange) - Uph(waverange)));
vw_wave = nansum(Vm_wave.*Wm_wave.*cos(Wph(waverange) - Vph(waverange)));

%sum wave components
uu_wave = nansum(Suu_wavecomp*df); 
vv_wave = nansum(Svv_wavecomp*df);
ww_wave = nansum(Sww_wavecomp*df);
pp_wave = nansum(Spp_wavecomp*df);

%compare to cw_wave
wc_wavebandpass = nansum(Spw(waverange)*df); 
 
%Full Reynolds stresses
uu = nansum(real(Suu)*df);
uu_alt = nansum(real(Amu.*conj(Amu))); %should match
uv = nansum(real(Suv)*df);
uv_alt = nansum(real(Amu.*conj(Amv))); % Reason this isn't nansum? 
uw = nansum(real(Suw)*df);
vv = nansum(real(Svv)*df);
vw = nansum(real(Svw)*df);
ww = nansum(real(Sww)*df);
wc= nansum(real(Spw)*df);
cc= nansum(real(Spp)*df);

%Turbulent reynolds stresses corrected
upup = uu - uu_wave;
vpvp = vv - vv_wave;
wpwp = ww - ww_wave;

upwp = uw - uw_wave;
upvp = uv - uv_wave;
vpwp = vw - vw_wave;

% fullfluxes = [uu vv ww cc uw uv vw wc];

TKEturb = 0.5*(upup + vpvp + wpwp)*10000; % in units cm2 s-2 

% phaseeffect = wc_wave/wc_wavebandpass;

% figure
% bar(Oxfluxes);
% xticks('manual');
% xticklabels({'wc', 'wc-wavebandpass','wc-wave', 'wc-turb'});
% ylabel('Oxygen Flux (mmol m^-^2 d^-^1)');
% pause

end