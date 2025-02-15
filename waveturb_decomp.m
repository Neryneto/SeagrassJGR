function [upwp, TKEturb] = waveturb_decomp(ib,Tp,dnewtl,dpnewtl)
%%=========================================================================
% This function provides code to separate wave and turbulence components of 
% Reynold Stresses and EC oxygen fluxes by the phase method 
% [Bricker and Monismith, 2007]

% This function was developed by C.E. Reimers for the manuscript published 
% in JGR-Oceans by Reimers and Fogaren (2021). 

%*** INPUTS:
% ib = burst number (-), not used in this function, may be used for data 
%   output purposes
% Tp = bottom sensed wave period (s) 
% Data columns for dnewtl and dpnewtl
%   dnewtl = data after rotation and time lag correction for oxygen 
%   dpnewtl = detrended data series 
%        col 1: time (s) 
%        col 2: horizontal velocity in dominate direction of waves (m/s)
%        col 3: horizontal velocity perpendicular to u (m/s)
%        col 4: vertical velocity after rotation to minimize wave           
%            components of w (m/s)
%        col 5: oxygen concentration (umol/L)

% *** OUTPUTS:
% Oxfluxes: (mmol m-2 d-1)
%        col 1: Oxygen Flux - Total Flux
%        col 2: Oxygen Flux - Bandpass wave component 
%        col 3: Oxygen Flux - Phase wave component 
%        col 4: Oxygen Flux - Turbulent Flux (Total - Phase wave component) 
% TKEturb: turbulent kinetic energy (cm2 s-2); calculated from variances
%   which removed wave components using the phase method 

%%=========================================================================

hz = 4; % sampling frequency (Hz) 
dt = 1/hz; % sampling interval (s)
t = dnewtl(:,1);  % time (s)

% Time series of velocity and oxygen measurements
U = dnewtl(:,2); % horizontal velocity in dominate direction of waves, cm/s units
V = dnewtl(:,3); % horizontal velocity perpendicular to u
W = dnewtl(:,4); % vertical velocity after rotation to minimize wave components of w
P = dnewtl(:,5); % oxygen concentration

% Time series of detrended velocity and oxygen measurements
up = U; 
vp = V;
wp = W;
cp = P;

if isnan(P(1))==1; P(1)=0; end
if isnan(P(end))==1; P(end)=0; end
nanp = isnan(P); tt = 1:numel(P);
cp(nanp) = interp1(tt(~nanp), P(~nanp), tt(nanp));

if isnan(U(1))==1; U(1)=0; end
if isnan(U(end))==1; U(end)=0; end
nanx = isnan(U); tt = 1:numel(U);
up(nanx) = interp1(tt(~nanx), U(~nanx), tt(nanx));

if isnan(V(1))==1; V(1)=0; end
if isnan(V(end))==1; V(end)=0; end
nany = isnan(V); tt = 1:numel(V);
vp(nany) = interp1(tt(~nany), V(~nany), tt(nany));

if isnan(W(1))==1; W(1)=0; end
if isnan(W(end))==1; W(end)=0; end
nanw = isnan(W); tt = 1:numel(W);
wp(nanw) = interp1(tt(~nanw), W(~nanw), tt(nanw));

df = 1./(dt*(length(wp)-1));
nt = size(wp,1);

f = (0:df:df*(length(wp)-1)); % frequency vector
nny = ceil((length(wp)./2));
%Fourier transform of cp 
Amc = fft(cp)./(nt-1); 
% Fourier transforms of up, vp, and wp. 
Amu = fft(up-nanmean(up))./(nt-1); 
Amv = fft(vp)./(nt-1);
Amw = fft(wp)./(nt-1);

Co = real(Amc.*conj(Amw));
% Total EC 02 flux, converted to mmol m-2 d-1
wpcp_alt = sum(Co)*86400;  

% Compute power densities
Sww = 2.*abs(Amw(1:nny)).^2/df;
Scc = 2.*abs(Amc(1:nny)).^2/df;
Suu = 2.*abs(Amu(1:nny)).^2/df;
Svv = 2.*abs(Amv(1:nny)).^2/df;
fm = f(1:nny);

% Compute phase, output in radians 
Uph = atan2(imag(Amu),real(Amu));  
Vph = atan2(imag(Amv),real(Amv));
Wph = atan2(imag(Amw),real(Amw));
Cph = atan2(imag(Amc),real(Amc));

% Computing the cross spectra
Suv = 2.*real(Amu(1:nny).*conj(Amv(1:nny)))/(df);
Suw = 2.*real(Amu(1:nny).*conj(Amw(1:nny)))/(df);
Svw = 2.*real(Amv(1:nny).*conj(Amw(1:nny)))/(df);
Swc = 2.*real(Amw(1:nny).*conj(Amc(1:nny)))/(df);
Suc = 2.*real(Amu(1:nny).*conj(Amc(1:nny)))/(df);
Svc= 2.*real(Amv(1:nny).*conj(Amc(1:nny)))/(df);

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

%check on peak assignments
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
% pause

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
% 
% figure
% loglog(fm,Svv,'k-',fWave,Svv_wave,'r-',fWave,y,'b-')
% title('Svv')
[g,h] = size(Svv_wave);
Svv_wavecomp = zeros(g,h);
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
Scc_inter = Scc(interprange);
Scc_wave = Scc(waverange);

%Linear interpolation over turbulent spectra
Sl = log10(Scc_inter);
[P(1:2),s] = polyfit(F,Sl,1);
    slope = P(1);
    interc = P(2);
logy = (slope*(log10(fWave)) + interc);   
y = 10.^(logy);

% figure
% loglog(fm,Scc,'k-',fWave,Scc_wave,'r-',fWave,y,'b-')
% title('Scc')
% ylabel('S_c_c (({\mu}mol L^-^1)^2 Hz^-^1)');
% xlabel('f (Hz)');
% hb = gca;
% set(hb,'position',[.15 .30 .7 .40]);

[g,h] = size(Scc_wave);
Scc_wavecomp = zeros(g,h);
y = y';
for b = 1:g
    Scc_wavecomp(b) = (Scc_wave(b) - y(b)); 
end
% wave fourier component
Amcc_wave = sqrt((Scc_wavecomp + 0j).*(df));
% pause   %use to reset ranges if needed

%Wave Magnitudes
Um_wave = sqrt(real(Amuu_wave).^2 + imag(Amuu_wave).^2);
Vm_wave = sqrt(real(Amvv_wave).^2 + imag(Amvv_wave).^2);
Wm_wave = sqrt(real(Amww_wave).^2 + imag(Amww_wave).^2);
Cm_wave = sqrt(real(Amcc_wave).^2 + imag(Amcc_wave).^2);
   
%Wave reynolds stresses
uw_wave = nansum(Um_wave.*Wm_wave.*cos(Wph(waverange) - Uph(waverange)));
uv_wave =  nansum(Um_wave.*Vm_wave.*cos(Vph(waverange) - Uph(waverange)));
vw_wave = nansum(Vm_wave.*Wm_wave.*cos(Wph(waverange) - Vph(waverange)));
wc_wave=nansum(Cm_wave.*Wm_wave.*cos(Wph(waverange) - Cph(waverange)));

%sum wave components
uu_wave = nansum(Suu_wavecomp*df); 
vv_wave = nansum(Svv_wavecomp*df);
ww_wave = nansum(Sww_wavecomp*df);
cc_wave = nansum(Scc_wavecomp*df);

%compare to cw_wave
wc_wavebandpass = nansum(Swc(waverange)*df); 
 
%Full Reynolds stresses
uu = nansum(real(Suu)*df);
uu_alt = nansum(real(Amu.*conj(Amu))); %should match
uv = nansum(real(Suv)*df);
uv_alt = nansum(real(Amu.*conj(Amv))); % Reason this isn't nansum? 
uw = nansum(real(Suw)*df);
vv = nansum(real(Svv)*df);
vw = nansum(real(Svw)*df);
ww = nansum(real(Sww)*df);
wc= nansum(real(Swc)*df);
cc= nansum(real(Scc)*df);

%Turbulent reynolds stresses corrected
upup = uu - uu_wave;
vpvp = vv - vv_wave;
wpwp = ww - ww_wave;
cpcp = cc - cc_wave;
upwp = uw - uw_wave;
upvp = uv - uv_wave;
vpwp = vw - vw_wave;
wpcp = wc - wc_wave; 

fullfluxes = [uu vv ww cc uw uv vw wc];
tfluxes = [upup vpvp wpwp cpcp upwp upvp vpwp wpcp];

TKEturb = 0.5*(upup + vpvp + wpwp); % in units m2 s-2 

wpcp_turbperc = wpcp/wc*100;
phaseeffect = wc_wave/wc_wavebandpass;
Oxfluxes = [wc*86400 wc_wavebandpass*86400 wc_wave*86400 wpcp*86400];

% figure
% bar(Oxfluxes);
% xticks('manual');
% xticklabels({'wc', 'wc-wavebandpass','wc-wave', 'wc-turb'});
% ylabel('Oxygen Flux (mmol m^-^2 d^-^1)');
% pause
end