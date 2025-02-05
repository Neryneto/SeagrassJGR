function [upup, vpvp, wpwp, upwp, uw_wave, tke] = calculate_BrickerMonismith (U,V,W,P,fs)

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
nfft = length(U);
jj = 0:nfft;

WIN = hanning(length(jj)-1);

%%     If using sig.welch/sig.csd
Amu = fft(U,nfft);

Amv = fft(V,nfft);

Amw = fft(W,nfft);

df = fs./(nfft);
nnyq = floor(nfft/2);
fm=(1:nnyq)*df;

%% Phase - eq 11

Uph = atan2(imag(Amu),real(Amu));
Vph = atan2(imag(Amv),real(Amv));
Wph = atan2(imag(Amw),real(Amw));

%% Computing full spectra

% [Suu,fm] = pwelch(nandetrend(U),WIN,.5,nfft);
% 
% [Svv,fm] = pwelch(nandetrend(V),WIN,.5,nfft);
% 
% [Suv,fm] = cpsd(nandetrend(U),nandetrend(V),WIN,.5,nfft);
% 
% [Sww,fm] = pwelch(nandetrend(W),WIN,.5,nfft);
% 
% [Suw,fm] = cpsd(nandetrend(U),nandetrend(W),WIN,.5,nfft);
% 
% [Svw,fm] = cpsd(nandetrend(V),nandetrend(W),WIN,.5,nfft);
% 
% [Spp,fm] = pwelch(nandetrend(P),WIN,.5,nfft);
% 
% [Spu,fm] = cpsd(nandetrend(P),nandetrend(U),WIN,.5,nfft);
% 
% [Spv,fm] = cpsd(nandetrend(P),nandetrend(V),WIN,.5,nfft);
% 
% [Spw,fm] = cpsd(nandetrend(P),nandetrend(W),WIN,.5,nfft);

% Suu = real(Amu.*conj(Amu))/(nnyq.*df); %more energetic - why??
% Suu = Suu(1:length(fm));
% 
% Svv = real(Amv.*conj(Amv))/(nnyq.*df); %more energetic - why??
% Svv = Svv(1:length(fm));
% 
% Sww = real(Amw.*conj(Amw))/(nnyq.*df); %more energetic - why??
% Sww = Sww(1:length(fm));
% 
% Suv = real(Amu.*conj(Amv))/(nnyq.*df); %more energetic - why??
% Suv = Suv(1:length(fm));
% 
% Suw = real(Amu.*conj(Amw))/(nnyq.*df); %more energetic - why??
% Suw = Suw(1:length(fm));
% 
% Svw = real(Amv.*conj(Amw))/(nnyq.*df); %more energetic - why??
% Svw = Svw(1:length(fm));
nfft = 2.^(nextpow2(length(U)));
df = fs/(nfft);
nyqFreq = fs/2;
freqRange = df:df:nyqFreq;
nnyq=floor(nfft/2)+1;
fm = (0:nnyq-1)*df;

Amu = fft(U-nanmean(U))./nfft;
Amv = fft(V-nanmean(V))./nfft;
Amw = fft(W-nanmean(W))./nfft;
Amp = fft(P-nanmean(P))./nfft;
 
Suu=[2*abs(Amu(2:nfft/2)).^2; abs(Amu(nfft/2+1)).^2] ;
Suu=Suu/df; 

Svv=[2*abs(Amv(2:nfft/2)).^2; abs(Amv(nfft/2+1)).^2] ;
Svv=Svv/df; 

Sww=[2*abs(Amw(2:nfft/2)).^2; abs(Amw(nfft/2+1)).^2] ;
Sww=Sww/df; % dividing by dn

Spp=[2*abs(Amp(2:nfft/2)).^2; abs(Amp(nfft/2+1)).^2] ;
Spp=Spp/df; % dividing by dn

SUV=cpsd(detrend(U),detrend(V),WIN,0.5,nfft-1,fs,'centered');   
Suv=abs(SUV);

SUW=cpsd(detrend(U),detrend(W),WIN,0.5,nfft-1,fs,'centered');   
Suw=abs(SUW);

SVW=cpsd(detrend(V),detrend(W),WIN,0.5,nfft-1,fs,'centered');   
Svw=abs(SVW);

SPU=cpsd(detrend(P),detrend(U),WIN,0.5,nfft-1,fs,'centered');   
Spu=abs(SPU);

SPV=cpsd(detrend(P),detrend(V),WIN,0.5,nfft-1,fs,'centered');   
Spv=abs(SPV);

SPW=cpsd(detrend(P),detrend(W),WIN,0.5,nfft-1,fs,'centered');   
Spw=abs(SPW);

fm=fm(2:end);
%% Filtering
offset = sum(fm<=0.1);

[c, umax] = max(Suu(fm>0.1 & fm < 0.7)); uumax = umax + offset;

widthratiolow = 2.333; 
widthratiohigh = 1.4; 
fmmax = fm(uumax);

waverange = uumax-floor((fmmax/widthratiolow)./df) : uumax+floor((fmmax/widthratiohigh)./df);

[lmin cmin] = nanmin(abs(fm-1));
interprange = 1:cmin;

interprangeW = 1:cmin;

interprange = interprange(interprange>=0 & interprange<nnyq);
waverange = waverange(waverange>=0 & waverange<nnyq);
interprangeW = interprangeW(interprangeW >= 0 & interprangeW < nnyq);

Suu_turb = Suu(interprange);
fmuu = fm(interprange);
try
Suu_turb (waverange) = []; Suu_turb (1) = [];
fmuu (waverange) = []; fmuu(1) = [];
Suu_turb = Suu_turb(fmuu>0);
fmuu = fmuu(fmuu>0);

Svv_turb = Svv(interprange);
fmvv = fm(interprange);
Svv_turb (waverange); Svv_turb (1) = [];
fmvv (waverange); fmvv(1) = [];
Svv_turb = Svv_turb(fmvv>0);
fmvv = fmvv(fmvv>0);

Sww_turb = Sww(interprangeW);
fmww = fm(interprangeW);
Sww_turb (waverange) = []; Sww_turb (1) = [];
fmww (waverange) = []; fmww(1) = [];
Sww_turb = Sww_turb(fmww>0);
fmww = fmww(fmww>0);


%% Linear interpolation over turbulent spectra

F = log(fmuu);
S = log(Suu_turb);
Puu = polyfit(F,S,1);
Puuhat = exp(polyval(Puu,log(fm)));clear F S

F = log(fmvv);
S = log(Svv_turb);
Pvv = polyfit(F,S,1);
Pvvhat = exp(polyval(Pvv,log(fm)));clear F S

F = log(fmww);
S = log(Sww_turb);
Pww = polyfit(F,S,1);
Pwwhat = exp(polyval(Pww,log(fm)));clear F S

%% Wave spectra

Suu_wave = Suu(waverange) - Puuhat(waverange)';
Svv_wave = Svv(waverange) - Pvvhat(waverange)';
Sww_wave = Sww(waverange) - Pwwhat(waverange)';

% Wave Fourier components
Amu_wave = sqrt((Suu_wave+0.*imag(Suu_wave))*(df));
Amv_wave = sqrt((Svv_wave+0.*imag(Svv_wave))*(df));
Amww_wave = sqrt((Sww_wave+0.*imag(Sww_wave))*(df));

% Wave Magnitudes
Um_wave = sqrt(real(Amu_wave).^2 + imag(Amu_wave).^2);
Vm_wave = sqrt(real(Amv_wave).^2 + imag(Amv_wave).^2);
wm_wave = sqrt(real(Amww_wave).^2 + imag(Amww_wave).^2);

% Wave reynolds stress
uw_wave = nansum(Um_wave.*wm_wave.*cos(Wph(waverange)-Uph(waverange)));
uv_wave = nansum(Um_wave.*Vm_wave.*cos(Vph(waverange)-Uph(waverange)));
vw_wave = nansum(Vm_wave.*wm_wave.*cos(Wph(waverange)-Vph(waverange)));

uu_wave = nansum(Suu_wave.*df);
vv_wave = nansum(Svv_wave.*df);
ww_wave = nansum(Sww_wave.*df);


% Full reynolds stresses
uu = nansum(real(Suu)*df);
uv = nansum(real(Suv)*df);
uw = nansum(real(Suw)*df);
vv = nansum(real(Svv)*df);
vw = nansum(real(Svw)*df);
ww = nansum(real(Sww)*df);

% Turbulent reynolds stresses

upup = uu - uu_wave;
vpvp = vv - vv_wave;
wpwp = ww - ww_wave;
upwp = uw - uw_wave;
upvp = uv - uv_wave;
vpwp = vw - vw_wave;

% Turbulent reynolds stresses
uw = upwp;
vw = vpwp;
uv = upvp;
uu = upup;
vv = vpvp;
ww = wpwp;


% Wave reynolds stresses
uw_wave = uw_wave;
vw_wave = vw_wave;
uv_wave = uv_wave;
uu_wave = uu_wave;
vv_wave = vv_wave;
ww_wave = ww_wave;

% tke=0.5*(uu.^2+vv.^2+ww.^2);
tke=0.5*(uu+vv+ww);
catch
upup=nan;
vpvp=nan;
wpwp=nan;
upwp = nan;
uw_wave = nan;
end