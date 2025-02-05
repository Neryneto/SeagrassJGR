function [SM1,f1]=calculatespectrum(Fs,eta,Nven,frequencies)

%Input:
% Fs: sampling frequency
% eta: free surface zero-crossing variation
% Nven: window size

T = 1/Fs;              

x=(eta); % my current data

Xd=detrend(x); %hourly data
Nover=round(Nven/2); %overlap - 50%

[SM1,f1] = pwelch(Xd(1:2^(nextpow2(length(Xd))-1)),Nven,Nover,frequencies,Fs); %2^nextpow2(length(Xd))

end