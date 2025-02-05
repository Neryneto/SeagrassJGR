% Linear dispersion relation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k=lindisp_explicit(f,h)

% Dispersion relation approximation based on Hunt (1979)
d1=0.6666666666;
d2=0.3555555555;
d3=0.1608465608;
d4=0.0632098765;
d5=0.0217540484;
d6=0.0065407983;

y=(2*pi*f).^2.*h/9.81;
num=1+d1*y+d2*y.^2+d3*y.^3+d4*y.^4+d5*y.^5+d6*y.^6;

k=1./h.*sqrt(y.^2+y./(num));