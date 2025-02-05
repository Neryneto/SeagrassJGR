% Calculation of the drag coefficient: MATLAB script
% The corresponding coding of the drag coefficient was numerically evaluated using the following script:
function CD=calculate_Cd (Lp,bv,hv,Nv,T,Ho,Hi,ho,hi)
%% General parameters
g=9.8; 			%out(g,'m/s/s, gravitational acceleration')
Nrdx=1000; 		%out (Nrdx,'-, # of integration steps')
%% Geometry, fluid and forcing parameters (from measurements)
TW=1;           %out(TW,'-, 1 for regular and 2 for irregular)
% Lp=1;			%out(Lp,'m, vegetation meadow length')
% bv=[0.0022];	%out(bv,'m, vegetation stem widths')
% hv=[0.172];		%out(hv,'m, vegetation stem heights')
% Nv=[2600];		%out(Nv,'-/m/m, vegetation stem densities')
rho=1025;		%out(rho,'kg/m/m/m, sea water density')
% T=0.4;			%out(T,'s, peak wave period')
% Ho=0.024;		%out(Ho,'m, Hrms, upwave')
% Hi=0.018;		%out(Hi,'m, Hrms, downwave')
% ho=  0.06;		%out(ho,'m, mean water depth, upwave')
% hi=  0.06;		%out(hi,'m, mean water depth, downwave')
if ~(numel(bv)==numel(Nv) && numel(Nv)==numel(hv))
    error('number of vegetation segments must be the same, check bv, Nv and hv')
end

%% Calculation grid definition
dx=Lp/Nrdx;		%out(dx,'m, integration distance')
Nr=Lp/dx+1;		%out(Nr,'-, Number of cells in calculations')
x=0:dx:Lp;
xc=x-dx/2;xc=xc(2:end);
h=interp1([0 Lp],[ho hi],x);
%% Wave parameters – calculating k from kh using deep water dispersion relation
koh=4*pi^2/g/T^2.*h;
kh=koh.*(1-exp(-1*koh.^(5/4))).^(-2/5);	% Goh 2002 to estimate kh
% Newton-Raphson method to solve for k
for j=1:20
    AA=tanh(kh);
    kh=kh-(kh.*AA-koh)./(AA+koh.*(1-AA.^2));
end
k=kh./h;		
kh1=kh(1);		%out(kh1,'-, wave number * depth')
L=2*pi./k;
c=L/T;
cg=1/2*c.*(1+2*kh./sinh(2*kh));		
cg1=cg(1);		%out(cg1,'m/s, group velocity')
omega=2*pi./T;

%% Vegetation dissipation model grid constants
A=1/8*rho*g*diff(cg)/dx;
C=dx/(rho*g/8)./cg;
hvs=[0 cumsum(hv)];
B=k*0;
for j=1:numel(bv)
    BB=sinh(k.*hvs(j+1)).^3+3*sinh(k.*hvs(j+1))-...
      (sinh(k.*hvs(j  )).^3+3*sinh(k.*hvs(j  )));
    B=B+bv(j)*Nv(j).*BB;
end
 B=2/3/pi*rho*(k*g/2/omega).^3/3./k./cosh(kh).^3.*B;
if TW==1 	   	 %  Regular waves
 
  %% model as is
elseif TW==2		%  Irregular  waves, modification from Mendez et al. 2004
    B=B * 3/4*sqrt(pi);
else
    error('Unknown wave type TW=%g',TW)
end

%% Searching for drag coefficient (solnHi calls ‘fminsearch’ function)
CD=1;
CD=solnHi(Ho,Hi,A,B,C,x,CD);	%out(CD,'-, Estimated CD') 
H2=zeros(size(x))+nan; 
H2(1)=Ho^2; 
for i=1:numel(x)-1 
    H2(i+1)=H2(i)-C(i)*(CD*B(i)*H2(i)^(3/2)+A(i)*H2(i)); 
end

%% fminsearch function 
function CDe=solnHi(Ho,Hi,A,B,C,x,iCD)
CDe=fminsearch(@calcHi,iCD);
    function errHi=calcHi(CD)
        H2=zeros(size(x))+nan;
        H2(1)=Ho^2;
        for i=1:numel(x)-1
            H2(i+1)=H2(i)-C(i)*(CD*B(i)*H2(i)^(3/2)+A(i)*H2(i));
        end
        errHi=(sqrt(H2(end))-Hi)^2;
    end
end
end
