%% write bathymetry

fulldata={};
minx=[];miny=[];
maxx=[];maxy=[];
for i=1:14
%     plot(fulldata{1,i}(:,1),fulldata{1,i}(:,2),'.')
%     hold on
    if i<10
        fulldata{i}=importdata(strcat('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Perth_LiDAR_2009\1448040',num2str(i),'_ahd.pts'));
        minx=[minx;min(fulldata{1,i}(:,1))];
        maxx=[maxx;max(fulldata{1,i}(:,1))];
        miny=[miny;min(fulldata{1,i}(:,2))];
        maxy=[maxy;max(fulldata{1,i}(:,2))];
    else
        fulldata{i}=importdata(strcat('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Perth_LiDAR_2009\144804',num2str(i),'_ahd.pts'));
        minx=[minx;min(fulldata{1,i}(:,1))];
        maxx=[maxx;max(fulldata{1,i}(:,1))];
        miny=[miny;min(fulldata{1,i}(:,2))];
        maxy=[maxy;max(fulldata{1,i}(:,2))];
    end
end
% lon = min (fulldata{1,5}(:,1)):100: max (fulldata{1,5}(:,1));
% lat = min (fulldata{1,5}(:,2)):100: max (fulldata{1,5}(:,2));
path='C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model'

%full grid
lon = 348000:10:365347.5;
lat = 6447715:10:6463551;

%bickley bay
lon2 = 363077.5:5:364767.5;
lat2 = 6456729.5:5:6458619.5;

%for map
lon3 = 3632007.5:5:3650197.5; 
lat3 = 6456500.5:5:6458570.5;

veg_lon = 363473:10:364241;
veg_lat = 6458441:5:6458619.5;

[xq,yq] = meshgrid(lon,lat);
[xq2,yq2] = meshgrid(lon2,lat2);

vq4 = griddata(fulldata{1,4}(:,1),fulldata{1,4}(:,2),fulldata{1,4}(:,3),xq,yq);
vq5 = griddata(fulldata{1,5}(:,1),fulldata{1,5}(:,2),fulldata{1,5}(:,3),xq,yq);

vq2 = griddata(fulldata{1,5}(:,1),fulldata{1,5}(:,2),fulldata{1,5}(:,3),xq2,yq2);

vq(vq>0)=nan;%vq(vq<-34)=-34;
contourf(lon,lat,vq)
hold on
contourf(lon,lat,vq4)

plot(349752,6448018,'.k','markersize',50)
dlmwrite('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\bat_bickley_coarse.dat',vq,' ');

bat_old=dlmread('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awac_forcing_spec_noFric_statio\bat_bickley_coarse.dat');
bat_fake=bat_old;

bat_fake(bat_old>=-0.04)=0;
bat_fake(bat_old<-0.04)=-3;

ind=bat_old==0;
bat_fake(ans)=-3;
dlmwrite('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awac_forcing_spec_noFric_statio\bat_bickley_fake_3m.dat',bat_fake,' ');

%% prepare vegetation grid
T = readtable('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Seagrass_assessment\results.xlsx','Sheet','SeagrassGridSimpl');
T_fake = readtable('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Seagrass_assessment\results.xlsx','Sheet','SeagrassGridSimplFake');

%bickley bay
lon2 = min(T_fake.Longitude):10:max(T_fake.Longitude);
lat2 = min(T_fake.Latitude):10:max(T_fake.Latitude);

[xq2,yq2] = meshgrid(lon2,lat2);

vq = griddata(T_fake.Longitude,T_fake.Latitude,T.density,xq2,yq2);
vq_fake = griddata(T_fake.Longitude,T_fake.Latitude,T_fake.density,xq2,yq2);
vq(isnan(vq))=0; vq_fake(isnan(vq_fake))=0;
dlmwrite('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\awacforcingVeg\veggrid.dat',floor(vq),' ');
dlmwrite('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVeg\veggrid.dat',floor(vq_fake),' ');
dlmwrite('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\fakehbv\awacforcingVegcd8\veggrid.dat',floor(vq_fake),' ');

%% write TPAR file Aquadopp
load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\Aquadopps\AQD2985\AQD2985.mat');
local='C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\aqdforcing';
fi2=fopen(fullfile(local,'AQD2985.txt'),'wt');

hm0=AQD2985.wav.Hm0;
nanx = isnan(AQD2985.wav.Hm0);
t    = 1:numel(AQD2985.wav.Hm0);
hm0(nanx) = interp1(t(~nanx), AQD2985.wav.Hm0(~nanx), t(nanx));

%monta o arquivo .bnd
fprintf(fi2,'%s\n','TPAR');
for i=1:length(AQD2985.wav.time)
    fprintf(fi2,'%s\t',datestr(AQD2985.wav.time(i),'yyyymmdd.HHMMSS')); 
    fprintf(fi2,['%1.4f\t%2.4f\t%3.4f\t%1.4f\n'],[hm0(i) AQD2985.wav.Tp(i) AQD2985.wav.DirTp(i) 2]) ;
end
fclose all

%
local='C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_8m';
fi2=fopen(fullfile(local,'AWAC8m.txt'),'wt');

hm0=AWAC_5523.wav.Hm0;
nanx = isnan(AWAC_5523.wav.Hm0);
t    = 1:numel(AWAC_5523.wav.Hm0);
hm0(nanx) = interp1(t(~nanx), AWAC_5523.wav.Hm0(~nanx), t(nanx));

%monta o arquivo .bnd
fprintf(fi2,'%s\n','TPAR');
for i=1:length(AWAC_5523.wav.time)
    fprintf(fi2,'%s\t',datestr(AWAC_5523.wav.time(i),'yyyymmdd.HHMMSS')); 
    fprintf(fi2,['%1.4f\t%2.4f\t%3.4f\t%1.4f\n'],[hm0(i) AWAC_5523.wav.Tp(i) AWAC_5523.wav.DirTp(i) 2]) ;
end
fclose all

%% write TPAR file RBR3

load('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\RBRs\Rottnest_RBRs_model_with_CD.mat');
local='C:\Users\22371812\Desktop\aqdforcing';
RBR.ST3.Tp(RBR.ST3.Tp>20)=nan;
hm0=RBR.ST3.Hsig_swell_out';
nanx = isnan(RBR.ST3.Hsig_swell_out');
t    = 1:numel(RBR.ST3.Hsig_swell_out');
hm0(nanx) = interp1(t(~nanx), RBR.ST3.Hsig_swell_out(~nanx)', t(nanx));

time_interp=datenum(2019,3,12,11,30,0):1/24/2:RBR.ST3.time_wave(end);
hm0_interp=interp1(RBR.ST3.time_wave,RBR.ST3.Hm0swell,time_interp);
tp_interp=interp1(RBR.ST3.time_wave,RBR.ST3.Tpswell,time_interp);

%monta o arquivo .txt
fi2=fopen(fullfile(local,'rbr3_input.txt'),'wt');
fprintf(fi2,'%s\n','TPAR');
for i=1:length(RBR.ST3.time_wave)
    fprintf(fi2,'%s\t',datestr(time_interp(i),'yyyymmdd.HHMMSS')); 
    fprintf(fi2,['%1.4f\t%2.4f\t%3.4f\t%1.4f\n'],[hm0_interp(i) tp_interp(i) 135 2]) ;
end
fclose all
model=load('C:\Users\22371812\Desktop\aqdforcing\rbr2.dat');
plot(datenum(num2str(model(:,1),'%f'),'yyyymmdd.HHMMSS'),model(:,2))
hold on
plot(RBR.ST2.time_wave,RBR.ST2.Hm0swell)

%% Prepare spectral AWAC
% Rotina desenvolvida para transformar os espectros de saída do AWAC em
% espectros 2D do SWAN
%%
local='C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\Instruments_data\1.Part_1\AWACs\AWAC_15m';
% latlonN='-40.193 -20.257';
% latlonE='-40.207 -20.339';
% latlonS='-40.275 -20.391';
declinacao_magnetica=-1.5; 

% Wave full directional spectra file *.wds
wds = load(fullfile(local,'AWACN101.wds'));
% Wave engery spectra file *.was
was = load(fullfile(local,'AWACN101.was'));
was(:,end)=[];
% whr: arquivo com os metadados
whr = load(fullfile(local,'AWACN101.sen'));
%% Abre o arquivo .whr e define o tempo inicial e o tempo final.
tempo=datenum(whr(:,3),whr(:,1),whr(:,2),whr(:,4),whr(:,5),whr(:,6));
tini=tempo(1);
tfim=tempo(end);
%% Organize data
was(was == -9.000000) = 0;
wds(wds == -9.000000) = 0;

w.freq = was(1,:)';
w.was= was(2:end,:);
w.dir = (0:4:356)+declinacao_magnetica;

w.dir=w.dir';
limite=max(find(w.dir<=0));
limite2=min(find(w.dir>=360));

w.dir(w.dir>=360)=w.dir(w.dir>=360)-360;
w.dir(w.dir<=0)=w.dir(w.dir<=0)+360;
w.dir=sort(w.dir);

w.wds = wds;

if ~isempty(limite)
    w.wds=[w.wds(:,limite+1:end) w.wds(:,1:limite)];
end
if ~isempty(limite2)
    w.wds=[w.wds(:,limite2:end) w.wds(:,1:limite2-1)];
end

nFreqs   = length(w.freq);
w.freq=w.freq(1:nFreqs,:);
nDirs    = length(w.dir);

nSamples = length(w.was);
w.wds = permute(reshape(w.wds', nDirs, nFreqs, nSamples), [ 2 1 3]);
w.t=tini:1/24/2:tfim;
w.t(end)=[];
% tempos=[tempos; w.t'];
% tempos=tempo;
%% Unnormalize data
% E(f,theta) = S(f) * d(f,theta)     ,   with S= Energy spectra (*.was)    and d=full directional spectra (*.wds)
E=zeros(nFreqs,nDirs,nSamples);
for tt=1:nSamples
    for f=1:nFreqs;
        for theta=1:nDirs;
            E(f,theta,tt) = (w.was(tt,f) * w.wds(f,theta,tt))/4;%divide por 4 para distribuir a energia na direção de cada bin
        end
        
    end
end
fi2=fopen(fullfile(local,'AWACN101_spc2.txt'),'wt');

%monta o arquivo .bnd
fprintf(fi2,'%s\n','SWAN 1 Swan standard spectral file, version');
fprintf(fi2,'%s\n','$ Data produced by SWAN version 40.51');
fprintf(fi2,'%s\n','$ Project:’Rotto’ ; run number:’spectral vegetation’');
fprintf(fi2,'%s\n','TIME                                    time-dependent data');
fprintf(fi2,'%s\n','     1                                  time coding option');
fprintf(fi2,'%s\n','LOCATIONS locations in x-y-space');
fprintf(fi2,'%s\n','1 	number of locations');
fprintf(fi2,'%s\n','364498 6456752');
fprintf(fi2,'%s\n','RFREQ relative frequencies in Hz');
fprintf(fi2,'%s\n',[num2str(length(w.freq)) ' 	number of frequencies']);
tamanho=size(w.was,1)*size(w.was,2);
for j=1:length(w.freq)
    fprintf(fi2,'%15.4f\n',w.freq(j));
end

fprintf(fi2,'%s\n','NDIR spectral Nautical directions in degr');
fprintf(fi2,'%s\n',[num2str(length(w.dir)) ' 	number of directions']);

for j=1:length(w.dir)
    fprintf(fi2,'%15.4f\n',(w.dir(j)));
end
fprintf(fi2,'%s\n','QUANT');
fprintf(fi2,'%s\n','     1                                  number of quantities in table');
fprintf(fi2,'%s\n','VaDens                                  energy densities in m2/Hz/degr');
fprintf(fi2,'%s\n','m2/Hz/degr                            unit');
fprintf(fi2,'%s\n','   -0.9900E+02                          exception value');

for i=1:size(E,3)
%     tamanho=length(arquivo);
    fprintf(fi2,'%s\n',datestr(w.t(i), 'yyyymmdd.HHMMSS'));
    fprintf(fi2,'%s\n','FACTOR');
    fprintf(fi2,'%s\n','1');
    fmt=repmat([repmat('  %7i',1,90) '\n'],1,98);
    for l=1:size(E,1)
        %         fprintf(fi2,fmt,(squeeze(E(l,:,i))))
        fprintf(fi2,fmt,(squeeze(E(l,:,i))));
    end
end
            
fclose all

%% RUN ALL CDS

for CD=1:8
    for h=1:4
%         fid=fopen(['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\cdh2\cd' num2str(CD) 'h' num2str(h) '\awacforcing.swn']);
% i = 1;
% tline = fgetl(fid);
% A{i} = tline;
% while ischar(tline)
%     i = i+1;
%     tline = fgetl(fid);
%     A{i} = tline;
% end
% A{40} = sprintf('VEGEtation 1 0.%d 0.015 1 0.%d',h,CD);
% % A{26} = sprintf('BOUN SIDE E VAR FILE 800.00 ''tpar_awac.txt'' 1'
% % A{27} = sprintf('BOUN SIDE S VAR FILE 910.00 ''tpar_awac2.txt'' 1'
% % A{28} = sprintf('BOUN SIDE W VAR FILE 910.00 ''tpar_awac3.txt'' 1'
% % A{69} = sprintf('%d',99);
% 
% fid = fopen(['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\cdh2\cd' num2str(CD) 'h' num2str(h) '\awac3bound.swn'], 'w');
% for i = 1:numel(A)
%     if A{i+1} == -1
%         fprintf(fid,'%s', A{i});
%         break
%     else
%         fprintf(fid,'%s\n', A{i});
%     end
% end
%         copyfile('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\cdh2\cd1h1\awacforcing.swn',...
%             ['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\cdh2\cd' num2str(CD) 'h' num2str(h)])
%          copyfile('C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\cdh2\cd1h1\tpar_awac3.txt',...
%             ['C:\Users\22371812\OneDrive - The University of Western Australia\Documents\PhD\chapters\Rottnest\model\cdh2\cd' num2str(CD) 'h' num2str(h)])

       cd(fullfile(['C:\Users\22371812\Desktop\cdh2\cd' num2str(CD) 'h' num2str(h)]))
       a=strcat({'swanrun awac3bound'});
       [s w]=dos(cell2mat(a));
    end
end
