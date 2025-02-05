function [ADP] = Get_Nortek_ADPHR2NetCDF(varargin)
%% PROCESS ADP-HR INSTRUMENT DATA
% This function processes Nortek Aquadopp Standard instrument data.  
%
% SYNTAX: [ADP] = Get_Nortek_ADPHR2NetCDF(varargin)
% 
% VARARGIN PAIRS e.g.('zinst',0.3):
% data_dir             [-]    Directory with data                (default = pwd)
% instname             [-]    Name of file(s) from instrument    (default = ADP01)
% ncfilename           [-]    NetCDF filename                    (required)
% schema               [-]    Schema for NetCDF file             (default = {})
% zhead                [m]    Height of head above bed           (default = 0)
% zpressure            [m]    Height of pressure sensor above bed(default = 0)
% zblank               [m]    Blanking distance                  (default = 0.1)
% cell_size            [m]    Cell size                          (default = 0.1)
% burstmode            [-]    Data is in burst mode              (default = 0)
%
% OPTIONAL SCHEMA DATA
% profile_interval     [s]    Profile interval                   (default = '')
% profile_avg_interval [s]    Profile interval                   (default = '')
%
% OUTPUT: instname.mat and ncfilename.nc
% tstamp            [-]      Time variable
% temperature       [deg C]  Temperature variable
% pressure          [dbar]   Pressure
% vx, vy, vz        [m/s]    Velocities in the X, Y and Z directions
% ve, vn, vu        [m/s]    Velocities in the east, north and up directions
% ax, ay, az        [counts] Signal strength in the X, Y and Z directions
%
% CODE CONTROL: A.Pomeroy (July 2014)

%% DEFINE DEFAULT VARIABLES WITHIN A VARIABLE STRUCTURE AND UPDATE
if isempty (varargin{1}), data_dir    = pwd; else, data_dir =  (varargin{1}); end
if isempty (varargin{2}), instname    = 'ADV';else, instname =  (varargin{2}); end
if isempty (varargin{3}), ncfilename  = '';else, ncfilename =  (varargin{3}); end
if isempty (varargin{4}), schema      = struct;else, schema =  (varargin{4}); end
if isempty (varargin{5}), zhead       = 0;else, zhead =  (varargin{5}); end
if isempty (varargin{6}), zpressure   = 0;else, zpressure =  (varargin{6}); end
if isempty (varargin{7}), zblank      = 0.1;else, zblank =  (varargin{7}); end
if isempty (varargin{8}), cell_size   = 0.1;else, cell_size =  (varargin{8}); end
if isempty (varargin{9}), burstmode   = 0;else, burstmode =  (varargin{9}); end
if isempty (varargin{10}), profile_interval    = 60;else, profile_interval =  (varargin{10}); end
if isempty (varargin{11}), profile_avg_interval = 58;else, profile_avg_interval =  (varargin{11}); end
if isempty (varargin{12}), sample_rate = 4; else, sample_rate =  (varargin{12}); end


%% -----------------------------------
%  WRITE DATA OUT TO NetCDF or MATLAB
%  -----------------------------------
% Update and write the file schema
Create_NetCDF_File(data_dir,ncfilename, schema, {...
                                'date_created',        datestr(now, 'yyyy-mm-ddTHH:MM:SS');
                                'file_version',        'Level 0 Raw Data';
                                'raw_datafile',        instname;
                                'zhead',               zhead;
                                'zpressure',           zpressure;
                                 });
                                
                                
%% PRELIMINARIES
% Turn on diary to record any processing information
% ncwriteschema([data_dir '\' ncfilename],schema);

% Add some more schema data
ncwriteatt(ncfilename, '/', 'instrument',                   'Nortek Vector');
ncwriteatt(ncfilename, '/', 'instrument_serial',            '3325');
ncwriteatt(ncfilename, '/', 'sampling_rate',                sprintf('%g Hz',sample_rate));
ncwriteatt(ncfilename, '/', 'Burst interval',               'CONTINUOUS');
ncwriteatt(ncfilename, '/', 'Samples per burst',            'N/A');   
ncwriteatt(ncfilename, '/', 'Nominal velocity range',       '1.00 m/s');
ncwriteatt(ncfilename, '/', 'Burst interval',               'CONTINUOUS');
ncwriteatt(ncfilename, '/', 'Samples per burst',            'N/A');
ncwriteatt(ncfilename, '/', 'Sampling volume',              '14.9 mm');
ncwriteatt(ncfilename, '/', 'transmit length',              '4.0 mm');
ncwriteatt(ncfilename, '/', 'receive length',               '0.01 m');
ncwriteatt(ncfilename, '/', 'heading',                      '');
ncwriteatt(ncfilename, '/', 'compass_calibration',          '1 degree error');
ncwriteatt(ncfilename, '/', 'compass_correlation',          '> 90 %');                                
ncwriteatt(ncfilename, '/', 'clock_reference',              'full memory. not applicable.');
ncwriteatt(ncfilename, '/', 'source',                       'observation');
ncwriteatt(ncfilename, '/', 'geospatial_lat',               'varies');
ncwriteatt(ncfilename, '/', 'geospatial_lat_units',         'degrees_north');
ncwriteatt(ncfilename, '/', 'geospatial_lon',               'varies');
ncwriteatt(ncfilename, '/', 'geospatial_lon_units',         'degrees_east');
ncwriteatt(ncfilename, '/', 'geospatial_vertical_velocity', '0.45');
ncwriteatt(ncfilename, '/', 'geospatial_vertical_pressure', 'N/A');
ncwriteatt(ncfilename, '/', 'geospatial_vertical_units',    'metres above bed');
ncwriteatt(ncfilename, '/', 'velocity_coordinates',         'XYZ');
ncwriteatt(ncfilename, '/', 'file_version',                 'raw data');
ncwriteatt(ncfilename, '/', 'comment',                      'deployment 1 of 2. position of frame varies see location movement metadata');

%% 3. DATA PROCESSING
% READ AND WRITE SENSOR DATA
ADV = struct;

% Sensor File (*.sen)
display 'Loading sensor file...'
sen_raw = load(fullfile(data_dir, [instname '.sen']));

% DETERMINE SENSOR QUANTITIES
% Time
ADV.tstamp_sensor = datenum(sen_raw(:,3), sen_raw(:,1), sen_raw(:,2), sen_raw(:,4), sen_raw(:,5), sen_raw(:,6));
    
% Temperature
ADV.temperature = sen_raw(:,14);

% Heading
ADV.heading = sen_raw(:,11);

% Pitch
ADV.pitch = sen_raw(:,12);

% Roll
ADV.roll = sen_raw(:,13);

% Determine the number of samples
sensorsamples = length(ADV.tstamp_sensor);

% Write out variables to NetCDF
VarName = 'TIME_SENSOR';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'double',... 
                                'Dimensions',          {'TIME_SENSOR', sensorsamples});
ncwriteatt(ncfilename, VarName, 'long_name',           'sensor measurement time stamp');
ncwriteatt(ncfilename, VarName, 'units',               'days since 0000-01-01T00:00:00+8');
ncwriteatt(ncfilename, VarName, 'local_time_zone',     'UTC+8');
ncwriteatt(ncfilename, VarName, 'valid_min',           ADV.tstamp_sensor(1));
ncwriteatt(ncfilename, VarName, 'valid_max',           ADV.tstamp_sensor(end));
ncwriteatt(ncfilename, VarName, 'calendar',            'gregorian');
ncwriteatt(ncfilename, VarName, 'time_start',          datestr(ADV.tstamp_sensor(1),'yyyy-mm-ddTHH:MM:SS'));
ncwriteatt(ncfilename, VarName, 'time_end',            datestr(ADV.tstamp_sensor(end),'yyyy-mm-ddTHH:MM:SS'));
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, ADV.tstamp_sensor);                             
                             
VarName = 'TEMPERATURE';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR', sensorsamples});
ncwriteatt(ncfilename, VarName, 'long_name',           'temperature measurement');
ncwriteatt(ncfilename, VarName, 'units',               'degrees celsius');
ncwriteatt(ncfilename, VarName, 'valid_min',           -30);
ncwriteatt(ncfilename, VarName, 'valid_max',           60);
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.temperature)); 

VarName = 'HEADING';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR', sensorsamples});
ncwriteatt(ncfilename, VarName, 'long_name',           'heading measurement');
ncwriteatt(ncfilename, VarName, 'units',               'degrees from north');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           360);
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.heading)); 
                             
VarName = 'PITCH';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR', sensorsamples});
ncwriteatt(ncfilename, VarName, 'long_name',           'pitch measurement');
ncwriteatt(ncfilename, VarName, 'units',               'degrees from horizontal');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           360);
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.pitch)); 

VarName = 'ROLL';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR', sensorsamples});
ncwriteatt(ncfilename, VarName, 'long_name',           'roll measurement');
ncwriteatt(ncfilename, VarName, 'units',               'degrees');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           360);
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.roll)); 

%% READ AND WRITE VELOCITY DATA
% Velocity Information (%.vhd)
display 'Loading velocity information...'
dat_info = load([instname '.vhd']);
time_raw = datenum(dat_info(:,3), dat_info(:,1), dat_info(:,2), dat_info(:,4), dat_info(:,5), dat_info(:,6));

% Velocity File (*.dat)
display 'Loading velocity files...'
dat_raw = load([instname '.dat']);

% Arrange data into burst
clc; display 'Processing velocity and amplitude data...'

% Determine the size of the data
bursts = max(dat_raw(:,1));

if bursts == 1
    display('Data is continuous')
    ADV.vx(1,:) = dat_raw(:,3);
    ADV.vy(1,:) = dat_raw(:,4);
    ADV.vz(1,:) = dat_raw(:,5);

    ADV.ax(1,:) = dat_raw(:,6);
    ADV.ay(1,:) = dat_raw(:,7);
    ADV.az(1,:) = dat_raw(:,8);

    ADV.snrx(1,:) = dat_raw(:,9);
    ADV.snry(1,:) = dat_raw(:,10);
    ADV.snrz(1,:) = dat_raw(:,11);

    ADV.cx(1,:) = dat_raw(:,12);
    ADV.cy(1,:) = dat_raw(:,13);
    ADV.cz(1,:) = dat_raw(:,14);

    ADV.pressure(1,:) = dat_raw(:,15);
else
    display('Data is in bursts')
    for bi = 1:size(dat_raw,1)
        ADV.vx(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,3);
        ADV.vy(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,4);
        ADV.vz(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,5);

        ADV.ax(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,6);
        ADV.ay(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,7);
        ADV.az(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,8);

        ADV.snrx(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,9);
        ADV.snry(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,10);
        ADV.snrz(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,11);

        ADV.cx(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,12);
        ADV.cy(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,13);
        ADV.cz(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,14);

        ADV.pressure(dat_raw(bi,1),dat_raw(bi,2)) = dat_raw(bi,15);
    end
end
        
% Determine the size of the data
bursts = size(ADV.vx,1);
samples = size(ADV.vx,2);

% Generate timeseries of time
% Correct time_raw length
if bursts < size(time_raw,1)
    time_raw = time_raw(1:end-1,:);
end

timestep = [1:1:samples].* datenum(0,0,0,0,0,1/sample_rate);
ADV.time = repmat(time_raw,1,samples)+ repmat(timestep,bursts,1);

% Write out variables to NetCDF
VarName = 'TIME';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'double',... 
                                'Dimensions',          {'BURSTS', bursts, 'TIME', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'measurement time stamp');
ncwriteatt(ncfilename, VarName, 'units',               'days since 0000-01-01T00:00:00+8');
ncwriteatt(ncfilename, VarName, 'local_time_zone',     'UTC+8');
ncwriteatt(ncfilename, VarName, 'valid_min',           ADV.time(1,1));
ncwriteatt(ncfilename, VarName, 'valid_max',           ADV.time(end,end));
ncwriteatt(ncfilename, VarName, 'calendar',            'gregorian');
ncwriteatt(ncfilename, VarName, 'deployment',          1);
ncwriteatt(ncfilename, VarName, 'time_start',          datestr(ADV.time(1,1),'yyyy-mm-ddTHH:MM:SS'));
ncwriteatt(ncfilename, VarName, 'time_end',            datestr(ADV.time(end,end),'yyyy-mm-ddTHH:MM:SS'));
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, ADV.time);

% Write "PRESSURE" variable and attributes
VarName = 'PRES';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'PRES', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'pressure measurements');
ncwriteatt(ncfilename, VarName, 'units',               'dBar');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           100);       
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.pressure));

% Write "VELOCITY_X" variable and attributes
VarName = 'VEL_X';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'VEL', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'velocity measurements in x direction');
ncwriteatt(ncfilename, VarName, 'units',               'm/s');
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.vx));

% Write "VELOCITY_Y" variable and attributes
VarName = 'VEL_Y';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'VEL', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'velocity measurements in y direction');
ncwriteatt(ncfilename, VarName, 'units',               'm/s');  
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.vy));

% Write "VELOCITY_Z" variable and attributes
VarName = 'VEL_Z';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'VEL', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'velocity measurements in z direction');
ncwriteatt(ncfilename, VarName, 'units',               'm/s');   
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.vz));

% Write "AMPLITUDE_X" variable and attributes
VarName = 'AMP_X';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'SAMPLES', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'amplitude measurements in x direction');
ncwriteatt(ncfilename, VarName, 'units',               'counts');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           256);       
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.ax));

% Write "AMPLITUDE_Y" variable and attributes
VarName = 'AMP_Y';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'SAMPLES', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'amplitude measurements in x direction');
ncwriteatt(ncfilename, VarName, 'units',               'counts');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           256);       
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.ay));

% Write "AMPLITUDE_Z" variable and attributes
VarName = 'AMP_Z';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'SAMPLES', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'amplitude measurements in x direction');
ncwriteatt(ncfilename, VarName, 'units',               'counts');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           256);       
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.az));

% Write "SN_RATIO_X" variable and attributes
VarName = 'SNR_X';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'SAMPLES', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'signal to noise ratio in x direction');
ncwriteatt(ncfilename, VarName, 'units',               'dB');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           50);       
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.snrx));

% Write "SN_RATIO_Y" variable and attributes
VarName = 'SNR_Y';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'SAMPLES', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'signal to noise ratio in y direction');
ncwriteatt(ncfilename, VarName, 'units',               'dB');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           50);       
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.snry));

% Write "SN_RATIO_Z" variable and attributes
VarName = 'SNR_Z';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'SAMPLES', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'signal to noise ratio in z direction');
ncwriteatt(ncfilename, VarName, 'units',               'dB');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           50);       
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.snrz));

% Write "CORRELATION_X" variable and attributes
VarName = 'CORR_X';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'SAMPLES', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'correlation in x direction');
ncwriteatt(ncfilename, VarName, 'units',               'percent');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           100);       
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.cx));

% Write "CORRELATION_Y" variable and attributes
VarName = 'CORR_Y';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'SAMPLES', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'correlation in y direction');
ncwriteatt(ncfilename, VarName, 'units',               'percent');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           100);       
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.cy));

% Write "CORRELATION_Z" variable and attributes
VarName = 'CORR_Z';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'BURSTS', bursts, 'SAMPLES', samples});
ncwriteatt(ncfilename, VarName, 'long_name',           'correlation in z direction');
ncwriteatt(ncfilename, VarName, 'units',               'percent');
ncwriteatt(ncfilename, VarName, 'valid_min',           0);
ncwriteatt(ncfilename, VarName, 'valid_max',           100);       
ncwriteatt(ncfilename, VarName, 'sample_rate',         sample_rate);
ncwriteatt(ncfilename, VarName, 'sample_rate_units',   'Hz');
ncwriteatt(ncfilename, VarName, 'comment',             'raw data');
ncwrite(ncfilename, VarName, single(ADV.cz));
