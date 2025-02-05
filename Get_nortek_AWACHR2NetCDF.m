function [AWAC] = Get_Nortek_AWACHR2NetCDF(varargin)
%% PROCESS AWAC-HR INSTRUMENT DATA
% This function processes Nortek Aquadopp Standard instrument data.  
%
% SYNTAX: [AWAC] = Get_Nortek_AWACHR2NetCDF(varargin)
% 
% VARARGIN PAIRS e.g.('zinst',0.3):
% data_dir             [-]    Directory with data                (default = pwd)
% instname             [-]    Name of file(s) from instrument    (default = AWAC01)
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
ncfilename  = varargin{6};
schema      = varargin{8};
zhead       = varargin{10};
zpressure   = varargin{12};
zblank      = varargin{14};
cell_size   = varargin{16};
burstmode   = varargin{18};

profile_interval    = '';
profile_avg_interval = '';

i_argin     = 1;

while i_argin<nargin - 0
switch lower ( varargin{i_argin})
    case 'data_dir';              i_argin=i_argin+1;       data_dir = varargin{i_argin};
    case 'instname';              i_argin=i_argin+1;       instname = varargin{i_argin}; 
    case 'ncfilename';            i_argin=i_argin+1;       ncfilename = varargin{i_argin};     
    case 'schema';                i_argin=i_argin+1;       schemaData = varargin{i_argin};
    case 'zhead';                 i_argin=i_argin+1;       zhead = varargin{i_argin};
    case 'zpressure';             i_argin=i_argin+1;       zpressure = varargin{i_argin};
    case 'zblank';                i_argin=i_argin+1;       zblank = varargin{i_argin};
    case 'cell_size';             i_argin=i_argin+1;       cell_size = varargin{i_argin};
    case 'profile_interval';      i_argin=i_argin+1;       profile_interval = varargin{i_argin};
    case 'profile_avg_interval';  i_argin=i_argin+1;       profile_avg_interval = varargin{i_argin};
    case 'burstmode';             i_argin=i_argin+1;       burstmode = varargin{i_argin};
otherwise
    error('Invalid string argument: %s.',varargin{i_argin});
end

i_argin = i_argin + 1;

end

if burstmode
    error 'Burst mode not implemented'
end

%% -----------------------------------
%  WRITE DATA OUT TO NetCDF or MATLAB
%  -----------------------------------
% Update and write the file schema
% schemaData=struct;
criar_NetCDF_File(data_dir, ncfilename, schemaData, {...
                                'date_created',        datestr(now, 'yyyy-mm-ddTHH:MM:SS'),...
                                'file_version',        'Level 0 Raw Data',...
                                'raw_datafile',        instname,...
                                'zhead',               zhead,...
                                'zpressure',           zpressure});
                             
%% PRELIMINARIES
% Turn on diary to record any processing information
delete([data_dir,'MatLabProcessing.log']);
diary([data_dir,'MatLabProcessing.log']);
AWAC = struct;

% Update data path to include instrument
datapath = [data_dir,'\',instname];

%% --------------------------
% READ AND WRITE SENSOR DATA
% ---------------------------
% Sensor File (*.sen)
display 'Loading sensor file...'
sen_raw = load([datapath '.sen']);

% DETERMINE SENSOR QUANTITIES
% Time
AWAC.tstamp_sensor = datenum(sen_raw(:,3), sen_raw(:,1), sen_raw(:,2), sen_raw(:,4), sen_raw(:,5), sen_raw(:,6));
    
% Temperature
AWAC.temperature = sen_raw(:,15);

% Heading
AWAC.heading = sen_raw(:,11);

% Pitch
AWAC.pitch = sen_raw(:,12);

% Roll
AWAC.roll = sen_raw(:,13);

% Pressure
AWAC.pressure = sen_raw(:,14);

% Determine the number of samples
sensorsamples = length(AWAC.tstamp_sensor);

% Write out variables to NetCDF
VarName = 'TIME_SENSOR_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'double',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'sensor measurement time stamp';
                                'units',               'days since 0000-01-01T00:00:00+8';
                                'local_time_zone',     'UTC+8';
                                'axis',                'T';
                                'valid min',           AWAC.tstamp_sensor(1);
                                'valid_max',           AWAC.tstamp_sensor(end);
                                'FillValue',           -999;
                                'calendar',            'gregorian';
                                'deployment',          1;
                                'time_start',          datestr(AWAC.tstamp_sensor(1),'yyyy-mm-ddTHH:MM:SS');
                                'time_end',            datestr(AWAC.tstamp_sensor(end),'yyyy-mm-ddTHH:MM:SS');
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, AWAC.tstamp_sensor);
                             
VarName = 'TEMPERATURE_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'temperature measurement';
                                'units',               'degrees celsius';
                                'valid min',           -30;
                                'valid_max',           60;
                                'FillValue',           -999;
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, single(AWAC.temperature)); 

VarName = 'HEADING_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'heading measurement';
                                'units',               'degrees from north';
                                'valid min',           0;
                                'valid_max',           360;
                                'FillValue',           -999;
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, single(AWAC.heading)); 
                             
VarName = 'PITCH_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'pitch measurement';
                                'units',               'degrees from horizontal';
                                'valid min',           0;
                                'valid_max',           360;
                                'FillValue',           -999;
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, single(AWAC.pitch)); 

VarName = 'ROLL_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'roll measurement';
                                'units',               'degrees';
                                'valid min',           0;
                                'valid_max',           360;
                                'FillValue',           -999;
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, single(AWAC.roll)); 

VarName = 'PRESSURE_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'pressure measurement';
                                'units',               'dBar';
                                'valid min',           0;
                                'valid_max',           5;
                                'FillValue',           -999;
                                'zpressure',           zpressure;
                                'comment',             'Raw data. Has not been adjusted for elevation above the bed'});
ncwrite(ncfilename, VarName, single(AWAC.pressure)); 


%% ----------------------------
% READ AND WRITE CONTINUOUS DATA
% -----------------------------

display 'Loading velocity file 1 of 3...'
AWAC.vx = load([datapath '.v1']);
AWAC.vx = single(AWAC.vx);
display 'Loading velocity file 2 of 3...'
AWAC.vy = load([datapath '.v2']);
AWAC.vy = single(AWAC.vy);
display 'Loading velocity file 3 of 3...'
AWAC.vz = load([datapath '.v3']);
AWAC.vz = single(AWAC.vz);
display 'Loading amplitude file 1 of 3...'
AWAC.ax = load([datapath '.a1']);
AWAC.ax = single(AWAC.ax);
display 'Loading amplitude file 2 of 3...'
AWAC.ay = load([datapath '.a2']);
AWAC.ay = single(AWAC.ay);
display 'Loading amplitude file 3 of 3...'
AWAC.az = load([datapath '.a3']);
AWAC.az = single(AWAC.az);
display 'Loading correlation file 1 of 3...'
% AWAC.cx = load([datapath '.c1']);
% AWAC.cx = single(AWAC.cx(:,3:end));
% display 'Loading correlation file 2 of 3...'
% AWAC.cy = load([datapath '.c2']);
% AWAC.cy = single(AWAC.cy(:,3:end));
% display 'Loading correlation file 3 of 3...'
% AWAC.cz = load([datapath '.c3']);
% AWAC.cz = single(AWAC.cz(:,3:end));

[velsamples, cells] = size(AWAC.vx);
AWAC.time = AWAC.tstamp_sensor;

% Height of bottom of first bin above seabed
zcell1 = zhead + zblank;

% Height of each bin from seabed (middle of bin)
AWAC.cells = [1:cells] * cell_size + zcell1;

% Write out variables to NetCDF
VarName = 'CELLS_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'CELLS_RAW', cells});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'middle of the cell location above bed';
                                'units',               'm';
                                'valid min',           0;
                                'valid_max',           100;
                                'FillValue',           -999;
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, single(AWAC.cells)); 

VarName = 'TIME_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'double',... 
                                'Dimensions',          {'TIME_RAW', velsamples});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'velocity measurement time stamp';
                                'units',               'days since 0000-01-01T00:00:00+8';
                                'local_time_zone',     'UTC+8';
                                'axis',                'T';
                                'valid min',           AWAC.time(1);
                                'valid_max',           AWAC.time(end);
                                'FillValue',           -999;
                                'calendar',            'gregorian';
                                'deployment',          1;
                                'time_start',          datestr(AWAC.time(1),'yyyy-mm-ddTHH:MM:SS');
                                'time_end',            datestr(AWAC.time(end),'yyyy-mm-ddTHH:MM:SS');
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, AWAC.time);

% Write "VELOCITY_X" variable and attributes
VarName = 'VEL_XE_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'Velocity of Beam 1 x or east (raw)';
                                'units',               'm/s';
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, AWAC.vx);

% Write "VELOCITY_Y" variable and attributes
VarName = 'VEL_YN_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'Velocity of Beam 2 y or north (raw)';
                                'units',               'm/s';  
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, AWAC.vy);

% Write "VELOCITY_Z" variable and attributes
VarName = 'VEL_ZU_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'Velocity of Beam 3 z or up (raw)';
                                'units',               'm/s';   
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, AWAC.vz);

% Write "AMPLITUDE_X" variable and attributes
VarName = 'AMP_XE_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'Amplitude of Beam 1 x or east (raw)';
                                'units',               'counts';
                                'valid min',           0;
                                'valid_max',           256;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, AWAC.ax);

% Write "AMPLITUDE_Y" variable and attributes
VarName = 'AMP_YN_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'Amplitude of Beam 2 y or north (raw)';
                                'units',               'counts';
                                'valid min',           0;
                                'valid_max',           256;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, AWAC.ay);

% Write "AMPLITUDE_Z" variable and attributes
VarName = 'AMP_ZU_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'Amplitude of Beam 3 z or up (raw)';
                                'units',               'counts';
                                'valid min',           0;
                                'valid_max',           256;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, AWAC.az);

% Write "CORRELATION_X" variable and attributes
VarName = 'CORR_XE_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'Correlation of Beam 1 x or east (raw)';
                                'units',               'percent';
                                'valid min',           0;
                                'valid_max',           100;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
% ncwrite(ncfilename, VarName, AWAC.cx);

% Write "CORRELATION_Y" variable and attributes
VarName = 'CORR_YN_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RWA', velsamples, 'CELLS_RAW', cells});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'Correlation of Beam 2 y or north (raw)';
                                'units',               'percent';
                                'valid min',           0;
                                'valid_max',           100;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
% ncwrite(ncfilename, VarName, AWAC.cy);

% Write "CORRELATION_Z" variable and attributes
VarName = 'CORR_ZU_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
escrever_NetCDF_atributos (ncfilename, VarName,...
                               {'long_name',           'Correlation of Beam 3 z or up (raw)';
                                'units',               'percent';
                                'valid min',           0;
                                'valid_max',           100;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
% ncwrite(ncfilename, VarName, AWAC.cz);
diary('off');
