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
if isempty (varargin{2}), instname    = 'ADP01';else, instname =  (varargin{2}); end
if isempty (varargin{3}), ncfilename  = '';else, ncfilename =  (varargin{3}); end
if isempty (varargin{4}), schema      = '';else, schema =  (varargin{4}); end
if isempty (varargin{5}), zhead       = 0;else, zhead =  (varargin{5}); end
if isempty (varargin{6}), zpressure   = 0;else, zpressure =  (varargin{6}); end
if isempty (varargin{7}), zblank      = 0.1;else, zblank =  (varargin{7}); end
if isempty (varargin{8}), cell_size   = 0.1;else, cell_size =  (varargin{8}); end
if isempty (varargin{9}), burstmode   = 0;else, burstmode =  (varargin{9}); end
if isempty (varargin{10}), profile_interval    = 60;else, profile_interval =  (varargin{10}); end
if isempty (varargin{11}), profile_avg_interval = 58;else, profile_avg_interval =  (varargin{11}); end


% i_argin     = 1;
% 
% while i_argin<nargin - 0
% switch lower ( varargin{i_argin})
%     case 'data_dir';              i_argin=i_argin+1;       data_dir = varargin{i_argin};
%     case 'instname';              i_argin=i_argin+1;       instname = varargin{i_argin}; 
%     case 'ncfilename';            i_argin=i_argin+1;       ncfilename = varargin{i_argin};     
%     case 'schema';                i_argin=i_argin+1;       schemaData = varargin{i_argin};
%     case 'zhead';                 i_argin=i_argin+1;       zhead = varargin{i_argin};
%     case 'zpressure';             i_argin=i_argin+1;       zpressure = varargin{i_argin};
%     case 'zblank';                i_argin=i_argin+1;       zblank = varargin{i_argin};
%     case 'cell_size';             i_argin=i_argin+1;       cell_size = varargin{i_argin};
%     case 'profile_interval';      i_argin=i_argin+1;       profile_interval = varargin{i_argin};
%     case 'profile_avg_interval';  i_argin=i_argin+1;       profile_avg_interval = varargin{i_argin};
%     case 'burstmode';             i_argin=i_argin+1;       burstmode = varargin{i_argin};
% otherwise
%     error('Invalid string argument: %s.',varargin{i_argin});
% end
% 
% i_argin = i_argin + 1;
% 
% end
% 
% if burstmode
%     error 'Burst mode not implemented'
% end

%% -----------------------------------
%  WRITE DATA OUT TO NetCDF or MATLAB
%  -----------------------------------
% Update and write the file schema
% Create_NetCDF_File(fullfile(data_dir,ncfilename), schema, {...
%                                 'date_created',        datestr(now, 'yyyy-mm-ddTHH:MM:SS');
%                                 'file_version',        'Level 0 Raw Data';
%                                 'raw_datafile',        instname;
%                                 'zhead',               zhead;
%                                 'zpressure',           zpressure;
%                                  });
                             
%% PRELIMINARIES
% Turn on diary to record any processing information
% delete([data_dir,'MatLabProcessing.log']);
% diary([data_dir,'MatLabProcessing.log']);
ADP = struct;

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
ADP.tstamp_sensor = datenum(sen_raw(:,3), sen_raw(:,1), sen_raw(:,2), sen_raw(:,4), sen_raw(:,5), sen_raw(:,6));
    
% Temperature
ADP.temperature = sen_raw(:,14);

% Heading
ADP.heading = sen_raw(:,11);

% Pitch
ADP.pitch = sen_raw(:,14);

% Roll
ADP.roll = sen_raw(:,15);

% Pressure
ADP.pressure = sen_raw(:,16);

% Determine the number of samples
sensorsamples = length(ADP.tstamp_sensor);

% Write out variables to NetCDF
VarName = 'TIME_SENSOR_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'double',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'sensor measurement time stamp';
                                'units',               'days since 0000-01-01T00:00:00+8';
                                'local_time_zone',     'UTC+8';
                                'axis',                'T';
                                'valid min',           ADP.tstamp_sensor(1);
                                'valid_max',           ADP.tstamp_sensor(end);
                                'FillValue',           -999;
                                'calendar',            'gregorian';
                                'deployment',          1;
                                'time_start',          datestr(ADP.tstamp_sensor(1),'yyyy-mm-ddTHH:MM:SS');
                                'time_end',            datestr(ADP.tstamp_sensor(end),'yyyy-mm-ddTHH:MM:SS');
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, ADP.tstamp_sensor);
                             
VarName = 'TEMPERATURE_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'temperature measurement';
                                'units',               'degrees celsius';
                                'valid min',           -30;
                                'valid_max',           60;
                                'FillValue',           -999;
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, single(ADP.temperature)); 

VarName = 'HEADING_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'heading measurement';
                                'units',               'degrees from north';
                                'valid min',           0;
                                'valid_max',           360;
                                'FillValue',           -999;
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, single(ADP.heading)); 
                             
VarName = 'PITCH_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'pitch measurement';
                                'units',               'degrees from horizontal';
                                'valid min',           0;
                                'valid_max',           360;
                                'FillValue',           -999;
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, single(ADP.pitch)); 

VarName = 'ROLL_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'roll measurement';
                                'units',               'degrees';
                                'valid min',           0;
                                'valid_max',           360;
                                'FillValue',           -999;
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, single(ADP.roll)); 

VarName = 'PRESSURE_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'TIME_SENSOR_RAW', sensorsamples});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'pressure measurement';
                                'units',               'dBar';
                                'valid min',           0;
                                'valid_max',           5;
                                'FillValue',           -999;
                                'zpressure',           zpressure;
                                'comment',             'Raw data. Has not been adjusted for elevation above the bed'});
ncwrite(ncfilename, VarName, single(ADP.pressure)); 


%% ----------------------------
% READ AND WRITE CONTINUOUS DATA
% -----------------------------

display 'Loading velocity file 1 of 3...'
ADP.vx = load([datapath '.v1']);
ADP.vx = single(ADP.vx(:,3:end));
display 'Loading velocity file 2 of 3...'
ADP.vy = load([datapath '.v2']);
ADP.vy = single(ADP.vy(:,3:end));
display 'Loading velocity file 3 of 3...'
ADP.vz = load([datapath '.v3']);
ADP.vz = single(ADP.vz(:,3:end));
display 'Loading amplitude file 1 of 3...'
ADP.ax = load([datapath '.a1']);
ADP.ax = single(ADP.ax(:,3:end));
display 'Loading amplitude file 2 of 3...'
ADP.ay = load([datapath '.a2']);
ADP.ay = single(ADP.ay(:,3:end));
display 'Loading amplitude file 3 of 3...'
ADP.az = load([datapath '.a3']);
ADP.az = single(ADP.az(:,3:end));
display 'Loading correlation file 1 of 3...'
% ADP.cx = load([datapath '.c1']);
% ADP.cx = single(ADP.cx(:,3:end));
% display 'Loading correlation file 2 of 3...'
% ADP.cy = load([datapath '.c2']);
% ADP.cy = single(ADP.cy(:,3:end));
% display 'Loading correlation file 3 of 3...'
% ADP.cz = load([datapath '.c3']);
% ADP.cz = single(ADP.cz(:,3:end));

[velsamples, cells] = size(ADP.vx);
ADP.time = ADP.tstamp_sensor;

% Height of bottom of first bin above seabed
zcell1 = zhead + zblank;

% Height of each bin from seabed (middle of bin)
ADP.cells = [1:cells] * cell_size + zcell1;

% Write out variables to NetCDF
VarName = 'CELLS_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',... 
                                'Dimensions',          {'CELLS_RAW', cells});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'middle of the cell location above bed';
                                'units',               'm';
                                'valid min',           0;
                                'valid_max',           100;
                                'FillValue',           -999;
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, single(ADP.cells)); 

VarName = 'TIME_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'double',... 
                                'Dimensions',          {'TIME_RAW', velsamples});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'velocity measurement time stamp';
                                'units',               'days since 0000-01-01T00:00:00+8';
                                'local_time_zone',     'UTC+8';
                                'axis',                'T';
                                'valid min',           ADP.time(1);
                                'valid_max',           ADP.time(end);
                                'FillValue',           -999;
                                'calendar',            'gregorian';
                                'deployment',          1;
                                'time_start',          datestr(ADP.time(1),'yyyy-mm-ddTHH:MM:SS');
                                'time_end',            datestr(ADP.time(end),'yyyy-mm-ddTHH:MM:SS');
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, ADP.time);

% Write "VELOCITY_X" variable and attributes
VarName = 'VEL_XE_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'Velocity of Beam 1 x or east (raw)';
                                'units',               'm/s';
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, ADP.vx);

% Write "VELOCITY_Y" variable and attributes
VarName = 'VEL_YN_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'Velocity of Beam 2 y or north (raw)';
                                'units',               'm/s';  
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, ADP.vy);

% Write "VELOCITY_Z" variable and attributes
VarName = 'VEL_ZU_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'Velocity of Beam 3 z or up (raw)';
                                'units',               'm/s';   
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, ADP.vz);

% Write "AMPLITUDE_X" variable and attributes
VarName = 'AMP_XE_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'Amplitude of Beam 1 x or east (raw)';
                                'units',               'counts';
                                'valid min',           0;
                                'valid_max',           256;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, ADP.ax);

% Write "AMPLITUDE_Y" variable and attributes
VarName = 'AMP_YN_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'Amplitude of Beam 2 y or north (raw)';
                                'units',               'counts';
                                'valid min',           0;
                                'valid_max',           256;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, ADP.ay);

% Write "AMPLITUDE_Z" variable and attributes
VarName = 'AMP_ZU_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'Amplitude of Beam 3 z or up (raw)';
                                'units',               'counts';
                                'valid min',           0;
                                'valid_max',           256;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
ncwrite(ncfilename, VarName, ADP.az);

% Write "CORRELATION_X" variable and attributes
VarName = 'CORR_XE_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'Correlation of Beam 1 x or east (raw)';
                                'units',               'percent';
                                'valid min',           0;
                                'valid_max',           100;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
% ncwrite(ncfilename, VarName, ADP.cx);

% Write "CORRELATION_Y" variable and attributes
VarName = 'CORR_YN_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RWA', velsamples, 'CELLS_RAW', cells});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'Correlation of Beam 2 y or north (raw)';
                                'units',               'percent';
                                'valid min',           0;
                                'valid_max',           100;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
% ncwrite(ncfilename, VarName, ADP.cy);

% Write "CORRELATION_Z" variable and attributes
VarName = 'CORR_ZU_RAW';
nccreate(ncfilename,   VarName, 'DeflateLevel',        9,...
                                'Datatype',            'single',...
                                'Dimensions',          {'TIME_RAW', velsamples, 'CELLS_RAW', cells});
Write_NetCDF_Attributes (ncfilename, VarName,...
                               {'long_name',           'Correlation of Beam 3 z or up (raw)';
                                'units',               'percent';
                                'valid min',           0;
                                'valid_max',           100;       
                                'FillValue',           -999;
                                'profile_interval',    profile_interval;
                                'profile_avg_interval',profile_avg_interval;
                                'profile_interval',    's';
                                'comment',             'raw data'});
% ncwrite(ncfilename, VarName, ADP.cz);
diary('off');
