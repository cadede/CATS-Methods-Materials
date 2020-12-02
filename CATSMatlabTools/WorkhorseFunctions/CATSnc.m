function		ncfile = CATSnc(prhfile,tagguideloc,ncloc)
% Creates NC file from Tag Guide and PRH file for CATS Tags
% Saves outputfile in same directory as PRH file
% Danuta Wisniewska, updated by James Fahlbusch
% Version 2.22.18
% Goldbogen Lab
% Stanford University
%
% format:
% ncfile = CATSnc(prhfile,tagguideloc, ncStructures);
% ncfile = CATSnc(prhfile,tagguideloc);
% ncfile = CATSnc(prhfile)
% ncfile = CATSnc()
%
% returns: name of ncfile created
% ncfile = CATSnc(...); 
%
% Workflow ToDo
% Tag Guide Updates: Make sure TagType and TagNum are consistent, 
% Add Life History Stage (Adult, Juv) Column
% Add Location of Tag on Body
% Add Provider_Name and Provider_Organization
% 
% nargin = 0; %uncomment this line if you are running from the code (not as a function)
if nargin <3;
    if nargin == 2; 
        % Load the NC Metadata Structures
        % look for Structures folder in same directory as this script
        p = mfilename('fullpath');
        i = strfind(p,'\');
        p = p(1:i(end));
        ncloc = [p 'Structures\'];
        clearvars p i;
        if 7~= exist(ncloc, 'dir')
            [~,ncloc]=uigetfile('*.*', 'Select any file in the Structures folder (check Utilities folder)');
        end
    end 
    if nargin == 1; 
        [filename,fileloc]=uigetfile('*.*', 'Select the Tag Guide');
        tagguideloc = [fileloc filename];
        % Load the NC Metadata Structures
        % look for Structures folder in same directory as this script
        p = mfilename('fullpath');
        i = strfind(p,'\');
        p = p(1:i(end));
        ncloc = [p 'Structures\'];
        clearvars p i;
        if 7~=exist(ncloc, 'dir')
            [~,ncloc]=uigetfile('*.*', 'Select any file in the Structures folder');
        end
    end
    if nargin == 0; 
        [prhfilename,prhfileloc]=uigetfile('*.*', 'Select the PRH File');
        prhfile = [prhfileloc prhfilename];
        [filename,fileloc]=uigetfile('*.*', 'Select the Tag Guide');
        tagguideloc = [fileloc filename];
        % Load the NC Metadata Structures
        % look for Structures folder in same directory as this script
        try
            p = mfilename('fullpath');
            i = strfind(p,'\');
            p = p(1:i(end));
            ncloc = [p 'Structures\'];
            clearvars p i;        
            if 7~=exist(ncloc, 'dir')
                [~,ncloc]=uigetfile('*.*', 'Select any file in the Structures folder');
            end
        catch
            [~,ncloc]=uigetfile('*.*', 'Select any file in the Structures folder');
        end

    end
end
%% 1. Load files
i = strfind(prhfile,'\');
prhfileloc = prhfile(1:i(end));
prhfilename = prhfile(i(end)+1:end);
ii = strfind(tagguideloc,'\');
fileloc = tagguideloc(1:ii(end));
filename = tagguideloc(ii(end)+1:end);
iii = strfind(prhfilename,' ');
whalename = prhfilename(1:iii-1);
clearvars i ii iii;
[num,txt,raw] = xlsread(fullfile(fileloc,filename),1);
ak=strmatch('ID',txt(3,:));
species=strmatch('Spec',txt(3,:));
row = strmatch(whalename,txt(:,ak));

disp('Loading PRH data, will take some time'); 
load([prhfileloc prhfilename(1:end-3) 'mat']); 

disp('Section 1 (Load Files) finished');
%% 2. Populate NC Info Fields from Tag Guide
%  Load default tag info from Setup CSV
prefix = txt{row,ak};
info=csv2struct([ncloc 'CATSsetup']);
if isfield(info,'field')
    info = rmfield(info,'field');
end
info.depid = prefix;

% Source File information
info.dtype_source = [filename ', ' prhfilename];
info.dtype_nfiles = num2str(2);
info.dtype_format = [upper(filename(strfind(filename,'.')+1:end)) ', ' upper(prhfilename(strfind(prhfilename,'.')+1:end))];
info.dtype_author = 'Dave Cade, davecade@stanford.edu';
dbfile = dir(fullfile(prhfileloc,prhfilename));
info.dtype_datetime_made = datestr(dbfile.date,'dd-mmm-yyyy HH:MM:SS');

% Populate fields from Tag Guide
kk = strmatch('Tag #',txt(3,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk}))
    info.device_serial = num2str(raw{row,kk}); 
else
    info.device_serial = 'Unknown';
end

kk = strmatch('Tag Type',txt(3,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk}))
    info.device_model = raw{row,kk};
else
    info.device_model = 'Unknown';
end

% Determine Sensors from PRH
disp('Determining sensors from PRH');
sensorList = '3 Axis Accelerometer, 3 Axis Magnetometer, 3 Axis Gyroscope, Pressure, Light-level, Temperature (External)';
if exist('INFO.TempInternal','var') == 1 || exist('Temp1', 'var') == 1
  sensorList = strcat(sensorList, ', Temperature (Internal)');
end
if exist('GPS','var') == 1
  sensorList = strcat(sensorList, ', Position (GPS)');
end
%perhaps a second check for camon to see if there is a camera, etc
info.sensors_list = sensorList;
fprintf(' Sensors: %s \n',sensorList);
clearvars sensorList

 
kk = strmatch('Biopsy',txt(3,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk})) & (~strcmpi('N',raw{row,kk}))
    if (strcmpi('Y',raw{row,kk}))
        info.animal_lifehist_sex = 'Biopsy taken, results unknown';
    else
        info.animal_lifehist_sex = raw{row,kk};
    end
else
    info.animal_lifehist_sex = 'Unknown';
end

kk = strmatch('UTC',txt(3,:));
if (~isempty(kk)) && (~isnan(raw{row,kk})) && (~strcmpi('NA',raw{row,kk}))
    info.dephist_device_tzone = raw{row,kk};
    info.dephist_utc2loc = raw{row,kk};
else
    info.dephist_device_tzone = 0;
    info.dephist_utc2loc = 0;
end

kk = strmatch('Study_Area',txt(3,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk})) & (~strcmpi('NA',raw{row,kk}))
    info.dephist_deploy_locality = raw{row,kk};
else
    info.dephist_deploy_locality = 'Unknown';
end 

%kk = strmatch('Attachment',txt(3,:));
%if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk})) & (~strcmpi('NA',raw{row,kk}))
%    info.dephist_deploy_method = raw{row,kk};
%else
%    info.dephist_deploy_method = 'Unknown';
%end

kk = strmatch('Animal',txt(3,:));
if ~isnan(raw{row,kk})
    if strcmpi(raw{row,kk},'U')| isempty(raw{row,kk})
        info.animal_id = 'Unknown';
    else
        info.animal_id = num2str(raw{row,kk});
    end
end
%Deployment Information
kk = strmatch('Lat_On',txt(3,:));
if ~isnan(raw{row,kk})
    if strcmpi(raw{row,kk},'U')| isempty(raw{row,kk})
        info.dephist_deploy_location_lat = 'Unknown';
    elseif ~isnumeric(raw{row,kk})
        disp('Check TagGuide Lat_On entry; using Unknown');
        info.dephist_deploy_location_lat = 'Unknown';
    else
        info.dephist_deploy_location_lat = raw{row,kk};
    end
end
kk = strmatch('Long_On',txt(3,:));
if ~isnan(raw{row,kk})
    if strcmpi(raw{row,kk},'U')| isempty(raw{row,kk})
        info.dephist_deploy_location_lon = 'Unknown';
    elseif ~isnumeric(raw{row,kk})
        disp('Check TagGuide Long_On entry; using Unknown');
        info.dephist_deploy_location_lon = 'Unknown';
    else
        info.dephist_deploy_location_lon = raw{row,kk};
    end
end
% Use PRH Data Start and Data End
info.dephist_device_datetime_start = datestr(DN(1)-(info.dephist_device_tzone/24),'dd-mmm-yyyy HH:MM:SS.fff');
info.dephist_device_datetime_end = datestr(DN(end)-(info.dephist_device_tzone/24),'dd-mmm-yyyy HH:MM:SS.fff');

%Check tagon and tagoff times against tag guide
kk = strmatch('Tag_On',txt(3,:));
if ~isnan(datenum(raw{row,kk})) && ~isempty(raw{row,kk}) 
    d = datenum(raw{row,kk});
    %convert to UTC
    d = d-(info.dephist_device_tzone/24);
    %Check to see if the Tag guide start is before the device start; if so, use tag guide value for tag_on
    if DN(1)-(info.dephist_device_tzone/24) > datenum(d)       
        disp(['Using TagGuide TagOnAnimal: ' datestr(d, 'dd-mmm-yyyy HH:MM:SS.fff') ' (UTC)']);
        info.dephist_deploy_datetime_start = datestr(d, 'dd-mmm-yyyy HH:MM:SS.fff');
        disp('TagOn before Device Start');
    else
        disp(['Using PRH TagOnAnimal: ' datestr(DN(find(tagon,1))-(info.dephist_device_tzone/24), 'dd-mmm-yyyy HH:MM:SS.fff') ' (UTC)']);
        info.dephist_deploy_datetime_start = datestr(DN(find(tagon,1))-(info.dephist_device_tzone/24),'dd-mmm-yyyy HH:MM:SS.fff');
    end   
else
    disp(['Using PRH TagOnAnimal (no entry in Tag Guide): ' datestr(DN(find(tagon,1))-(info.dephist_utc2loc/24), 'dd-mmm-yyyy HH:MM:SS.fff') ' (UTC)']);
    info.dephist_deploy_datetime_start = datestr(DN(find(tagon,1))-(info.dephist_utc2loc/24),'dd-mmm-yyyy HH:MM:SS.fff');
end
kk = strmatch('Tag_Off',txt(3,:));
if ~isnan(datenum(raw{row,kk})) && ~isempty(raw{row,kk}) 
    d = datenum(raw{row,kk});
    %convert to UTC
    d = d-(info.dephist_utc2loc/24);
    if DN(end)-(info.dephist_device_tzone/24) < datenum(d)       
        disp(['Using TagGuide TagOffAnimal: ' datestr(d, 'dd-mmm-yyyy HH:MM:SS.fff') ' (UTC)']);
        info.dephist_deploy_datetime_end = datestr(d, 'dd-mmm-yyyy HH:MM:SS.fff');
        disp('TagOff after Device End');
    else
        disp(['Using PRH TagOffAnimal: ' datestr(DN(find(tagon,1,'last'))-(info.dephist_device_tzone/24), 'dd-mmm-yyyy HH:MM:SS.fff') ' (UTC)']);
        info.dephist_deploy_datetime_end = datestr(DN(find(tagon,1,'last'))-(info.dephist_device_tzone/24),'dd-mmm-yyyy HH:MM:SS.fff');
    end
else
    disp(['Using PRH TagOff Animal (no entry in Tag Guide): ' datestr(DN(find(tagon,1,'last'))-(info.dephist_device_tzone/24), 'dd-mmm-yyyy HH:MM:SS.fff') ' (UTC)']);
    info.dephist_deploy_datetime_end = datestr(DN(find(tagon,1,'last'))-(info.dephist_device_tzone/24),'dd-mmm-yyyy HH:MM:SS.fff');
end

% Recovery Information
kk = strmatch('Recover_Lat',txt(3,:));
if ~isnan(raw{row,kk})
    if strcmpi(raw{row,kk},'U')| isempty(raw{row,kk})
        info.dephist_recov_location_lat = 'Unknown';
    elseif ~isnumeric(raw{row,kk})
        disp('Check TagGuide Recover_Lat entry; using Unknown');
        info.dephist_recov_location_lat = 'Unknown';
    else
        info.dephist_recov_location_lat = raw{row,kk};
    end
end
kk = strmatch('Recover_Long',txt(3,:));
if ~isnan(raw{row,kk})
    if strcmpi(raw{row,kk},'U')| isempty(raw{row,kk})
        info.dephist_recov_location_lon = 'Unknown';
    elseif ~isnumeric(raw{row,kk})
        disp('Check TagGuide Recover_Long entry; using Unknown');
        info.dephist_recov_location_lon = 'Unknown';
    else
        info.dephist_recov_location_lon = raw{row,kk};
    end
end
kk = strmatch('Recover_Time',txt(3,:));
if ~isnan(raw{row,kk})
    if isnumeric(datenum(raw{row,kk}))
        d = datenum(raw{row,kk});
        %convert to UTC
        d = d-(info.dephist_device_tzone/24);       
        info.dephist_recov_datetime_start = datestr(d, 'dd-mmm-yyyy HH:MM:SS.fff');
        clearvars d;
    else
        disp(['Recovery entry: ' raw{row,kk} ' not a date or in the wrong format. Recovery Date set to Unknown']);
        info.dephist_recov_datetime_start = 'Unknown';
    end
else
    disp('No entry in Tag Guide for Recovery Date; set to Unknown');
    info.dephist_recov_datetime_start = 'Unknown';
end
 
% Populate species names
[S,~]=read_csv([ncloc 'tagged_species.csv'],1) ;
k = strncmpi(txt{row,species},{S.abbrev},2) ;
if all(k==0),
    fprintf(' Warning: unknown species %s. Set metadata manually\n', prefix) ;
    return ;
end
if sum(k)>1,
    fprintf(' More than one species match "%s". Retry with a longer species name.\n',prefix) ;
    return
end
info.animal_species_code = S(k).abbrev;
info.animal_species_common = S(k).common_name;
info.animal_species_science = S(k).science_name;

% Add User Defined Metadata (Drone, Video)
kk = strmatch('Drone',txt(3,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk})) & (~strcmpi('n/a',raw{row,kk}))
   info.udm_drone = raw{row,kk};
else
   info.udm_drone = 'Not Applicable';
end
kk = strmatch('Primary Camera',txt(3,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk})) & (~strcmpi('n/a',raw{row,kk}))
   info.udm_camera_primaryview = raw{row,kk};
else
   info.udm_camera_primaryview = 'Not Applicable';
end
kk = strmatch('Total Video',txt(3,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk})) & (~strcmpi('NO VID',raw{row,kk}))
    if isnumeric(raw{row,kk})
        info.udm_camera_videotime = datestr(raw{row,kk}, 'HH:MM:SS');
    else
        disp('Check Total Video Time entry in TagGuide');
        info.udm_camera_videotime = 'Unknown';
    end
else
   info.udm_camera_videotime = 'Not Applicable';
end

% Add Project and Provider details from Tag Guide
kk = strmatch('Project    _',txt(3,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk}))
    info.project_name = raw{row,kk};
else
    info.project_name = 'Unknown';
    fprintf(' Warning: project %s. Set metadata manually\n', info.project_name) ;
end
kk = strmatch('PI Contact',txt(3,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk})) & (~strcmpi('n/a',raw{row,kk}))
    info.provider_email = raw{row,kk};
else
    info.provider_email = 'Unknown';
end
kk = strmatch('Project Dates',txt(3,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk})) & (~strcmpi('n/a',raw{row,kk}))
    p = raw{row,kk};
    i = strfind(p,'-');
    info.project_datetime_start = strtrim(p(1:i-1));
    info.project_datetime_end = strtrim(p(i+1:end));
    clearvars p i;
else
    info.project_datetime_start = '';
    info.project_datetime_end = '';
end

%Get Provider details from DataProviders tab of TagGuide
[num,txt,raw] = xlsread(fullfile(fileloc,filename), 'DataProviders');
ak=strmatch('provider email',txt(1,:));
row = strmatch(info.provider_email,txt(:,ak));
if isempty(row)
    disp(['Check DataProviders Tab in TagGuide, ' info.provider_email ' not found']); 
end
kk = strmatch('provider name',txt(1,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk})) & (~strcmpi('n/a',raw{row,kk}))
    info.provider_name = raw{row,kk};
else
    info.provider_name = 'Unknown';
end
kk = strmatch('provider details',txt(1,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk})) & (~strcmpi('n/a',raw{row,kk}))
    info.provider_details = raw{row,kk};
else
    info.provider_details = 'Unknown';
end
kk = strmatch('provider license',txt(1,:));
if (~isempty(kk)) & (~isempty(row)) & (~isnan(raw{row,kk})) & (~strcmpi('n/a',raw{row,kk}))
    info.provider_license = raw{row,kk};
else
    info.provider_license = 'Unknown';
end

info.creation_date = datestr(now,'dd-mmm-yyyy HH:MM:SS');

disp('Section 2 (Create Info Structure) finished');
%% 3. Create NC file 
% Converts the PRH File into NC format and combines with info struct
ncfile = convert_prh_CATS([prhfileloc prhfilename],info ,prhfileloc);
disp('Section 3 (Create NC file) finished');

%% 4. Add Jiggle Speed to NC file

speedJJ = csv2struct([ncloc 'nc_addSpeedJiggle']);
%remove 1 extraneous line from import 
if isfield(speedJJ,'field')
    speedJJ = rmfield(speedJJ,'field');
end
% Populate the sensor struct with Speed data from PRH
speedJJ.depid = prefix ;
speedJJ.creation_date = datestr(now,'dd-mmm-yyyy HH:MM:SS');
speedJJ.data = speed.JJ;
speedJJ.sampling_rate = fs;
try speedJJ.speed_calibration_periods_end_times = round(speedstats.sections_end_index/fs); 
catch; speedJJ.speed_calibration_periods_end_times = INFO.tagslip.SpeedPeriods/fs;
end
try speedJJ.fit_r2 = speedstats.r2used.JJr2;
catch; for ii = 1:size(speedstats.Fit,1); speedJJ.fit_r2(ii) = speedstats.Fit{ii,1}.rsquare; end
end
speedJJ.speed_jiggleRMS = speed.JRMS;
speedJJ.start_offset = 0;
speedJJ.history = 'sens_struct';
speedJJ.author = 'Dave Cade, davecade@stanford.edu';

add_nc([prhfileloc ncfile],speedJJ) ;
disp('Section 4 (Add Jiggle Speed) finished');

%% 5. Add Speed from FN

speedFN = csv2struct([ncloc 'nc_addSpeedFN']);
%remove 1 extraneous line from import 
if isfield(speedFN,'field')
    speedFN = rmfield(speedFN,'field');
end

speedFN.depid = prefix ;
speedFN.creation_date = datestr(now,'dd-mmm-yyyy HH:MM:SS');
speedFN.data = speed.FN;
speedFN.sampling_rate = fs;
try speedFN.speed_calibration_periods_end_times = round(speedstats.sections_end_index/fs); 
catch; speedFN.speed_calibration_periods_end_times = INFO.tagslip.SpeedPeriods/fs;
end
try speedFN.fit_r2 = speedstats.r2used.FNr2;
catch; for ii = 1:size(speedstats.Fit,1); speedJJ.fit_r2(ii) = speedstats.Fit{ii,2}.rsquare; end
end
speedFN.speed_flownoiseRMS = flownoise;
speedFN.start_offset = 0;
speedFN.history = 'sens_struct';
speedFN.author = 'Dave Cade, davecade@stanford.edu';

add_nc([prhfileloc ncfile],speedFN) ;
disp('Section 5 (Add Speed from FN) finished');

%% 6. Plot NC
clearvars -except prhfileloc ncfile
oi = load_nc([prhfileloc ncfile]);
load_nc([prhfileloc ncfile]);
pitch = oi.pitch;
[ax,h]=plott(P, pitch, roll, head, speedJJ);
a=findobj(ax,'type','axe');
b = get( get(a(1),'YLabel') );
set( get(a(1),'YLabel'), 'String', 'Depth' );
b = get( get(a(2),'YLabel') );
set( get(a(2),'YLabel'), 'String', 'Pitch' );
b = get( get(a(3),'YLabel') );
set( get(a(3),'YLabel'), 'String', 'Roll' );
b = get( get(a(4),'YLabel') );
set( get(a(4),'YLabel'), 'String', 'Heading' );
b = get( get(a(5),'YLabel') );
set( get(a(5),'YLabel'), 'String', 'Speed' );

disp('Section 6 (Plot) finished');
