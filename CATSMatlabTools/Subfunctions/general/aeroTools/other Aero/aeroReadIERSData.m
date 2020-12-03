function file =  aeroReadIERSData(folder,varargin)
%   AEROREADIERSDATA Create MAT-file containing current International
%   Astronomical Union (IAU) 2000A Earth orientation data
%
%   FILE = AEROREADIERSDATA(FOLDER) creates a MAT-file (aeroiersdataYYYYMMDD.mat),
%   based on IAU 2000A Earth orientation data from the International Earth
%   Rotation and Reference Systems Service (IERS). This function saves the
%   file to FOLDER. The MAT-file name has the format
%   aeroiersdataYYYYMMDD.mat, where:
%
%   * YYYY - Year
%   * MM - Month
%   * DD - Day
%
%   FILE = AEROREADIERSDATA(FOLDER,URL,'URLADDRESS') creates a MAT-file
%   (aeroiersdataYYYMMDD.mat), based on Earth orientation data from a
%   specified website.
%
%   Input arguments:
%
%   FOLDER: Folder in which the generated MAT-file is saved. Before
%   running this function, create the folder with write permission.
%
%   URL: Optional website or file that contains the IAU 2000A Earth
%   orientation data, for example,
%   'https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt'.
%   The default is 'http://maia.usno.navy.mil/ser7/finals2000A.data'.
%
%   Output arguments:
%   FILE: Character array containing the full path location of the new
%   MAT-file.
%
%   Examples:
%
%   Generate the aeroiersdataYYYMMDD.mat file for the current day, in the
%   current folder, from http://maia.usno.navy.mil/ser7/finals2000A.data.
%
%   file = aeroReadIERSData(pwd)
%
%   Generate the aeroiersdataYYYMMDD.mat file for the current day, in the
%   current folder, from the alternate website
%   https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt. 
%
%   file = aeroReadIERSData(pwd,'url',...
%       'https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt')
%
%   Generate the aeroiersdataYYYMMDD.mat file for the current day, in the
%   current folder, from finals2000A.data obtained from the IERS website
%   file located at C:\Documents.
%
%   file = aeroReadIERSData(pwd,'url','file:///C:\Documents\finals2000A.data')
%
%   See also DCMECI2ECEF, LLA2ECI, ECI2LLA, ECI2AER, DELTAUT1, MJULIANDATE,
%   DELTACIP, POLARMOTION.

%   Copyright 2017-2018 The MathWorks, Inc.

%   References: 
%   [1] IERS website. http://www.iers.org/ 
%   [2] United States Naval Observatory website. http://www.usno.navy.mil/

% Validate entry data
narginchk(1,3);
p = inputParser;
addRequired(p,'folder',@(x) validateattributes(x,{'char','string'},{'nonempty'}));
addParameter(p,'url','http://maia.usno.navy.mil/ser7/finals2000A.data',...
    @(x) validateattributes(x,{'char','string'},{'nonempty'}));
parse(p,folder,varargin{:});

% Folder write permissions check
[status,msg,~] =  fileattrib(folder);
if ~status || msg.UserWrite~=1
    error(message('aero:aerout1utc:invalidFolder',char(folder)));
end

% Read the URL
t = urlread(p.Results.url,'Timeout',30);

% Extract data strings
formatSpec = '%*6s%9s%2s%10s%9s%10s%9s%3s%10s%10s%8s%7s%3s%10s%9s%10s%9s%[^\n\r]';
readText=textscan(t, formatSpec, inf,'Delimiter', '', 'WhiteSpace', '', ...
    'TextType', 'string', 'EmptyValue', NaN,'ReturnOnError', false, ...
    'EndOfLine', '\r\n');

% Process text data
mjdStr = strtrim(readText{1});
PMipStr = strtrim(readText{2});
PMxStr = strtrim(readText{3});
PMxErrStr = strtrim(readText{4});
PMyStr = strtrim(readText{5});
PMyErrStr = strtrim(readText{6});
ut1ipStr = strtrim(readText{7});
ut1Str = strtrim(readText{8});
ut1ErrStr = strtrim(readText{9});
dXYipStr = strtrim(readText{12});
dXStr = strtrim(readText{13});
dXErrStr = strtrim(readText{14});
dYStr = strtrim(readText{15});
dYErrStr = strtrim(readText{16});

% Do a pre-check of the data
if isempty(mjdStr) || isempty(PMipStr) || isempty(PMxStr) || isempty(PMyStr) || ...
   isempty(ut1ipStr) || isempty(ut1Str) || isempty(dXYipStr) || ...
   isempty(dXStr) || isempty(dYStr) || isempty(PMxErrStr) || isempty(PMyErrStr) ||...
   isempty(ut1ErrStr) || isempty(dXErrStr) || isempty(dYErrStr)
   error(message('aero:aerout1utc:invalidFormat',p.Results.url));
end
if ~strcmp(PMipStr{1},'I') || ~strcmp(ut1ipStr{1},'I') || ~strcmp(dXYipStr{1},'I')
   error(message('aero:aerout1utc:invalidFormat',p.Results.url));
end
for k = 1:length(mjdStr)
    if ~isempty(PMipStr{k}) || ~isempty(ut1ipStr{k}) || ~isempty(dXYipStr{k})
        % Modified Julian Date
        mjd(k,1) = str2double(mjdStr{k}); %#ok<*AGROW,*NASGU>
        
        % Polar Motion [PMx PMy]
        if ~isempty(PMipStr{k})
            pmip{k,1} = PMipStr{k}; 
            pm(k,:) = [str2double(PMxStr{k}) str2double(PMyStr{k})];
            pmerr(k,:) = [str2double(PMxErrStr{k}) str2double(PMyErrStr{k})];
        end
        
        % UT1-UTC
        if ~isempty(ut1ipStr{k})
            ut1utcip{k,1} = ut1ipStr{k};
            ut1utc(k,1) = str2double(ut1Str{k});
            ut1utcerr(k,1) = str2double(ut1ErrStr{k});
        end
        
        % Nutation/dCIP [dX dY]
        if ~isempty(dXYipStr{k})
            dxyip{k,1} = dXYipStr{k};
            dxy(k,:) = [str2double(dXStr{k}) str2double(dYStr{k})];
            dxyerr(k,:) = [str2double(dXErrStr{k}) str2double(dYErrStr{k})];
        end
    else
        break;
    end
end

% Validate format of the data
for k=1:10
    % Check that the first 10 values for the read data have "I" for the I/P
    % characteristic. This is to make sure the format is correct.
    if (~strcmp(ut1utcip(k),'I') || ~strcmp(dxyip(k),'I')||...
            ~strcmp(pmip(k),'I'))
       error(message('aero:aerout1utc:invalidFormat',p.Results.url));
    end
end
% Save data
datestr = strrep(char(datetime('now','Format','uuuu-MM-dd')),'-','');
file =fullfile(p.Results.folder,['aeroiersdata' datestr '.mat']); 
try
    save(file,'mjd','ut1utc','ut1utcip','pmip','pm','dxyip','dxy','ut1utcerr',...
        'pmerr','dxyerr');
catch ME
    if strcmp(ME.identifier,'MATLAB:save:permissionDenied')
        error(message('aero:aerout1utc:invalidFolder',char(folder)));
    else
        throw(ME);
    end
end