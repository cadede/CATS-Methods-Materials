function varargout = moonLibration(ephemerisTime,varargin)
%   MOONLIBRATION Implement JPL development ephemeris for Moon libration
%   Euler angles
%   [angles,rates] = MOONLIBRATION(EPHEMERISTIME, MODEL, ACTION)
%   implements the mathematical representation of the moon librations based
%   on the Chebyshev coefficients calculated by NASA's Jet Propulsion
%   Laboratory (JPL). It calculates the Euler angles for the Moon's
%   librations.
%   
%   Inputs required by MOONLIBRATION are:
%   EPHEMERISTIME  :Julian dates for which the librations are calculated. It
%                   can be a scalar or column vector with M elements.
%                   Alternatively, the date can be given as a 2 element
%                   vector or M by 2 matrix. The first element can be a
%                   fixed Julian date and the second element can be the
%                   elapsed time between the fixed Julian date and the
%                   ephemeris time.
%   MODEL           :String defining the Development Ephemeris (DE) by
%                   JPL. It can be: '405', '421', '423', '430' or '432t'.
%                   If it is not defined, the default model is '405'.
%   ACTION          :String to determine action for out of range
%                   EPHEMERISTIME. Specify if out of range input invokes a
%                   'Warning', 'Error', or no action ('None'). The default
%                   is 'Error'.
%
%   Outputs calculated by MOONLIBRATION are:
%   angles  :Numeric array describing the three Euler angles for the 
%            librations. Its size is M by 3,with M rows, one for each
%            ephemeris time.The Euler angles are calculated in radians. 
%   rates   :Numeric array describing the three Euler angular rates for the 
%            librations. Its size is M by 3,with M rows, one for each
%            ephemeris time.The Euler angles are calculated in radians/day. 
%
%   Examples:
%
%   Calculate the libration angles of the Moon for December 1st 1990 with
%   DE405:
%       angles = moonLibration(juliandate(1990,12,1))
%
%   Calculate the libration angles and rates for the Moon for noon January 
%   1st, 2000 using DE421:
%       [angles,rates] = moonLibration([2451544.5 0.5],'421')
%
%   Note: This function uses ephemeris data that can be obtained using the
%   aeroDataPackage command.
%
%   See also PLANETEPHEMERIS, EARTHNUTATION.

%   Copyright 2012-2018 The MathWorks, Inc.

%   References:
%   [1]Folkner, W. M., J.G. Williams, D.H. Boggs, "The Planetary and
%      Lunar Ephemeris DE 421", JPL Interplanetary Network Progress Report
%   24-178, 2009.
%   [2]Ma, C. et al., "The International Celestial Reference Frame as
%      Realized by Very Long Baseline Interferometry", Astronomical Journal,
%      Vol. 116, 1998.
%   [3]Vallado, D. A., Fundamentals of Astrodynamics and Applications,
%      McGraw-Hill, New York, 1997.

%% Input check
narginchk(1,3);
nargoutchk(0,2);
% Valid inputs
validDatabases = {'405','421','423','430','432t'};
validError = {'Error','Warning','None'};

% Validate units and databases depending on number of inputs
switch nargin
    case 1
        de = '405';
        errorFlag = 'error';
    case 2
        de = varargin{1};
        % Validate database
        if isnumeric(de)
            de = num2str(de);
        end
        de = lower(validatestring(de,validDatabases));
        errorFlag = 'error';
    case 3
        de = varargin{1};
        errorFlag = varargin{2};
        % Validate database
        if isnumeric(de)
            de = num2str(de);
        end
        de = lower(validatestring(de,validDatabases));
        % Validate error flag
        errorFlag = lower(validatestring(errorFlag,validError));
end

% Julian Date
validateattributes(ephemerisTime,{'numeric'},{'nonempty','nonnan','finite'})
if size(ephemerisTime,2)>2
    error(message('aero:aeroephemerides:julianDateVector'));
end

%% Data loading
% Load Constants
try	 load(['ephConstants' de]);
 catch err	 
     if strcmp(err.identifier,'MATLAB:load:couldNotReadFile')	 
         error(message('aero:aeroephemerides:unavailableDatabase',...	 
                       '<a href="matlab:aeroDataPackage">aeroDataPackage</a>'));	 
     else	 
         rethrow(err)	 
     end
 end
% Load Objects
objLib = matfile(['ephLibrations' de '.mat']);
bufLib = [];
% Initialize warning flags
warnFlag = [0 0];
%% Output Initialization
switch nargout
    case {0,1} % Only position requested
        angles = zeros(size(ephemerisTime,1),3);
        velFlag = false;
    case 2 % Position and velocity requested
        angles = zeros(size(ephemerisTime,1),3);
        rates = zeros(size(ephemerisTime,1),3);
        velFlag = true;
end

%% Main

% Initialize loadFlag
loadFlag = [];
registerNumberFlag = [];

% Iterations for et inputs
for m = 1:size(ephemerisTime,1)
    
    % Register calculation
    [registerNumber,~,t,warnFlag] = ephRegister(ephemerisTime(m,1:end),'AU',ephConstants,...
        de,errorFlag,warnFlag);
    
    % Check if the register number has changed from a previous run
    [loadFlag,registerNumberFlag] = ephLoadFlag(loadFlag,registerNumberFlag,...
                                    registerNumber);
    
    % Calculation of angular position and rate
    if loadFlag
        bufLib = objLib.ephData(1:end,registerNumber);
    end
    pvec = ephInterp(bufLib,t,ephConstants.POINTERS(2,13),3,...
                     ephConstants.POINTERS(3,13),velFlag);
    % Output
    if velFlag
        angles(m,:)=pvec(:,1)';
        rates(m,:)=pvec(:,2)';
    else
        angles(m,:)=pvec(:,1)';
    end
end

if velFlag
    varargout{1} = angles;
    varargout{2} = rates;
else
    varargout{1} = angles;
end
%[EOF]
