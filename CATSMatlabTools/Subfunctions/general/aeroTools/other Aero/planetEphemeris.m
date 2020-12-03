function varargout = planetEphemeris(ephemerisTime,center,target,varargin)
%   PLANETEPHEMERIS Implement JPL development ephemeris for specified
%   bodies of the Solar System
%   [pos,vel] = PLANETEPHEMERIS(EPHEMERISTIME, CENTER, TARGET, MODEL, UNITS, ACTION) 
%   implements the mathematical representation of the planetary ephemeris
%   based on the Chebyshev coefficients calculated by NASA's Jet Propulsion
%   Laboratory (JPL). It calculates position and velocity of the planets
%   and celestial bodies that belong to the Solar System.
%   
%   
%   Inputs required by PLANETEPHEMERIS are:
%   EPHEMERISTIME  :Julian dates for which the positions are calculated. It
%                   can be a scalar or column vector with M elements.
%                   Alternatively, the date can be given as a 2 element
%                   vector or M by 2 matrix. The first element can be a
%                   fixed Julian date and the second element can be the
%                   elapsed time between the fixed Julian date and the
%                   ephemeris time.
%   CENTER         :String defining the reference body or point from which
%                   the position  and velocity of the TARGET planet is
%                   calculated. It can be: 'Sun', 'Mercury', 'Venus',
%                   'Earth', 'Moon', 'Mars', 'Jupiter','Saturn', 'Uranus',
%                   'Neptune','Pluto', 'SolarSystem' (for the Solar System
%                   barycenter) and 'EarthMoon' (for the Earth-Moon
%                   barycenter).
%   TARGET          :String defining the planet or point of interest for
%                   which the position and velocity is calculated. It can
%                   be: 'Sun', 'Mercury', 'Venus', 'Earth', 'Moon', 'Mars',
%                   'Jupiter', 'Saturn', 'Uranus', 'Neptune','Pluto',
%                   'SolarSystem' (for the Solar System barycenter) and
%                   'EarthMoon' (for the Earth-Moon barycenter).
%   MODEL           :String defining the Development Ephemeris (DE) by
%                   JPL. It can be: '405', '421', '423', '430' or '432t'.
%                   Ephemerides '405' and '421' are calculated with respect
%                   to the International Celestial Reference Frame version
%                   1.0, adopted in 1998. Ephemeris '423', '430' and '432t'
%                   are calculated with respect to the International
%                   Celestial Reference Frame version 2.0, adopted in 2010.
%                   If it is not defined, the default model is '405'.
%   UNITS           :String defining the output distance units. It can be
%                   'km' for position given in km and velocity in km/s or
%                   'AU' for position given in AU and velocity in AU/day.
%                   If it is not defined the units are set to 'km'.
%   ACTION          :String to determine action for out of range
%                   EPHEMERISTIME. Specify if out of range input invokes a
%                   'Warning', 'Error', or no action ('None'). The default
%                   is 'Error'.
%
%   Outputs calculated by PLANETEPHEMERIS are:
%   pos     :Numeric array describing the X, Y and Z coordinates of the
%            barycenter of the selected TARGET planet (or point of 
%            reference) with respect to the CENTER. Its size is M by 3,
%            with M rows, one for each ephemeris time. Its units are
%            defined by the UNITS input, and they can be either km or AU.
%   vel     :Numeric array describing the velocity along the X, Y and Z 
%            coordinates of the barycenter of the selected TARGET planet
%            (or point of reference) with respect to the CENTER. Its size
%            is M by 3, with M rows, one for each ephemeris time. Its units
%            are defined by the UNITS input, and they can be either km/s or
%            AU/day.
%
%   Examples:
%
%   Calculate the position of the Moon with respect to the Earth for
%   December 1st 1990 with DE405:
%       pos = planetEphemeris(juliandate(1990,12,1),'Earth','Moon')
%
%   Calculate the position and velocity for Saturn with respect to the 
%   Solar System barycenter for noon January 1st, 2000 using DE421 and AU 
%   units:
%       [pos,vel] = planetEphemeris([2451544.5 0.5],'SolarSystem','Saturn','421','AU')
%
%   Note: This function uses ephemeris data that can be obtained using the
%   aeroDataPackage command.
%
%   See also MOONLIBRATION, EARTHNUTATION.

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
narginchk(3,6);
nargoutchk(0,2);

% Valid inputs
validPlanets = {'Sun','Mercury','Venus','Earth','Moon','Mars','Jupiter',...
    'Saturn','Uranus','Neptune','Pluto','EarthMoon','SolarSystem'};
validDatabases = {'405','421','423','430','432t'};
validError = {'Error','Warning','None'};
validUnits = {'AU','km'};

% Planet validation
center = lower(validatestring(center,validPlanets));
target = lower(validatestring(target,validPlanets));

% Validate units and databases depending on number of inputs
switch nargin
    case 3
        de = '405';
        units = 'km';
        errorFlag = 'error';       
    case 4
        de = varargin{1};
        % Validate database
        if isnumeric(de)
            de = num2str(de);
        end
        de = lower(validatestring(de,validDatabases));
        units = 'km';
        errorFlag = 'error';
    case 5
        de = varargin{1};
        units = varargin{2};
        % Validate database
        if isnumeric(de)
            de = num2str(de);
        end
        de = lower(validatestring(de,validDatabases));
        % Validate units
        units = lower(validatestring(units,validUnits));
        errorFlag = 'error';
    case 6
        de = varargin{1};
        units = varargin{2};
        errorFlag = varargin{3};
        % Validate database
        if isnumeric(de)
            de = num2str(de);
        end
        de = lower(validatestring(de,validDatabases));
        % Validate units
        units = lower(validatestring(units,validUnits));
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
try	 load(['ephConstants' de])
 catch err	 
     if strcmp(err.identifier,'MATLAB:load:couldNotReadFile')	 
         error(message('aero:aeroephemerides:unavailableDatabase',...	 
                       '<a href="matlab:aeroDataPackage">aeroDataPackage</a>'));	 
     else	 
         rethrow(err)	 
     end	 
 end
% Load Objects
% Load center object
[objCent1,objCent2,calCent,indCent] = ephSelect(center,de);
% Initialize buffers
bufCent1 = [];
bufCent2 = [];
% Load target object
[objTarg1,objTarg2,calTarg,indTarg] = ephSelect(target,de);
% Initialize buffers
bufTarg1 = [];
bufTarg2 = [];
% Initialize warning flags
warnFlag = [0 0];

%% Output Initialization
switch nargout
    case {0,1} % Only position requested
        pos = zeros(size(ephemerisTime,1),3);
        velFlag = false;
    case 2 % Position and velocity requested
        pos = zeros(size(ephemerisTime,1),3);
        vel = zeros(size(ephemerisTime,1),3);
        velFlag = true;
end

%% Main
% Initialize loadFlag
loadFlag = [];
registerNumberFlag = [];

% Iterations for et inputs
for m = 1:size(ephemerisTime,1)
    
    % Register calculation
    [registerNumber,aufac,t,warnFlag] = ephRegister(ephemerisTime(m,1:end),units,...
        ephConstants,de,errorFlag,warnFlag);
    
    % Check that the target and center position are not the same, if so,
    % return zeros.
    if ~strcmp(target,center)
        % Check if the register number has changed from a previous run
        [loadFlag,registerNumberFlag] = ephLoadFlag(loadFlag,registerNumberFlag,...
            registerNumber);
        
        % Calculation of position and velocity
        if strcmpi(center,'Earth') && strcmpi(target,'Moon')
            if loadFlag
                bufTarg1 = objTarg1.ephData(1:end,registerNumber);
            end
            pvec = aufac*ephInterp(bufTarg1,t,ephConstants.POINTERS(2,10),3,...
                ephConstants.POINTERS(3,10),velFlag);
        elseif strcmpi(target,'Earth') && strcmpi(center,'Moon')
            if loadFlag
                bufCent1 = objCent1.ephData(1:end,registerNumber);
            end
            pvec = -aufac*ephInterp(bufCent1,t,ephConstants.POINTERS(2,10),3,...
                ephConstants.POINTERS(3,10),velFlag);
        else
            [pvtarg,bufTarg1,bufTarg2] = ephCalc(calTarg,objTarg1,objTarg2,bufTarg1,...
                bufTarg2,t,ephConstants,indTarg,velFlag,registerNumber,loadFlag);
            [pvcent,bufCent1,bufCent2] = ephCalc(calCent,objCent1,objCent2,bufCent1,...
                bufCent2,t,ephConstants,indCent,velFlag,registerNumber,loadFlag);
            pvec = aufac*(pvtarg-pvcent);
        end
        
        % Output
        if velFlag
            pos(m,:)=pvec(:,1)';
            vel(m,:)=pvec(:,2)';
        else
            pos(m,:)=pvec(:,1)';
        end
    end
end
if velFlag
    varargout{1} = pos;
    varargout{2} = vel;
else
    varargout{1} = pos;
end
%[EOF]
