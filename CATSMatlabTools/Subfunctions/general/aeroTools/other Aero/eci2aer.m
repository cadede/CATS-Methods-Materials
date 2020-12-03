function aer = eci2aer(POSITION,UTC,LLA0,varargin)
%   ECI2AER Convert Earth-centered Inertial (ECI) to local azimuth,
%   elevation, and slant range coordinates.
%   AER = ECI2AER( POSITION,UTC,LLA0,REDUCTION,DELTAAT,DELTAUT1,POLARMOTION,'ADDPARAMNAME',ADDPARAMVALUE )
%   converts a set of ECI Cartesian coordinates, POSITION, to local azimuth,
%   elevation and slant range coordinates.
%
%   Inputs arguments for ECI2AER are:
%   POSITION:       M-by-3 array of ECI coordinates in meters.
%   UTC:            Array of Universal Coordinated Time (UTC) in year,
%                   month, day, hour, minutes, and seconds for which the
%                   function calculates the coordinate conversion. Define
%                   the array as one of the following: Array with 1 row and
%                   6 columns, or M-by-6 array for M transformation
%                   matrices, one for each UTC date. Values for year,
%                   month, day, hour, and minutes should be whole numbers.
%   LLA0:           M-by-3 array with the geodetic coordinates of the local
%                   reference with latitude, longitude and ellipsoidal
%                   altitude in degrees, degrees and meters respectively.
%   REDUCTION:      String indicating the reduction process through which
%                   the function calculates the coordinate conversion.
%                   It can either be IAU-76/FK5 (which uses the IAU 1976
%                   Precession Model and the IAU 1980 Theory of Nutation,
%                   which is no longer current but some programs still use
%                   this reduction) or IAU-2000/2006 (which uses the P03
%                   precession model). The reduction method you select
%                   determines the ADDPARAMNAME parameter pair
%                   characteristics. The IAU-76/FK5 method returns a
%                   coordinate conversion that is not orthogonal due to the
%                   polar motion approximation. The default value is
%                   IAU-2000/2006.
%   DELTAAT:        Difference, in seconds, between the International Atomic
%                   Time (TAI) and UTC. It can be defined as either a
%                   scalar or a one-dimensional array with M elements (if M
%                   UTC dates are defined), for equal number ECI
%                   coordinates. The default value is an M-by-1 null array.
%   DELTAUT1:       Difference, in seconds, between UTC and Universal Time
%                   (UT1). It can be defined as either a scalar or a one-
%                   dimensional array with M elements (if M UTC dates 
%                   defined) for equal number ECI coordinates. The default
%                   value is an M-by-1 null array.
%   POLARMOTION:    Polar displacement due to the motion of the Earth's
%                   crust, in radians, along the x- and y-axis. It can be
%                   defined as either a 1-by-2 array or an M-by-2 (if M UTC
%                   dates defined) for equal number of ECI coordinates. The
%                   default value is an M-by-2 null array.
%
%   Input parameter Name/Value pairs for 'ADDPARAMNAME' and ADDPARAMVALUE
%   are: 
%   'DNUTATION':    (IAU-76/FK5 reduction only) M-by-2 array for the
%                   adjustment in radians to the longitude (dDeltaPsi) and
%                   obliquity (dDeltaEpsilon). The default value is an
%                   M-by-2 null array. 
%   'DCIP':         (IAU-2000/2006 reduction only) M by 2 array for the
%                   adjustment in radians to the location of the Celestial
%                   Intermediate Pole (CIP) along the x (dDeltaX) and y
%                   (dDeltaY) axis. The default value is an M-by-2 null
%                   array.
%   'FLATTENING':   Earth's flattening. The default value is the one
%                   defined by WGS84 (1/298.257223563).
%   'RE':           Earth's equatorial radius. The default value is the one
%                   defined by WGS84 (6378137 meters).
%   'PSI0':         Angular direction of the local reference system (degrees
%                   clockwise from north). The default value is an M-by-1
%                   null array.
%
%   For historical values for DNUTATION and DCIP, see the International
%   Earth Rotation and Reference Systems Service (IERS) website
%   (http://www.iers.org) under the 'Earth Orientation Data' product. 
%
%   Output calculated by ECI2AER is:
%   AER:    M-by-3 array with the local reference coordinates azimuth,
%           elevation, and slant range in degrees, degrees, and meters,
%           respectively. Azimuth is defined as the angle measured
%           clockwise from true north and varies between 0 and 360 deg.
%           Elevation is defined as the angle between a plane perpendicular
%           to the ellipsoid's surface at LL0 and the line that goes from
%           the local reference to the object's position and varies between
%           -90 and 90 degrees. Slant range is defined as the straight line
%           distance between the local reference and the object.
%
%   This method has higher accuracy over small distances from the local
%   geodetic reference frame (LLA0) in latitude and longitude, and nearer to
%   the equator.
%
%   Example:
%   Estimate the position of the azimuth, elevation and range for an object
%   with Earth Centered Inertial position 1e08*[-3.8454 -0.5099 -0.3255]
%   meters for the date 1969/7/20 21:17:40 UTC at 28.4 deg North, 80.5 deg
%   West and 2.7 meters altitude.
%
%   AER = eci2aer(1e08*[-3.8454,-0.5099,-0.3255],[1969,7,20,21,17,40], ...
%                 [28.4,-80.5,2.7])
% 
%   See also LLA2ECI, ECI2LLA, DCMECI2ECEF, FLAT2LLA, LLA2FLAT, DELTAUT1, 
%   DELTACIP, POLARMOTION.

%   Copyright 2014-2018 The MathWorks, Inc.
%% Validate i/o

% Validate outputs
nargoutchk(0,1)

% Validate date
validateattributes(UTC,{'numeric'},{'ncols',6,'real','finite','nonnan'})

% Validate latitude longitude altitude
validateattributes(POSITION,{'numeric'},{'ncols',3,'real','finite','nonnan'})

% Assign date vectors
year = UTC(:,1);
month = UTC(:,2);
day = UTC(:,3);
hour = UTC(:,4);
min = UTC(:,5);
sec = UTC(:,6);

% Validate vectors
if any(year<1)
    error(message('aero:dcmeci2ecef:invalidYear'));
end
if any(month<1) || any(month>12)
    error(message('aero:dcmeci2ecef:invalidMonth'));
end
if any(day<1) || any(day>31)
    error(message('aero:dcmeci2ecef:invalidDay'));
end
if any(hour<0) || any(hour>24)
    error(message('aero:dcmeci2ecef:invalidHour'));
end
if any(min<0) || any(min>60)
    error(message('aero:dcmeci2ecef:invalidMin'));
end
if any(sec<0) || any(sec>60)
    error(message('aero:dcmeci2ecef:invalidSec'));
end
len = length(year);

% Parse and validate the inputs
ob = inputParser;
validReduction = {'IAU-2000/2006','IAU-76/FK5'};
addRequired(ob,'POSITION',@(x) validateattributes(x,{'numeric'},{'ncols',3,'real',...
    'finite','nonnan','size',[len,3]}));
addRequired(ob,'UTC',@(x) validateattributes(x,{'numeric'},{'ncols',6,'real',...
    'finite','nonnan'}));
addRequired(ob,'LLA0',@(x) validateattributes(x,{'numeric'},{'ncols',3,'real',...
    'finite','nonnan','size',[len,3]}));
addOptional(ob,'reduction','iau-2000/2006',@(x) validateReduction(x,validReduction));
addOptional(ob,'deltaAT',zeros(len,1),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,1]}));
addOptional(ob,'deltaUT1',zeros(len,1),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,1]}));
addOptional(ob,'polarMotion',zeros(len,2),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,2]}));
addParameter(ob,'dNutation',zeros(len,2),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,2]}));
addParameter(ob,'dCIP',zeros(len,2),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,2]}));
addParameter(ob,'flattening',1/298.257223563,@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[1 1]}));
addParameter(ob,'RE',6378137,@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[1 1]}));
addOptional(ob,'PSI0',zeros(len,1),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,1]}));

% Parse input object
parse(ob,POSITION,UTC,LLA0,varargin{:});

% Validate reduction
reduction = ob.Results.reduction;
reduction = lower(validatestring(reduction,validReduction));

%Validate that the additional parameter matches the reduction method
if any(strcmp(ob.UsingDefaults,'dNutation')) && ~any(strcmp(ob.UsingDefaults,'dCIP')) &&...
        ~strcmp(reduction,'iau-2000/2006')
    error(message('aero:dcmeci2ecef:invalidDNutation'));
elseif ~any(strcmp(ob.UsingDefaults,'dNutation')) && any(strcmp(ob.UsingDefaults,'dCIP')) && ...
        ~strcmp(reduction,'iau-76/fk5')
    error(message('aero:dcmeci2ecef:invalidDCIP'));
end

%Validate that both reduction methods are not defined (just one should be
%defined)
if ~any(strcmp(ob.UsingDefaults,'dNutation')) && ~any(strcmp(ob.UsingDefaults,'dCIP'))
    error(message('aero:dcmeci2ecef:invalidDNutationDCIP'))
end

%% Calculate DCM ECI to ECEF
switch reduction
    case 'iau-76/fk5'
        dcm = dcmeci2ecef(ob.Results.reduction,ob.Results.UTC,ob.Results.deltaAT,...
            ob.Results.deltaUT1,ob.Results.polarMotion,'dNutation',ob.Results.dNutation);
    case 'iau-2000/2006'
        dcm = dcmeci2ecef(ob.Results.reduction,ob.Results.UTC,ob.Results.deltaAT,...
            ob.Results.deltaUT1,ob.Results.polarMotion,'dCIP',ob.Results.dCIP);
end
%% Calculate position in ECEF coordinates
tmp = arrayfun(@(k) (dcm(:,:,k)*POSITION(k,:)'),1:len,'UniformOutput',false);
ecefObject = cell2mat(tmp)';

%% Calculate position in NED coordinates
% Wrap latitude and longitude for the local geodetic reference system
LLA0 = ob.Results.LLA0;
[~, LLA0(:,1), LLA0(:,2)] = wraplatitude( ob.Results.LLA0(:,1), ob.Results.LLA0(:,2),'deg');
[~, LLA0(:,2)] = wraplongitude( LLA0(:,2), 'deg', '180' );
% Use local function for NED transform
tmp2 = arrayfun(@(k) (nedCalc(ecefObject(k,:),LLA0(k,:),ob.Results.flattening,ob.Results.RE)),1:len,'UniformOutput',false);
nedObject = cell2mat(tmp2');

%% Calculate Azimuth, Elevation, Range
hypotxy = hypot(nedObject(:,1),nedObject(:,2));
r = hypot(hypotxy,-nedObject(:,3));
elev = atan2d(-nedObject(:,3),hypotxy);
az = atan2d(nedObject(:,2),nedObject(:,1))-ob.Results.PSI0;
az = mod(az,360);
aer = [az elev r];

end

function validateReduction(reduction,validReduction)
validatestring(reduction,validReduction);
end

function ned = nedCalc(posEcef,LLA0,f,Re)
% This helper function calculates the NED coordinates for the object given
% in ECEF coordinates given a local geodetic position LLA0. It requires the
% flattening and equatorial radius. It calculates the ENU coordinates and
% finally rotates it for NED.

refPosEcef = lla2ecef(LLA0,f,Re);
dPos = posEcef - refPosEcef;

% ENU position
enu = [ -sind(LLA0(2))                  cosd(LLA0(2))                   0;...
        -cosd(LLA0(2))*sind(LLA0(1))    -sind(LLA0(1))*sind(LLA0(2))    cosd(LLA0(1)); ...
        cosd(LLA0(1))*cosd(LLA0(2))     cosd(LLA0(1))*sind(LLA0(2))     sind(LLA0(1))]*...
        dPos';
    
% NED position
ned = [enu(2) enu(1) -enu(3)];

end
