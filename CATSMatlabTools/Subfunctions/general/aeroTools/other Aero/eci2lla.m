function lla = eci2lla(P,UTC,varargin)
%   ECI2LLA Convert Earth-centered Inertial (ECI) to geodetic coordinates.
%   LLA = ECI2LLA( P,UTC,REDUCTION,DELTAAT,DELTAUT1,POLARMOTION,'ADDPARAMNAME',ADDPARAMVALUE )
%   converts a set of ECI coordinates, P, to geodetic coordinates
%   (latitude, longitude and altitude), LLA.
%
%   Inputs arguments for ECI2LLA are:
%   P:              M by 3 array of ECI coordinates in meters.
%   UTC:            Array of Universal Coordinated Time (UTC) in year,
%                   month, day, hour, minutes and seconds for which the
%                   function calculates the coordinate conversion. Define
%                   the array as one of the following: Array with 1 row and
%                   6 columns or M by 6 array for M transformation
%                   matrices, one for each UTC date. Values for year,
%                   month, day, hour and minutes should be whole numbers.
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
%   DELTAAT:        Difference in seconds between the International Atomic
%                   Time (TAI) and UTC. It can be defined as either a
%                   scalar or a one dimensional array with M elements (if M
%                   UTC dates defined) for equal number of ECI coordinates.
%                   The default value is an M by 1 null array.
%   DELTAUT1:       Difference, in seconds, between UTC and Universal Time
%                   (UT1). It can be defined as either a scalar or a one
%                   dimensional array with M elements (if M UTC dates 
%                   defined) for equal number ECI coordinates. The default
%                   value is an M by 1 null array.
%   POLARMOTION:    Polar displacement due to the motion of the Earth's
%                   crust, in radians, along the x and y axis. It can be
%                   defined as either a 1 by 2 array or an M by 2 (if M UTC
%                   dates defined) for equal number of ECI coordinates. The
%                   default value is an M by 2 null array.
%
%   Input parameter Name/Value pairs for 'ADDPARAMNAME' and ADDPARAMVALUE
%   are: 
%   'DNUTATION':    (IAU-76/FK5 reduction only) M by 2 array for the
%                   adjustment in radians to the longitude (dDeltaPsi) and
%                   obliquity (dDeltaEpsilon). The default value is an M by
%                   2 null array.
%   'DCIP':         (IAU-2000/2006 reduction only) M by 2 array for the
%                   adjustment in radians to the location of the Celestial
%                   Intermediate Pole (CIP) along the x (dDeltaX) and y
%                   (dDeltaY) axis. The default value is an M by 2 null
%                   array.
%   'FLATTENING':   Earth's flattening. The default value is the one
%                   defined by WGS84 (1/298.257223563).
%   'RE':           Earth's equatorial radius. The default value is the one
%                   defined by WGS84 (6378137 meters).
%
%   For historical values for DNUTATION and DCIP, see the International
%   Earth Rotation and Reference Systems Service (IERS) website
%   (http://www.iers.org) under the 'Earth Orientation Data' product. 
%
%   Output calculated by ECI2LLA is:
%   LLA:    M by 3 array with the geodetic coordinates latitude,
%           longitude and altitude in degrees, degrees and meters
%           respectively.
%
%   Examples:
%   Calculate the latitude, longitude and altitude for the ECI coordinates
%   [-6.07 -1.28 0.66]*1e6 at 01/17/2010 10:20:36 UTC.
%
%   LLA = eci2lla([-6.07 -1.28 0.66]*1e6,[2010 1 17 10 20 36])
%
%   Calculate the latitude, longitude and altitude for the ECI coordinates
%   [-1.1 3.2 -4.9]*1e4 at 01/12/2000 4:52:12.4 UTC, with a difference
%   of 32 seconds between TAI and UTC, and 0.234 seconds between UTC and
%   UT1; use the IAU-76/FK5 reduction, polar motion [-0.0682e-5 0.1616e-5]
%   radians and nutation angles [-0.2530e-6 -0.0188e-6] for an ellipsoid
%   with a flattening of 1/290 and an equatorial radius of 60000 meters.
%
%   LLA = eci2lla([-1.1 3.2 -4.9]*1e4,[2000 1 12 4 52 12.4],'IAU-76/FK5',32,...
%           0.234,[-0.0682e-5 0.1616e-5],'dNutation',[-0.2530e-6 -0.0188e-6],...
%           'flattening',1/290,'RE',60000)
%
%   See also LLA2ECI, DCMECI2ECEF, LLA2ECEF, ECEF2LLA, GEOC2GEOD, GEOD2GEOC,
%   DELTAUT1, DELTACIP, POLARMOTION.

%   Copyright 2013-2018 The MathWorks, Inc.
%% Validate i/o

% Validate outputs
nargoutchk(0,1)

% Validate date
validateattributes(UTC,{'numeric'},{'ncols',6,'real','finite','nonnan'})

% Validate latitude longitude altitude
validateattributes(P,{'numeric'},{'ncols',3,'real','finite','nonnan'})

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
p = inputParser;
validReduction = {'IAU-2000/2006','IAU-76/FK5'};
addRequired(p,'P',@(x) validateattributes(x,{'numeric'},{'ncols',3,'real',...
    'finite','nonnan','size',[len,3]}));
addRequired(p,'UTC',@(x) validateattributes(x,{'numeric'},{'ncols',6,'real',...
    'finite','nonnan'}));
addOptional(p,'reduction','iau-2000/2006',@(x) validateReduction(x,validReduction));
addOptional(p,'deltaAT',zeros(len,1),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,1]}));
addOptional(p,'deltaUT1',zeros(len,1),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,1]}));
addOptional(p,'polarMotion',zeros(len,2),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,2]}));
addParameter(p,'dNutation',zeros(len,2),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,2]}));
addParameter(p,'dCIP',zeros(len,2),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,2]}));
addParameter(p,'flattening',1/298.257223563,@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[1 1]}));
addParameter(p,'RE',6378137,@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[1 1]}));

% Parse input object
parse(p,P,UTC,varargin{:});

% Validate reduction
reduction = p.Results.reduction;
reduction = lower(validatestring(reduction,validReduction));

%Validate that the additional parameter matches the reduction method
if any(strcmp(p.UsingDefaults,'dNutation')) && ~any(strcmp(p.UsingDefaults,'dCIP')) &&...
        ~strcmp(reduction,'iau-2000/2006')
    error(message('aero:dcmeci2ecef:invalidDNutation'));
elseif ~any(strcmp(p.UsingDefaults,'dNutation')) && any(strcmp(p.UsingDefaults,'dCIP')) && ...
        ~strcmp(reduction,'iau-76/fk5')
    error(message('aero:dcmeci2ecef:invalidDCIP'));
end

%Validate that both reduction methods are not defined (just one should be
%defined)
if ~any(strcmp(p.UsingDefaults,'dNutation')) && ~any(strcmp(p.UsingDefaults,'dCIP'))
    error(message('aero:dcmeci2ecef:invalidDNutationDCIP'))
end

%% Calculate DCM ECI to ECEF
switch reduction
    case 'iau-76/fk5'
        dcm = dcmeci2ecef(p.Results.reduction,p.Results.UTC,p.Results.deltaAT,...
            p.Results.deltaUT1,p.Results.polarMotion,'dNutation',p.Results.dNutation);
    case 'iau-2000/2006'
        dcm = dcmeci2ecef(p.Results.reduction,p.Results.UTC,p.Results.deltaAT,...
            p.Results.deltaUT1,p.Results.polarMotion,'dCIP',p.Results.dCIP);
end
%% Calculate position in ECEF coordinates
tmp = arrayfun(@(k) (dcm(:,:,k)*P(k,:)'),1:len,'UniformOutput',false);
ecef = cell2mat(tmp)';
%% Calculate LLA coordinates
lla = ecef2lla(ecef,p.Results.flattening,p.Results.RE);

function validateReduction(reduction,validReduction)
validatestring(reduction,validReduction);
