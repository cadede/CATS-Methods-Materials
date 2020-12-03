function gc  = geod2geoc( gd, h, varargin )
%  GEOD2GEOC Convert geodetic latitude to geocentric latitude.
%   GC = GEOD2GEOC( GD, H ) converts an array of M geodetic latitudes,
%   GD, and an array of mean sea level altitude, H, into an array of M
%   geocentric latitudes, GC. Both GC and GD are in degrees. H is in meters.    
%
%   GC = GEOD2GEOC( GD, H, MODEL ) is an alternate method for converting
%   from geodetic latitude, for a specific ellipsoid planet.  Currently only 'WGS84' is
%   supported for MODEL.
%
%   GC = GEOD2GEOC( GD, H, F, RE ) is another alternate method for
%   converting from geodetic latitude for a custom ellipsoid planet defined by
%   flattening, F, and the equatorial radius, RE in meters. 
%
%   Limitation:
%
%   This implementation generates a geocentric latitude that lies between +/-
%   90 degrees.
%
%   Examples:
%
%   Determine geocentric latitude at a geodetic latitude, and altitude: 
%      gc = geod2geoc( 45, 1000 )
%
%   Determine geocentric latitude at multiple geodetic latitude, and
%   altitude specifying WGS84 ellipsoid model:
%      gc = geod2geoc( [ 0 45 90 ], [1000 0 2000], 'WGS84' )
%
%   Determine geocentric latitude at multiple geodetic latitude, and
%   altitude specifying custom ellipsoid model:
%      f = 1/196.877360;
%      Re = 3397000;
%      gc = geod2geoc([ 0 45 90 ], 2000,  f, Re )
%
%   See also ECEF2LLA, GEOC2GEOD, LLA2ECEF, LLA2ECI, ECI2LLA.

%   Copyright 2000-2018 The MathWorks, Inc.

%   References:
%   [1] Stevens, B. L., and F. L. Lewis, Aircraft Control and Simulation, John
%   Wiley & Sons, New York, NY, 1992.  

narginchk(2, 4);

if ~isnumeric( gd )
    error(message('aero:geod2geoc:notNumericPosition'));
end

if ~isnumeric( h )
    error(message('aero:geod2geoc:notNumericAltitude'));
end

[f,R] = worldparams(nargin - 1, varargin);

[~, gd, ~] = wraplatitude( gd, zeros(size(gd)), 'deg');

try
    % Determine radial distance from polar axis (rho) and signed distance
    % from the equator in a spheroid-centric (ECEF) cylindrical coordinate
    % system.
    [rho,z] = map.geodesy.internal.geodetic2cylindrical(gd,h,R,f,true);
    % Determine geocentric latitude.
    gc = atan2d(z,rho);
catch wrongdim
    newExc = MException('aero:geod2geoc:wrongDimension', ...
        getString(message('aero:geod2geoc:wrongDimension')));
    newExc = newExc.addCause(wrongdim);
    throw(newExc);
end