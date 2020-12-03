function p = lla2ecef( lla, varargin )
%  LLA2ECEF Convert geodetic coordinates to Earth-centered Earth-fixed
%  (ECEF)  coordinates.
%   P = LLA2ECEF( LLA ) converts an M-by-3 array of geodetic coordinates
%   (latitude, longitude and altitude), LLA, to an M-by-3 array of ECEF
%   coordinates, P.  LLA is in [degrees degrees meters].  P is in meters.
%   The default ellipsoid planet is WGS84. 
%
%   P = LLA2ECEF( LLA, MODEL ) is an alternate method for converting
%   the coordinates, for a specific ellipsoid planet.  Currently only 'WGS84' is
%   supported for MODEL.
%
%   P = LLA2ECEF( LLA, F, RE ) is another alternate method for
%   converting the coordinates for a custom ellipsoid planet defined by
%   flattening, F, and the equatorial radius, RE in meters. 
%
%   Examples:
%
%   Determine ECEF coordinates at a latitude, longitude and altitude: 
%      p = lla2ecef( [ 0 45 1000 ] )
%
%   Determine ECEF coordinates at multiple latitude, longitude and altitude
%   specifying WGS84 ellipsoid model:
%      p = lla2ecef( [  0 45 1000; 45 90 2000 ], 'WGS84' )
%
%   Determine ECEF coordinates at multiple latitude, longitude and altitude
%   specifying custom ellipsoid model:
%      f = 1/196.877360;
%      Re = 3397000;
%      p = lla2ecef( [  0 45 1000; 45 90 2000 ],  f, Re )
%
%   See also ECEF2LLA, GEOC2GEOD, GEOD2GEOC, LLA2ECI, ECI2LLA.

%   Copyright 2000-2013 The MathWorks, Inc. 

narginchk(1, 3);

if ~isnumeric( lla )
    error(message('aero:lla2ecef:notNumeric'));
end

if (size( lla, 2) ~= 3)
    error(message('aero:lla2ecef:wrongDimension'));
end

[f,R] = worldparams(nargin, varargin);

% Determine radial distance from polar axis (rho) and signed distance
% from the equator in a spheroid-centric (ECEF) cylindrical coordinate
% system.
[rho,z] = map.geodesy.internal.geodetic2cylindrical(lla(:,1),lla(:,3),R,f,true);
% Convert from cylindrical to Cartesian coordinates.
x = rho.*cosd(lla(:,2));
y = rho.*sind(lla(:,2));
p = [x y z];
