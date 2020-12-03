function lla = ecef2lla( p, varargin )
%  ECEF2LLA Convert Earth-centered Earth-fixed (ECEF) coordinates to
%  geodetic coordinates. 
%   LLA = ECEF2LLA( P ) converts the M-by-3 array of ECEF coordinates, P, to
%   an M-by-3 array of geodetic coordinates (latitude, longitude and
%   altitude), LLA.  LLA is in [degrees degrees meters].  P is in meters.
%   The default ellipsoid planet is WGS84. 
%
%   LLA = ECEF2LLA( P, MODEL ) is an alternate method for converting
%   the coordinates, for a specific ellipsoid planet.  Currently only 'WGS84' is
%   supported for MODEL.
%
%   LLA = ECEF2LLA( P, F, RE ) is another alternate method for
%   converting the coordinates for a custom ellipsoid planet defined by
%   flattening, F, and the equatorial radius, RE in meters. 
%
%   Examples:
%
%   Determine latitude, longitude and altitude at a coordinate: 
%      lla = ecef2lla( [ 4510731 4510731 0 ] )
%
%   Determine latitude, longitude and altitude at multiple coordinates
%   specifying WGS84 ellipsoid model:
%      lla = ecef2lla( [ 4510731 4510731 0; 0 4507609 4498719 ], 'WGS84' )
%
%   Determine latitude, longitude and altitude at multiple coordinates
%   specifying custom ellipsoid model:
%      f = 1/196.877360;
%      Re = 3397000;
%      lla = ecef2lla( [ 4510731 4510731 0; 0 4507609 4498719 ],  f, Re )
%
%   See also GEOC2GEOD, GEOD2GEOC, LLA2ECEF, LLA2ECI, ECI2LLA.

%   Copyright 2000-2013 The MathWorks, Inc.

narginchk(1, 3);

if ~isnumeric( p )
    error(message('aero:ecef2lla:notNumeric1'));
end

if (size( p, 2) ~= 3)
    error(message('aero:ecef2lla:wrongDimension'));
end

[f,R] = worldparams(nargin, varargin);

% Determine longitude
lambda = atan2d(p(:,2),p(:,1));
% Determine radial distance from polar axis
rho = hypot(p(:,1),p(:,2));
% Determine geodetic latitude and ellipsoidal height
[phi,h] = map.geodesy.internal.cylindrical2geodetic(rho,p(:,3),R,f,true);

lla = [phi lambda h];
