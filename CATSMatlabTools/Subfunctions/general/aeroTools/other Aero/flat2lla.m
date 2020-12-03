function lla = flat2lla( p, ll0, psi0, href, varargin )
%  FLAT2LLA Estimate geodetic latitude, longitude, and altitude from flat
%    Earth position. 
%
%   LLA = FLAT2LLA ( P,  LL0,  PSI0, HREF ) estimates an M-by-3 array of
%   flat Earth  coordinates, P, to an M-by-3 array of geodetic coordinates
%   (latitude, longitude and altitude), LLA, with respect to a reference
%   location defined by initial latitude and longitude for the origin of
%   the estimation and the origin of the flat Earth coordinate system, LL0,
%   angular direction of the flat Earth x-axis (degrees clockwise from
%   north), PSI0, and reference height from the surface of the Earth to the
%   flat Earth frame w.r.t flat Earth frame, HREF, in meters.  P is in meters.
%   LL0 is in [degrees degrees].  LLA is in [degrees degrees meters].  The
%   default ellipsoid planet is WGS84.
%
%   LLA = FLAT2LLA ( P, LL0, PSI0, HREF, MODEL ) is an alternate method for
%   converting the coordinates, for a specific ellipsoid planet.  Currently
%   only 'WGS84' is supported for MODEL. 
%
%   LLA = FLAT2LLA (P, LL0, PSI0, HREF, F, RE ) is another alternate method
%   for converting the coordinates for a custom ellipsoid planet defined by
%   flattening, F, and the equatorial radius, RE in meters.
%
%   Limitations:
%
%   This estimation method assumes the flight path and bank angle are zero.
%
%   This estimation method assumes the flat Earth z-axis is normal to the
%   Earth at the initial geodetic latitude and longitude only. This method
%   has higher accuracy over small distances from the initial geodetic
%   latitude and longitude, and nearer to the equator. The longitude will
%   have higher accuracy the smaller the variations in latitude.
%   Additionally, longitude is singular at the poles. 
%
%   Examples:
%
%   Estimate latitude, longitude and altitude at a coordinate: 
%      lla = flat2lla( [ 4731 4511 120 ], [0 45], 5, -100)
%
%   Estimate latitude, longitude and altitude at multiple coordinates
%   specifying WGS84 ellipsoid model:
%      lla = flat2lla( [ 4731 4511 120; 0 5074 4498 ], [0 45], 5, -100, 'WGS84' )
%
%   Estimate latitude, longitude and altitude at multiple coordinates
%   specifying custom ellipsoid model:
%      f = 1/196.877360;
%      Re = 3397000;
%      lla = flat2lla( [ 4731 4511 120; 0 5074 4498 ], [0 45], 5, -100,  f, Re )
%
%   See also LLA2FLAT.

%   Copyright 2010-2014 The MathWorks, Inc.

%   References:
%   [1] Etkin, B., Dynamics of Atmospheric Flight, John Wiley & Sons, New
%   York, 1972.
%   [2] Stevens, B. L., and F. L. Lewis, Aircraft Control and Simulation,
%   Second Edition, John Wiley & Sons, New York, 2003. 

narginchk(4, 6);

if ~isnumeric( p ) || ~isreal( p ) || any(any(~isfinite( p )))
    error(message('aero:flat2lla:notNumericP'));
end

if (size( p, 2) ~= 3)
    error(message('aero:flat2lla:wrongDimensionP'));
end

if ~isnumeric( ll0 ) || ~isreal( ll0 ) || any(any(~isfinite( ll0 )))
    error(message('aero:flat2lla:notNumericLL0'));
end

if (length( ll0 ) ~= 2)
    error(message('aero:flat2lla:wrongDimensionLL0'));
end

if ~isnumeric( psi0 ) || ~isscalar( psi0 ) || ~isreal( psi0 ) || ~isfinite( psi0 )
    error(message('aero:flat2lla:notScalarPsi0'));
end

if ~isnumeric( href ) || ~isscalar( href ) || ~isreal( href ) || ~isfinite( href )
    error(message('aero:flat2lla:notScalarHRef'));
end

[f,R] = worldparams(nargin - 3, varargin);

% wrap latitude and longitude if needed
[~, ll0(1), ll0(2)] = wraplatitude( ll0(1), ll0(2), 'deg' );

% check and fix angle wrapping in longitude
[~, ll0(2)] = wraplongitude( ll0(2), 'deg', '180' );

rll0 = convang(ll0,'deg','rad');
rpsi0 = convang(psi0,'deg','rad');

dNorth = cos(rpsi0)*p(:,1) - sin(rpsi0)*p(:,2);
dEast  = sin(rpsi0)*p(:,1) + cos(rpsi0)*p(:,2);
lla(:,3) = -p(:,3) - href;

Rn = R/sqrt(1-(2*f-f*f)*sin(rll0(1))*sin(rll0(1)));
Rm = Rn*((1-(2*f-f*f))/(1-(2*f-f*f)*sin(rll0(1))*sin(rll0(1))));

dLat = dNorth.*atan2(1,Rm);
dLon = dEast.*atan2(1,Rn*cos(rll0(1)));

lla(:,1) = convang(dLat,'rad','deg') + ll0(1);
lla(:,2) = convang(dLon,'rad','deg') + ll0(2);

% wrap latitude and longitude if needed
[~, lla(:,1), lla(:,2)] = wraplatitude( lla(:,1)', lla(:,2)', 'deg' );

% check and fix angle wrapping in longitude
[~, lla(:,2)] = wraplongitude( lla(:,2), 'deg', '180' );
