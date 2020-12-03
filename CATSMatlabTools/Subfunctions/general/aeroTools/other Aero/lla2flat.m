function p = lla2flat( lla, ll0, psi0, href, varargin )
%  LLA2FLAT Estimate flat Earth position from geodetic latitude, longitude,
%   and altitude. 
%   P = LLA2FLAT( LLA,  LL0,  PSI0, HREF ) estimates an M-by-3 array of
%   geodetic coordinates (latitude, longitude and altitude), LLA, to an
%   M-by-3 array of flat Earth  coordinates, P with respect to a reference
%   location defined by initial latitude and longitude for the origin of
%   the estimation and the origin of the flat Earth coordinate system, LL0,
%   angular direction of the flat Earth x-axis (degrees clockwise from
%   north), PSI0, and reference height from the surface of the Earth to the
%   flat Earth frame w.r.t flat Earth frame, HREF, in meters.  LLA is in
%   [degrees degrees meters].  LL0 is in [degrees degrees].  P is in
%   meters.   The default ellipsoid planet is WGS84.
%
%   P = LLA2FLAT( LLA, LL0, PSI0, HREF, MODEL ) is an alternate method for
%   converting the coordinates, for a specific ellipsoid planet.  Currently
%   only 'WGS84' is supported for MODEL. 
%
%   P = LLA2FLAT( LLA, LL0, PSI0, HREF, F, RE ) is another alternate method
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
%   Estimate coordinates at a latitude, longitude and altitude: 
%      p = lla2flat( [ 0.1 44.95 1000 ], [0 45], 5, -100 )
%
%   Estimate coordinates at multiple latitude, longitude and altitude
%   specifying WGS84 ellipsoid model:
%      p = lla2flat( [ 0.1 44.95 1000; -0.05 45.3 2000 ], [0 45], 5, -100, 'WGS84' )
%
%   Estimate coordinates at multiple latitude, longitude and altitude
%   specifying custom ellipsoid model:
%      f = 1/196.877360;
%      Re = 3397000;
%      p = lla2flat( [ 0.1 44.95 1000; -0.05 45.3 2000 ], [0 45], 5, -100,  f, Re )
%
%   See also FLAT2LLA.

%   Copyright 2010-2014 The MathWorks, Inc.

%   References:
%   [1] Etkin, B., Dynamics of Atmospheric Flight, John Wiley & Sons, New
%   York, 1972.
%   [2] Stevens, B. L., and F. L. Lewis, Aircraft Control and Simulation,
%   Second Edition, John Wiley & Sons, New York, 2003. 

narginchk(4, 6);

if ~isnumeric( lla ) || ~isreal( lla ) || any(any(~isfinite( lla )))
    error(message('aero:lla2flat:notNumericLLA'));
end

if (size( lla, 2) ~= 3)
    error(message('aero:lla2flat:wrongDimensionLLA'));
end

if ~isnumeric( ll0 ) || ~isreal( ll0 ) || any(any(~isfinite( ll0 )))
    error(message('aero:lla2flat:notNumericLL0'));
end

if (length( ll0 ) ~= 2)
    error(message('aero:lla2flat:wrongDimensionLL0'));
end

if ~isnumeric( psi0 ) || ~isscalar( psi0 ) || ~isreal( psi0 ) || ~isfinite( psi0 )
    error(message('aero:lla2flat:notScalarPsi0'));
end

if ~isnumeric( href ) || ~isscalar( href ) || ~isreal( href ) || ~isfinite( href )
    error(message('aero:lla2flat:notScalarHRef'));
end

[f,R] = worldparams(nargin - 3, varargin);

% wrap latitude and longitude if needed
[~, lla(:,1), lla(:,2)] = wraplatitude( lla(:,1), lla(:,2), 'deg' );
[~, ll0(1), ll0(2)] = wraplatitude( ll0(1), ll0(2), 'deg' );

% check and fix angle wrapping in longitude
[~, lla(:,2)] = wraplongitude( lla(:,2), 'deg', '180' );
[~, ll0(2)] = wraplongitude( ll0(2), 'deg', '180' );

dLat = lla(:,1) - ll0(1);
dLon = lla(:,2) - ll0(2);

% wrap latitude and longitude if needed
[~, dLat, dLon] = wraplatitude( dLat', dLon', 'deg' );

% check and fix angle wrapping in longitude
[~, dLon] = wraplongitude( dLon, 'deg', '180' );

rll0 = convang(ll0,'deg','rad');
rpsi0 = convang(psi0,'deg','rad');

Rn = R/sqrt(1-(2*f-f*f)*sin(rll0(1))*sin(rll0(1)));
Rm = Rn*((1-(2*f-f*f))/(1-(2*f-f*f)*sin(rll0(1))*sin(rll0(1))));

dNorth = convang(dLat,'deg','rad')./atan2(1,Rm);
dEast = convang(dLon,'deg','rad')./atan2(1,Rn*cos(rll0(1)));

p(:,1) = cos(rpsi0)*dNorth + sin(rpsi0)*dEast;
p(:,2) = -sin(rpsi0)*dNorth + cos(rpsi0)*dEast;
p(:,3) = -lla(:,3) - href;
