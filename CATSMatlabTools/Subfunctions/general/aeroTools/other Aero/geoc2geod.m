function gd = geoc2geod( gc , r, varargin )
%  GEOC2GEOD Convert geocentric latitude to geodetic latitude.
%   GD = GEOC2GEOD( GC, R ) converts an array of M geocentric latitudes,
%   GC, and an array of radii from the center of the planet, R, into an
%   array of M geodetic latitudes, GD. Both GC and GD are in degrees. R is in meters.   
%
%   GD = GEOC2GEOD( GC, R, MODEL ) is an alternate method for converting
%   from geocentric latitude, for a specific ellipsoid planet.  Currently only 'WGS84' is
%   supported for MODEL.
%
%   GD = GEOC2GEOD( GC, R, F, RE ) is another alternate method for
%   converting from geocentric latitude for a custom ellipsoid planet defined by
%   flattening, F, and the equatorial radius, RE in meters. 
%
%   Geometric relationships are used to calculate the geodetic latitude in
%   this non-iterative method.
%
%   Limitation:
%
%   This implementation generates a geodetic latitude that lies between +/-
%   90 degrees. 
%
%   Examples:
%
%   Determine geodetic latitude at a geocentric latitude, and radius: 
%      gd = geoc2geod( 45, 6379136 )
%
%   Determine geodetic latitude at multiple geocentric latitude, and radius
%   specifying WGS84 ellipsoid model:
%      gd = geoc2geod( [ 0 45 90 ], 6379136, 'WGS84' )
%
%   Determine geodetic latitude at multiple geocentric latitude, and radius
%   specifying custom ellipsoid model:
%      f = 1/196.877360;
%      Re = 3397000;
%      gd = geoc2geod([ 0 45 90 ], 6379136,  f, Re )
%
%   See also ECEF2LLA, GEOD2GEOC, LLA2ECEF, LLA2ECI, ECI2LLA.

%   Copyright 2000-2013 The MathWorks, Inc.

%   References:
%   [1] Jackson, E. B., Manual for a Workstation-based Generic Flight
%   Simulation Program (LaRCsim) Version 1.4, NASA TM 110164, April, 1995. 
%   [2] Hedgley, D. R., Jr., "An Exact Transformation from Geocentric to
%   Geodetic Coordinates for Nonzero Altitudes," NASA TR R-458, March, 1976. 
%   [3] Clynch, J. R., "Radius of the Earth - Radii Used in Geodesy," Naval
%   Postgraduate School, 2002,
%   http://www.oc.nps.navy.mil/oc2902w/geodesy/radiigeo.pdf. 
%   [4] Stevens, B. L., and F. L. Lewis, Aircraft Control and Simulation, John
%   Wiley & Sons, New York, NY, 1992.  
%   [5] Edwards, C. H., and D. E. Penny, Calculus and Analytical Geometry 2nd
%   Edition, Prentice-Hall, Englewood Cliffs, NJ, 1986. 

narginchk(2, 4);

if ~isnumeric( gc )
    error(message('aero:geoc2geod:notNumeric1'));
end

if ~isnumeric( r )
    error(message('aero:geoc2geod:notNumeric2'));
end

[~, gc, ~] = wraplatitude( gc, zeros(size(gc)), 'deg');

[f,R] = worldparams(nargin - 1, varargin);

lambda = convang( gc ,'deg','rad');

% By rearranging the equation for an ellipse determine the horizontal coordinate
xa = (1 - f ).* R ./ sqrt(tan(lambda).*tan(lambda) + ( 1 - f ).*( 1 - f ));

% Define vertical coordinate in terms of horizontal coordinate
% ya = sqrt(R.*R - xa.*xa).*( 1 - f );

% Use relationship between geocentric latitude at the planet's surface and
% geodetic latitude and tan(gc) = ya/xa; 
mua = atan2(sqrt(R.*R - xa.*xa),( 1 - f ).*xa);

for k = 1:length(lambda)
    if lambda(k) < 0
        mua(k) = -mua(k);
    end
end

% Determine radius from the center of the planet to the surface of the
% planet by using trigonometric relationship
ra = xa./cos(lambda);

% Determine length from the surface of the planet to the body 
l = r - ra;

% The angular difference between geocentric latitude and geodetic latitude at
% surface of the planet
dlambda = mua - lambda;

% Estimate mean sea-level altitude
h = l.*cos(dlambda);

% Find the radius of curvature in the Meridian
den = 1 - (2.*f - f.*f).*sin(mua).*sin(mua);
rhoa=(R.*( 1 - f ).*( 1 - f ))./sqrt(den.*den.*den);

% the angular difference between geodetic latitude at the surface of the
% planet and geodetic latitude at the body 
dmu = atan2(l.*sin(dlambda),rhoa + h);

gd = convang(mua - dmu,'rad','deg');
