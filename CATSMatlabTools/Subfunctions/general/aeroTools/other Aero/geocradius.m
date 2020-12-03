function radius = geocradius( lambda, varargin )
%  GEOCRADIUS Estimate radius of ellipsoid planet at geocentric latitude.
%   R = GEOCRADIUS( LAMBDA ) estimates the radius, R, of an ellipsoid
%   planet at a particular geocentric latitude, LAMBDA.  LAMBDA is in
%   degrees.  R is in meters.  The default ellipsoid planet is WGS84.
%
%   R = GEOCRADIUS( LAMBDA, MODEL ) is an alternate method for estimating
%   the radius, for a specific ellipsoid planet.  Currently only 'WGS84' is
%   supported for MODEL.
%
%   R = GEOCRADIUS( LAMBDA, F, RE ) is another alternate method for
%   estimating the radius for a custom ellipsoid planet defined by
%   flattening, F, and the equatorial radius, RE in meters. 
%
%   Examples:
%
%   Determine radius at 45 degrees latitude:
%      r = geocradius( 45 )
%
%   Determine radius at multiple latitudes:
%      r = geocradius( [0 45 90] )
%
%   Determine radius at multiple latitudes specifying WGS84 ellipsoid
%   model:
%      r = geocradius( [0 45 90], 'WGS84' )
%
%   Determine radius at multiple latitudes specifying custom ellipsoid
%   model:
%      f = 1/196.877360;
%      Re = 3397000;
%      r = geocradius( [0 45 90],  f, Re )
%
%   See also GEOC2GEOD, GEOD2GEOC.

%   Copyright 2000-2011 The MathWorks, Inc.

%   References: 
%   [1] Stevens, B. L., and F. L. Lewis, "Aircraft Control and Simulation,"
%   John Wiley & Sons, New York, NY, 1992.
%   [2] Zipfel, P. H., "Modeling and Simulation of Aerospace Vehicle
%   Dynamics," AIAA Education Series, Reston, VA, 2000.

narginchk(1, 3);

if ~isnumeric( lambda )
    error(message('aero:geocradius:notNumeric'));
end

[f R] = worldparams(nargin, varargin);

radlam = convang( lambda , 'deg','rad');
sinlam = sin( radlam );
radius  = sqrt(( R.^2 )./( 1 + (1/(( 1 - f ).^2) - 1).*sinlam.^2 ));
