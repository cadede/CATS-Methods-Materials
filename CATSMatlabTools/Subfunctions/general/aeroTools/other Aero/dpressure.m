function q = dpressure( vel, rho )
%  DPRESSURE Compute dynamic pressure using velocity and density.
%   Q = DPRESSURE( V, R ) computes M dynamic pressures, Q, from an M-by-3
%   array of velocities, V, and an array of M densities, R.  V
%   and R must have the same length units.
%
%   Examples:
%
%   Determine dynamic pressure for velocity in feet per
%   second and density in slugs per feet cubed:
%      q = dpressure( [84.3905  33.7562  10.1269], 0.0024 )
%
%   Determine dynamic pressure for velocity  in meters per
%   second and density in kilograms per meters cubed:
%      q = dpressure( [25.7222 10.2889 3.0867], [ 1.225  0.3639 ])
%
%   Determine dynamic pressure for velocity in meters per second and
%   density in kilograms per meters cubed: 
%      q = dpressure( [ 50 20 6; 5 0.5 2], [1.225  0.3639])
%
%   See also AIRSPEED, MACHNUMBER.

%   Copyright 2000-2011 The MathWorks, Inc.

if ~isnumeric( vel )
    error(message('aero:dpressure:notNumericVelocity'));
end

if ~isnumeric( rho )
    error(message('aero:dpressure:notNumericDensity'));
end

% airspeed checks velocity dimension
v = airspeed(vel);

if ~( isscalar(v) || isscalar(rho) || (length(v)== length(rho)) )
    error(message('aero:dpressure:wrongDimension'));
end

q = 0.5*rho(:).*v.^2;

