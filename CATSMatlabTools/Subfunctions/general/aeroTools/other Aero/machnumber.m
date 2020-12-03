function mach = machnumber( vel, a )
%  MACHNUMBER Compute Mach number using velocity and speed of sound.
%   MACH = MACHNUMBER( V, A ) computes M Mach numbers, MACH, from an M-by-3
%   array of velocities, V and a array of M speeds of sound, A.  V and A
%   must have the same velocity units.  
%
%   Examples:
%
%   Determine mach number for velocity and speed of sound in feet per
%   second:
%      mach = machnumber( [84.3905  33.7562  10.1269], 1116.4505 )
%
%   Determine mach number for velocity and speed of sound in meters per
%   second:
%      mach = machnumber( [25.7222 10.2889 3.0867], [ 340.2941  295.0696 ])
%
%   Determine mach number for velocity and speed of sound in knots:
%      mach = machnumber( [ 50 20 6; 5 0.5 2], [ 661.4789  573.5694])
%
%   See also AIRSPEED, ALPHABETA, DPRESSURE.

%   Copyright 2000-2011 The MathWorks, Inc.

if ~isnumeric( vel )
    error(message('aero:machnumber:notNumeric1'));
end

if ~isnumeric( a )
    error(message('aero:machnumber:notNumeric2'));
end

if (size(vel,2)~=3)
    error(message('aero:machnumber:wrongDimensionVelocity'));
end

if ( size(vel,1)==size(a,2) )
    a = a';
elseif  ~( isscalar(a) || (size(vel,1)==size(a,1)) )
    error(message('aero:machnumber:wrongDimensionSoS'));
end

mach = airspeed(vel)./a;
