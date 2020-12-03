function as = airspeed( vel )
%  AIRSPEED Compute airspeed from velocity.
%   AS = AIRSPEED( V ) computes M airspeeds, AS, from an M-by-3 array of
%   velocities, V.   
%
%   Examples:
%
%   Determine airspeed for velocity in feet per second:
%      as = airspeed( [84.3905  33.7562  10.1269] )
%
%   Determine airspeed for velocity in knots:
%      as = airspeed( [ 50 20 6; 5 0.5 2])
%
%   See also ALPHABETA, CORRECTAIRSPEED, DPRESSURE, MACHNUMBER. 

%   Copyright 2000-2010 The MathWorks, Inc.

if ~isnumeric( vel )
    error(message('aero:airspeed:notNumeric'));
end

if (size(vel,2)==3)
    as = sqrt(vel(:,1).^2 + vel(:,2).^2 + vel(:,3).^2);
else
    error(message('aero:airspeed:wrongDimension'));
end
