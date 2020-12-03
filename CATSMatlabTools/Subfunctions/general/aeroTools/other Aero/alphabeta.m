function [alpha, beta] = alphabeta( vel )
%  ALPHABETA Compute incidence and sideslip angles.
%   [A, B] = ALPHABETA( V ) computes M incidence and sideslip angles, A
%   and B, between the velocity vector and the body from an M-by-3 or array
%   of velocities, V in body-axes.  A and B are in radians. 
%
%   Examples:
%
%   Determine incidence and sideslip angles for velocity in feet per
%   second: 
%      [alpha, beta] = alphabeta( [84.3905  33.7562  10.1269] )
%
%   Determine incidence and sideslip angles for velocity in knots: 
%      [alpha, beta] = alphabeta( [ 50 20 6; 5 0.5 2] )
%
%   See also AIRSPEED, MACHNUMBER.

%   Copyright 2000-2018 The MathWorks, Inc.

if ~isnumeric(vel)
    error(message('aero:alphabeta:notNumeric'));
end

if (size(vel,2)==3)
    alpha  = atan2(vel(:,3), vel(:,1));
    beta = asin(vel(:,2)./airspeed(vel));
else
    error(message('aero:alphabeta:wrongDimension'));
end

