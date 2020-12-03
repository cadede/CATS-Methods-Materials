function out = convangvel( in, uin, uout )
%  CONVANGVEL Convert from angular velocity units to desired angular
%  velocity units. 
%   N = CONVANGVEL( V, UI, UO ) computes the conversion factor from
%   specified input angular velocity units, UI, to specified output angular
%   velocity units, UO, and applies the conversion factor to the input, V,
%   to produce the output, N, in the desired units.  V and N are floating
%   point vectors of length M.	All of the values in V must have the same
%   unit conversions from UI to UO.  UI and UO are strings.
%
%   Supported unit strings are:
%      'deg/s'     :degrees per second        
%      'rad/s'     :radians per second 
%      'rpm'       :revolutions per minute
%
%   Example:
%
%   Convert three angular velocities from degrees per second to radians per
%   second: 
%      a = convangvel([0.3 0.1 0.5], 'deg/s','rad/s')
%
%   See also CONVACC, CONVANG, CONVANGACC, CONVDENSITY, CONVFORCE,
%   CONVLENGTH, CONVMASS, CONVPRES, CONVTEMP, CONVVEL. 

%   Copyright 2000-2010 The MathWorks, Inc.

if ~isfloat( in )
    error(message('aero:convangvel:notFloat'));
end

slope = unitconversion('angular velocity conversion',uin,uout,0);

out = in.*slope;
