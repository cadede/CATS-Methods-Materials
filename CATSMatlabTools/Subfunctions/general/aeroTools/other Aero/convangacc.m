function out = convangacc( in, uin, uout )
%  CONVANGACC Convert from angular acceleration units to desired angular
%  acceleration units.
%   A = CONVANGACC( V, UI, UO ) computes the conversion factor from
%   specified input angular acceleration units, UI, to specified output
%   angular acceleration units, UO, and applies the conversion factor to
%   the input, V, to produce the output, A, in the desired units.  V and A
%   are floating point arrays of size M-by-N.	All of the values in V must
%   have the same unit conversions from UI to UO.  UI and UO are strings. 
%
%   Supported unit strings are:
%      'deg/s^2'     :degrees per second squared       
%      'rad/s^2'     :radians per second squared
%      'rpm/s'       :revolutions per minute per second
%
%   Example:
%
%   Convert three angular accelerations from degrees per second squared to
%   radians per second squared:
%      a = convangacc([0.3 0.1 0.5], 'deg/s^2','rad/s^2')
%
%   See also CONVACC, CONVANG, CONVANGVEL, CONVDENSITY, CONVFORCE, 
%   CONVLENGTH, CONVMASS, CONVPRES, CONVTEMP, CONVVEL.

%   Copyright 2000-2010 The MathWorks, Inc.

if ~isfloat( in )
    error(message('aero:convangacc:notFloat'));
end

slope = unitconversion('angular acceleration conversion',uin,uout,0);

out = in.*slope;
