function out  = convacc( in, uin, uout )
%  CONVACC Convert from acceleration units to desired acceleration units.
%   A = CONVACC( V, UI, UO ) computes the conversion factor from
%   specified input acceleration units, UI, to specified output
%   acceleration units, UO, and applies the conversion factor to the input,
%   V, to produce the output,  A, in the desired units.  V and A are
%   floating point arrays of size M-by-N.	All of the values in V must
%   have the same unit conversions from UI to UO.  UI and UO are strings. 
%
%   Supported unit strings are:
%      'ft/s^2'      :feet per second squared
%      'm/s^2'       :meters per second squared
%      'km/s^2'      :kilometers per second squared
%      'in/s^2';     :inches per second squared
%      'km/h-s';     :kilometers per hour per second
%      'mph/s'       :miles per hour per second
%      'G''s'        :gees
%
%   Example:
%
%   Convert three accelerations from feet per second squared to meters per
%   second squared:
%      a = convacc([3 10 20], 'ft/s^2','m/s^2')
%
%   See also CONVANG, CONVANGACC, CONVANGVEL, CONVDENSITY, CONVFORCE,
%   CONVLENGTH, CONVMASS, CONVPRES, CONVTEMP, CONVVEL. 

%   Copyright 2000-2010 The MathWorks, Inc.

if ~isfloat( in )
    error(message('aero:convacc:notFloat'));
end

slope = unitconversion('acceleration conversion',uin,uout,0);

out = in.*slope;
