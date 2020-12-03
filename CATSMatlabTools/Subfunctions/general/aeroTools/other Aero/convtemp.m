function out = convtemp( in, uin, uout )
%  CONVTEMP Convert from temperature units to desired temperature units.
%   A = CONVTEMP( V, UI, UO ) computes the conversion factor from
%   specified input temperature units, UI, to specified output temperature
%   units, UO, and applies the conversion factor to the input, V, to
%   produce the output, A, in the desired units.  V and A are floating
%   point arrays of size M-by-N.	All of the values in V must have the
%   same unit conversions from UI to UO.  UI and UO are strings.
%
%   Supported unit strings are:
%      'K'           :kelvin       
%      'F'           :degrees Fahrenheit
%      'C'           :degrees Celsius
%      'R'           :degrees Rankine
%
%   Example:
%
%   Convert three temperatures from degrees Celsius to degrees Fahrenheit: 
%      a = convtemp([0 100 15], 'C','F')
%
%   See also CONVACC, CONVANG, CONVANGACC, CONVANGVEL, CONVDENSITY,
%   CONVFORCE, CONVLENGTH, CONVMASS, CONVPRES, CONVVEL.

%   Copyright 2000-2010 The MathWorks, Inc.

if ~isfloat( in )
    error(message('aero:convtemp:notFloat'));
end

[slope bias] = unitconversion('temperature conversion',uin,uout,0);

out = in.*slope + bias;
