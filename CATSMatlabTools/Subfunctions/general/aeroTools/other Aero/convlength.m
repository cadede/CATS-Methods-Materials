function out = convlength( in, uin, uout )
%  CONVLENGTH Convert from length units to desired length units.
%   A = CONVLENGTH( V, UI, UO ) computes the conversion factor from
%   specified input length units, UI, to specified output length units, UO, and
%   applies the conversion factor to the input, V, to produce the output,
%   A, in the desired units.  V and A are floating point arrays of size
%   M-by-N.	All of the values in V must have the same unit conversions from
%   UI to UO.  UI and UO are strings.
%
%   Supported unit strings are:
%      'ft'        :feet       
%      'm'         :meters
%      'km'        :kilometers
%      'in'        :inches
%      'mi'        :miles  
%      'naut mi'   :nautical miles
%
%   Example:
%
%   Convert three lengths from feet to meters:
%      a = convlength([3 10 20], 'ft','m')
%
%   See also CONVACC, CONVANG, CONVANGACC, CONVANGVEL, CONVDENSITY,
%   CONVFORCE, CONVMASS, CONVPRES, CONVTEMP, CONVVEL.

%   Copyright 2000-2010 The MathWorks, Inc.

if ~isfloat( in )
    error(message('aero:convlength:notFloat'));
end

slope = unitconversion('length conversion',uin,uout,0);

out = in.*slope;

