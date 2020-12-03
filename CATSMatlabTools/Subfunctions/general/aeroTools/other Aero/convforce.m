function out = convforce( in, uin, uout )
%  CONVFORCE Convert from force units to desired force units.
%   A = CONVFORCE( V, UI, UO ) computes the conversion factor from
%   specified input force units, UI, to specified output force units, UO, and
%   applies the conversion factor to the input, V, to produce the output,
%   A, in the desired units.  V and A are floating point arrays of size
%   M-by-N.	All of the values in V must have the same unit conversions from
%   UI to UO.  UI and UO are strings.
%
%   Supported unit strings are:
%      'lbf'         :pound force     
%      'N'           :newton
%
%   Example:
%
%   Convert three forces from pound force to newtons: 
%      a = convforce([120 1 5], 'lbf','N')
%
%   See also CONVACC, CONVANG, CONVANGACC, CONVANGVEL, CONVDENSITY,
%   CONVLENGTH, CONVMASS, CONVPRES, CONVTEMP, CONVVEL.

%   Copyright 2000-2010 The MathWorks, Inc.

if ~isfloat( in )
    error(message('aero:convforce:notFloat'));
end

slope = unitconversion('force conversion',uin,uout,0);

out = in.*slope;
