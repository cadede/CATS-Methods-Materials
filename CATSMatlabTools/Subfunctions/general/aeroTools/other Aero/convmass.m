function out = convmass( in, uin, uout )
%  CONVMASS Convert from mass units to desired mass units.
%   A = CONVMASS( V, UI, UO ) computes the conversion factor from
%   specified input mass units, UI, to specified output mass units, UO, and
%   applies the conversion factor to the input, V, to produce the output,
%   A, in the desired units.  V and A are floating point arrays of size
%   M-by-N.	All of the values in V must have the same unit conversions from
%   UI to UO.  UI and UO are strings.
%
%   Supported unit strings are:
%      'lbm'    :pound mass     
%      'kg'     :kilograms
%      'slug'   :slugs
%
%   Example:
%
%   Convert three masses from pound mass to kilograms: 
%      a = convmass([3 1 5], 'lbm','kg')
%
%   See also CONVACC, CONVANG, CONVANGACC, CONVANGVEL, CONVDENSITY,
%   CONVFORCE, CONVLENGTH, CONVPRES, CONVTEMP, CONVVEL.

%   Copyright 2000-2010 The MathWorks, Inc.

if ~isfloat( in )
    error(message('aero:convmass:notFloat'));
end

slope = unitconversion('mass conversion',uin,uout,0);

out = in.*slope;
