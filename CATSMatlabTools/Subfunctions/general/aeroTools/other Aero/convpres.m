function out = convpres( in, uin, uout )
%  CONVPRES Convert from pressure units to desired pressure units.
%   A = CONVPRES( V, UI, UO ) computes the conversion factor from
%   specified input pressure units, UI, to specified output pressure units,
%   UO, and applies the conversion factor to the input, V, to produce the
%   output, A, in the desired units.  V and A are floating point arrays of
%   size M-by-N.	All of the values in V must have the same unit
%   conversions from UI to UO.  UI and UO are strings.
%
%   Supported unit strings are:
%      'psi'         :pound force per square inch      
%      'Pa'          :pascal
%      'psf'         :pound force per square foot
%      'atm'         :atmosphere
%
%   Example:
%
%   Convert two pressures from pound force per square inch to atmospheres:  
%      a = convpres([14.696  35], 'psi','atm')
%
%   See also CONVACC, CONVANG, CONVANGACC, CONVANGVEL, CONVDENSITY,
%   CONVFORCE, CONVLENGTH, CONVMASS, CONVTEMP, CONVVEL.

%   Copyright 2000-2010 The MathWorks, Inc.

if ~isfloat( in )
    error(message('aero:convpres:notFloat'));
end

slope = unitconversion('pressure conversion',uin,uout,0);

out = in.*slope;
