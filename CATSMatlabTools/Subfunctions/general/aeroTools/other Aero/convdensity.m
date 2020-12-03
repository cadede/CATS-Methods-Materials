function out = convdensity( in, uin, uout )
%  CONVDENSITY Convert from density units to desired density units.
%   A = CONVDENSITY( V, UI, UO ) computes the conversion factor from
%   specified input density units, UI, to specified output density units,
%   UO, and applies the conversion factor to the input, V, to produce the
%   output, A, in the desired units.  V and A are floating point arrays of
%   size M-by-N.	All of the values in V must have the same unit
%   conversions from UI to UO.  UI and UO are strings.
%
%   Supported unit strings are:
%      'lbm/ft^3'    :pound mass per feet cubed      
%      'kg/m^3'      :kilograms per meters cubed
%      'slug/ft^3'   :slugs per feet cubed
%      'lbm/in^3'    :pound mass per inches cubed
%
%   Example:
%
%   Convert three densities from pound mass per feet cubed to kilograms per
%   meters cubed: 
%      a = convdensity([0.3 0.1 0.5], 'lbm/ft^3','kg/m^3')
%
%   See also CONVACC, CONVANG, CONVANGACC, CONVANGVEL, CONVFORCE,
%   CONVLENGTH, CONVMASS, CONVPRES, CONVTEMP, CONVVEL. 

%   Copyright 2000-2010 The MathWorks, Inc.

if ~isfloat( in )
    error(message('aero:convdensity:notFloat'));
end

slope = unitconversion('density conversion',uin,uout,0);

out = in.*slope;

