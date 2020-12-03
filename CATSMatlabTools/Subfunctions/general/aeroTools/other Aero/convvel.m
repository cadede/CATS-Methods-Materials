function out = convvel( in, uin, uout )
%  CONVVEL Convert from velocity units to desired velocity units.
%   A = CONVVEL( V, UI, UO ) computes the conversion factor from
%   specified input velocity units, UI, to specified output velocity units,
%   UO, and applies the conversion factor to the input, V, to produce the
%   output, A, in the desired units.  V and A are floating point arrays of
%   size M-by-N.	All of the values in V must have the same unit
%   conversions from UI to UO.  UI and UO are strings.
%
%   Supported unit strings are:
%      'ft/s'        :feet per second      
%      'm/s'         :meters per second
%      'km/s'        :kilometers per second
%      'in/s'        :inches per second
%      'km/h'        :kilometers per hour
%      'mph'         :miles per hour
%      'kts'         :knots
%      'ft/min'      :feet per minute      
%
%   Example:
%
%   Convert three velocities from feet per minute to meters per second:
%      a = convvel([30 100 250], 'ft/min','m/s')
%
%   See also CONVACC, CONVANG, CONVANGACC, CONVANGVEL, CONVDENSITY,
%   CONVFORCE, CONVLENGTH, CONVMASS, CONVPRES, CONVTEMP.

%   Copyright 2000-2010 The MathWorks, Inc.

if ~isfloat( in )
    error(message('aero:convvel:notFloat'));
end

slope = unitconversion('velocity conversion',uin,uout,0);

out = in.*slope;
