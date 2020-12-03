function [ XYZ, H, DEC, DIP, F, DXDYDZ, DH, DDEC, DDIP, DF] = igrf11magm(height,lat,lon,dyear)
% IGRF11MAGM Use 11th generation of International Geomagnetic Reference Field
%  [XYZ, H, DEC, DIP, F, DXDYDZ, DH, DDEC, DDIP, DF] = 
%                                     IGRF11MAGM( HEIGHT, LAT, LON, DYEAR )
%  calculates the Earth's magnetic field and the secular variation at a specific 
%  location and time using the 11th generation of International Geomagnetic 
%  Reference Field (IGRF-11). 
%
%  Inputs required by IGRF-11 are:
%   HEIGHT :a scalar value in meters. 
%   LAT    :a scalar geodetic latitude in degrees where north latitude is 
%          positive, and south latitude is negative.
%   LON    :a scalar geodetic longitude in degrees where east longitude 
%          is positive, west is negative.
%   DYEAR  :a scalar decimal year.  Decimal year is the desired year in 
%          a decimal format to include any fraction of the year that has 
%          already passed.
%
%  Output calculated for the Earth's magnetic field and secular variation include:
%   XYZ    :magnetic field vector in nanotesla (nT). Z is vertical component (+ve down)
%   H      :horizontal intensity in nanotesla (nT).
%   DEC    :declination in degrees. (+ve east)
%   DIP    :inclination in degrees. (+ve down)
%   F      :total intensity in nanotesla (nT).
%   DXDYDZ :secular variation in magnetic field vector in nT/year. Z is
%          vertical component (+ve down) 
%   DH     :secular variation in horizontal intensity in nT/year.
%   DDEC   :secular variation in declination in minutes/year. (+ve east)
%   DDIP   :secular variation in inclination in minutes/year. (+ve down)
%   DF     :secular variation in total intensity in nT/year.
%
%   Limitations:
%
%   The igrf11magm function is no longer supported. A version of this 
%   function is supplied for backward compatiblity. This function will be 
%   removed in a future release. Use IGRFMAGM instead.
%
%   This function is valid between the heights of -1000 meters to 600000
%   meters. 
%
%   This function is valid between the years of 1900 and 2015.
%
%   This function has the limitations of the International Geomagnetic
%   Reference Field (IGRF). For more information see the IGRF web site,
%   http://www.ngdc.noaa.gov/IAGA/vmod/igrfhw.html.   
%
%   Example:
%
%   Calculate the magnetic model 1000 meters over Natick, Massachusetts on 
%   July 4, 2005 using IGRF-11:
%      [XYZ, H, DEC, DIP, F] = igrf11magm(1000, 42.283, -71.35, decyear(2005,7,4))
%
%   See also DECYEAR, WRLDMAGM, IGRFMAGM.

%   Copyright 2010-2015 The MathWorks, Inc.

%   Reference:
%   [1] The IGRF-11 can be found on the web at 
%       http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
%   [2] Blakely, R. J., "Potential Theory in Gravity & Magnetic Applications", 
%       Cambridge University Press, Cambridge UK, 1996. 
narginchk(4, 4);
[ XYZ, H, DEC, DIP, F, DXDYDZ, DH, DDEC, DDIP, DF] = igrfmagm(height,lat,lon,dyear,11);