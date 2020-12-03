function [T, a, P, rho] = atmosisa( h )
%  ATMOSISA Use International Standard Atmosphere Model.
%   [T, A, P, RHO] = ATMOSISA(H) implements the mathematical representation
%   of the International Standard Atmosphere values for ambient temperature, 
%   pressure, density, and speed of sound for the input geopotential altitude. 
%
%   Input required by ATMOSISA is:
%   H      :a numeric array of M geopotential height in meters. 
%
%   Output calculated for the International Standard Atmosphere are: 
%   T      :a numeric array of M temperature in kelvin.
%   a      :a numeric array of M speed of sound in meters per second.
%   P      :a numeric array of M air pressure in pascal.
%   rho    :a numeric array of M air density in kilograms per meter cubed.
%
%   Limitation: 
%
%   Below the geopotential altitude of 0 km and above the geopotential
%   altitude of the tropopause, temperature and pressure values are held.
%   Density and speed of sound are calculated using a perfect gas
%   relationship.
%
%   Examples:
%
%   Calculate the International Standard Atmosphere at 1000 meters:
%      [T, a, P, rho] = atmosisa(1000)
%
%   Calculate the International Standard Atmosphere at 1000, 11000 and 20000
%   meters:
%      [T, a, P, rho] = atmosisa([1000 11000 20000])
%
%   See also ATMOSCIRA, ATMOSCOESA, ATMOSLAPSE, ATMOSNONSTD.

%   Copyright 2000-2010 The MathWorks, Inc.

%   Reference:  U.S. Standard Atmosphere, 1976, U.S. Government Printing 
%   Office, Washington, D.C.

[T, a, P, rho] = atmoslapse(h, 9.80665, 1.4, 287.0531, 0.0065, 11000, 20000, ...
    1.225, 101325, 288.15 );
