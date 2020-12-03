function [T, a, P, rho] = atmoslapse(h, g, gamma, R, L, hts, htp, rho0, P0, T0, varargin)
%   ATMOSLAPSE Use Lapse Rate Atmosphere Model. 
%   [T, A, P, RHO] = ATMOSLAPSE(H, G, GAMMA, R, L, HTS, HTP, RHO0, P0, T0) 
%   implements the mathematical representation of the lapse rate atmospheric 
%   equations for ambient temperature, pressure, density, and speed of sound 
%   for the input geopotential altitude. To customize this atmospheric model, 
%   specify the atmospheric properties in the function input.
%
%   Below the geopotential altitude of 0 m and above the geopotential 
%   altitude of the tropopause, temperature and pressure values are held. 
%   Density and speed of sound are calculated using a perfect gas 
%   relationship.
%
%   [T, A, P, RHO] = ATMOSLAPSE(H, G, GAMMA, R, L, HTS, HTP, RHO0, P0, T0, H0)  
%   indicates that the values for ambient temperature, pressure, density, and 
%   speed of sound are required for below mean sea level geopotential 
%   altitudes. 
%
%   Below the geopotential altitude of H0 m and above the geopotential 
%   altitude of the tropopause, temperature and pressure values are held. 
%   Density and speed of sound are calculated using a perfect gas
%   relationship.
% 
%   Inputs required by ATMOSLAPSE are: 
%   H       :Numeric array of M geopotential height in meters. 
%   G       :Scalar of acceleration due to gravity in meters per second squared. 
%   GAMMA   :Scalar of specific heat ratio. 
%   R       :Scalar of characteristic gas constant joules per kilogram-kelvin. 
%   L       :Scalar of lapse rate in kelvin per meter. 
%   HTS     :Scalar of height of troposphere in meters. 
%   HTP     :Scalar of height of tropopause in meters. 
%   RHO0    :Scalar of air density at mean sea level in kilograms per meter cubed. 
%   P0      :Scalar of static pressure at mean sea level in pascal. 
%   T0      :Scalar of absolute temperature at mean sea level in kelvin. 
%   H0      :Scalar of minimum sea level altitude in m. 
% 
%   Outputs calculated for the lapse rate atmosphere are: 
%   T       :Numeric array of M temperature in kelvin. 
%   A       :Numeric array of M speed of sound in meters per second. 
%   P       :Numeric array of M air pressure in pascal. 
%   RHO     :Numeric array of M air density in kilograms per meter cubed. 
%
%   Example: 
% 
%   Calculate the atmosphere at 1000 meters with the International Standard 
%   Atmosphere input values: 
%       [T, A, P, RHO] = atmoslapse(1000, 9.80665, 1.4, 287.0531, 0.0065, ... 
%           11000, 20000, 1.225, 101325, 288.15 ) 
% 
%   See also ATMOSCIRA, ATMOSCOESA, ATMOSISA, ATMOSNONSTD. 

%   Copyright 2000-2012 The MathWorks, Inc. 

%   Reference: U.S. Standard Atmosphere, 1976, U.S. Government Printing 
%   Office, Washington, D.C. 


narginchk(10,11);

if nargin==11
    H0 = varargin{1};
else
    H0 = 0;
end

if ~(isscalar(g)&&isscalar(gamma)&&isscalar(R)&&isscalar(L)&&isscalar(hts) ...
     &&isscalar(htp)&&isscalar(rho0)&&isscalar(P0)&&isscalar(T0)&&isscalar(H0))
    error(message('aero:atmoslapse:nonScalar'));
end

if ~isnumeric( h )
    error(message('aero:atmoslapse:notNumeric'));
end

for i = length(h): -1 :1
    if ( h(i) > htp )
        h(i) = htp;
    end

    if ( h(i) <  H0 )
        h(i) = H0;
    end

    if ( h(i) > hts )
        T(i) = T0 - L*hts; 
        expon(i) = exp(g/(R*T(i))*(hts - h(i))); 
    else
        T(i) = T0 - L*h(i); 
        expon(i) = 1.0;  
    end
end

a = sqrt(T*gamma*R);

theta = T/T0;

P = P0*theta.^(g/(L*R)).*expon;
rho = rho0*theta.^((g/(L*R))-1.0).*expon;
