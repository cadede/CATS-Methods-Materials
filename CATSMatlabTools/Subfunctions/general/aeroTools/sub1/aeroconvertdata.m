 function ConvertToMKS = aeroconvertdata( mtype )
% AEROCONVERTDATA internal function containing data for unit conversion.

%  All conversion factors convert specified units to metric  
  
%   Copyright 2000-2010 The MathWorks, Inc.

%
% Create structure for conversions
%

switch lower(mtype)
case 'temperature conversion'
    
    ConvertToMKS.type = mtype; % to deg C
    ConvertToMKS.tdata(1).unit  = 'K';
    ConvertToMKS.tdata(1).slope = 1.0;
    ConvertToMKS.tdata(1).bias  = -273.15;
    ConvertToMKS.tdata(2).unit  = 'F';
    ConvertToMKS.tdata(2).slope = 5.0/9.0;
    ConvertToMKS.tdata(2).bias  = -32.0*ConvertToMKS.tdata(2).slope;
    ConvertToMKS.tdata(3).unit  = 'C';
    ConvertToMKS.tdata(3).slope = 1.0;
    ConvertToMKS.tdata(3).bias  = 0.0;
    ConvertToMKS.tdata(4).unit  = 'R';
    ConvertToMKS.tdata(4).slope = 5.0/9.0;
    ConvertToMKS.tdata(4).bias  = -273.15;
    
    %
    % Create structure for MKS conversion
    % -- accurate to 9 decimal places
    %
case 'angle conversion'
    
    ConvertToMKS.type = mtype; % to rad
    ConvertToMKS.tdata(1).unit  = 'deg';
    ConvertToMKS.tdata(1).slope = pi/180.0;
    ConvertToMKS.tdata(2).unit  = 'rad';
    ConvertToMKS.tdata(2).slope = 1.0;
    ConvertToMKS.tdata(3).unit  = 'rev';
    ConvertToMKS.tdata(3).slope = 2*pi;
    
case 'length conversion'
    
    ConvertToMKS.type = mtype; % to m
    ConvertToMKS.tdata(1).unit = 'ft';
    ConvertToMKS.tdata(1).slope = 0.3048;
    ConvertToMKS.tdata(2).unit = 'm';
    ConvertToMKS.tdata(2).slope = 1.0;
    ConvertToMKS.tdata(3).unit = 'km';
    ConvertToMKS.tdata(3).slope = 1000.0;
    ConvertToMKS.tdata(4).unit = 'in';
    ConvertToMKS.tdata(4).slope = 0.0254;
    ConvertToMKS.tdata(5).unit = 'mi';
    ConvertToMKS.tdata(5).slope = 1609.344;
    ConvertToMKS.tdata(6).unit = 'naut mi';
    ConvertToMKS.tdata(6).slope = 1852.0;
    
case 'angular velocity conversion'
    
    ConvertToMKS.type = mtype; % to rad/s
    ConvertToMKS.tdata(1).unit = 'deg/s';
    ConvertToMKS.tdata(1).slope = pi/180.0;
    ConvertToMKS.tdata(2).unit = 'rad/s';
    ConvertToMKS.tdata(2).slope = 1.0;
    ConvertToMKS.tdata(3).unit = 'rpm';
    ConvertToMKS.tdata(3).slope = 0.10471975512;
    
case 'velocity conversion'
    
    ConvertToMKS.type = mtype; % to m/s
    ConvertToMKS.tdata(1).unit = 'ft/s';
    ConvertToMKS.tdata(1).slope = 0.3048;
    ConvertToMKS.tdata(2).unit = 'm/s';
    ConvertToMKS.tdata(2).slope = 1.0;
    ConvertToMKS.tdata(3).unit = 'km/s';
    ConvertToMKS.tdata(3).slope = 1000.0;
    ConvertToMKS.tdata(4).unit = 'in/s';
    ConvertToMKS.tdata(4).slope = 0.0254;
    ConvertToMKS.tdata(5).unit = 'km/h';
    ConvertToMKS.tdata(5).slope = 0.2777777777778;
    ConvertToMKS.tdata(6).unit = 'mph';
    ConvertToMKS.tdata(6).slope = 0.44704;
    ConvertToMKS.tdata(7).unit = 'kts';
    ConvertToMKS.tdata(7).slope =  0.514444444444;
    ConvertToMKS.tdata(8).unit = 'ft/min';
    ConvertToMKS.tdata(8).slope = 0.3048/60;

case 'angular acceleration conversion'
    
    ConvertToMKS.type = mtype; % to rad/s^2
    ConvertToMKS.tdata(1).unit = 'deg/s^2';
    ConvertToMKS.tdata(1).slope = pi/180.0;
    ConvertToMKS.tdata(2).unit = 'rad/s^2';
    ConvertToMKS.tdata(2).slope = 1.0;
    ConvertToMKS.tdata(3).unit = 'rpm/s';
    ConvertToMKS.tdata(3).slope = 0.10471975512;
    
case 'acceleration conversion'
    
    ConvertToMKS.type = mtype; % to m/s^2
    ConvertToMKS.tdata(1).unit = 'ft/s^2';
    ConvertToMKS.tdata(1).slope = 0.3048;
    ConvertToMKS.tdata(2).unit = 'm/s^2';
    ConvertToMKS.tdata(2).slope = 1.0;
    ConvertToMKS.tdata(3).unit = 'km/s^2';
    ConvertToMKS.tdata(3).slope = 1000.0;
    ConvertToMKS.tdata(4).unit = 'in/s^2';
    ConvertToMKS.tdata(4).slope = 0.0254;
    ConvertToMKS.tdata(5).unit = 'km/h-s';
    ConvertToMKS.tdata(5).slope = 0.2777777777778;
    ConvertToMKS.tdata(6).unit = 'mph/s';
    ConvertToMKS.tdata(6).slope = 0.44704;
    ConvertToMKS.tdata(7).unit = 'G''s';
    ConvertToMKS.tdata(7).slope = 9.80665;
    
case 'mass conversion'
    
    ConvertToMKS.type = mtype; % to kg
    ConvertToMKS.tdata(1).unit = 'lbm';
    ConvertToMKS.tdata(1).slope = 0.45359237;
    ConvertToMKS.tdata(2).unit = 'kg';
    ConvertToMKS.tdata(2).slope = 1.0;
    ConvertToMKS.tdata(3).unit = 'slug';
    ConvertToMKS.tdata(3).slope = 14.5939029;
    
case 'force conversion'
    
    ConvertToMKS.type = mtype; % to N
    ConvertToMKS.tdata(1).unit = 'lbf';
    ConvertToMKS.tdata(1).slope = 4.448222;
    ConvertToMKS.tdata(2).unit = 'N';
    ConvertToMKS.tdata(2).slope = 1.0;
    
case 'pressure conversion'
    
    ConvertToMKS.type = mtype; % to Pa
    ConvertToMKS.tdata(1).unit = 'psi';
    ConvertToMKS.tdata(1).slope = 6894.757889516;
    ConvertToMKS.tdata(2).unit = 'Pa';
    ConvertToMKS.tdata(2).slope = 1.0;
    ConvertToMKS.tdata(3).unit = 'psf';
    ConvertToMKS.tdata(3).slope = 47.8802631216;
    ConvertToMKS.tdata(4).unit = 'atm';
    ConvertToMKS.tdata(4).slope = 101325.0;
    
case 'density conversion'
    
    ConvertToMKS.type = mtype; % to kg/m^3
    ConvertToMKS.tdata(1).unit = 'lbm/ft^3';
    ConvertToMKS.tdata(1).slope = 16.018463374;
    ConvertToMKS.tdata(2).unit = 'kg/m^3';
    ConvertToMKS.tdata(2).slope = 1.0;
    ConvertToMKS.tdata(3).unit = 'slug/ft^3';
    ConvertToMKS.tdata(3).slope = 515.378817079;
    ConvertToMKS.tdata(4).unit = 'lbm/in^3';
    ConvertToMKS.tdata(4).slope = 27679.904710203;
    
    otherwise
        
    error(message('aero:aeroconvertdata:unknownUnit', mtype));
end
