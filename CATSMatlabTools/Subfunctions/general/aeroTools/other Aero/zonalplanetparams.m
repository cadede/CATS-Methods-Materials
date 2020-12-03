function [ Re, GM, Jvalues, omega ] = zonalplanetparams( model )
% ZONALPLANETPARAMS internal function containing planetary parameters for
% zonal harmonic and centrifugal gravity.
  
%   Copyright 2009 The MathWorks, Inc.

%   References:  
%   [1] Vallado, D. A., "Fundamentals of Astrodynamics and Applications",
%       McGraw-Hill, New York, 1997.  

switch lower( model )
    case 'moon'
        Re      = 1738.0e3;   % m
        GM      = 4902.799e9; % m^3/s^2
        Jvalues = 0.0002027;
        siderealRotRate = 27.32166; % days 
    case 'mercury'
        Re      = 2439.0e3;   % m
        GM      = 2.2032e13;  % m^3/s^2
        Jvalues = 0.00006;
        siderealRotRate = 58.6462; % days 
    case 'venus'
        Re      = 6052.0e3;   % m
        GM      = 3.257e14;   % m^3/s^2
        Jvalues = 0.000027;
        siderealRotRate = -243.01; % days 
    case 'mars'
        Re      = 3397.2e3;   % m
        GM      = 4.305e13;   % m^3/s^2
        Jvalues = [ 0.001964 0.000036 ];
        siderealRotRate = 1.02595675; % days 
    case 'jupiter'
        Re      = 71492.e3;   % m
        GM      = 1.268e17;   % m^3/s^2
        Jvalues = [0.01475 0 -0.00058];
        siderealRotRate = 0.41354; % days
    case 'saturn'
        Re      = 60268.e3;   % m
        GM      = 3.794e16;   % m^3/s^2
        Jvalues = [0.01645 0 -0.001];
        siderealRotRate = 0.4375; % days
    case 'uranus'
        Re      = 25559.e3;   % m
        GM      = 5.794e15;   % m^3/s^2
        Jvalues = 0.012;
        siderealRotRate = -0.65; % days
    case 'neptune'
        Re      = 24764e3;    % m
        GM      = 6.809e15;   % m^3/s^2
        Jvalues = 0.004;
        siderealRotRate = 0.768; % days
%     case 'pluto'
%         Re      = 1151e3; % m
%         GM      = 9e11;   % m^3/s^2
%         Jvalues = [];
%         siderealRotRate = -6.3867; % days
    otherwise
        % case 'earth'
        % JGM-2 gravity model
        Re      = 6378.1363e3;    % m
        GM      = 3.986004415e14; % m^3/s^2
        Jvalues = [ 0.0010826269 -0.0000025323 -0.0000016204 ];
        siderealRotRate = 0.99726968; % days
end

omega   = 2*pi/siderealRotRate/24/3600; % rad/sec
