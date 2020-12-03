function [gx, gy, gz] =  gravityzonal( p, varargin )
%  GRAVITYZONAL Implement a zonal harmonic representation of planetary gravity.
%   [GX, GY, GZ] = GRAVITYZONAL( P ) implements the mathematical representation
%   of zonal harmonic planetary gravity based on planetary gravitational
%   potential. Using P, a M-by-3 array of Planet-Centered Planet-Fixed
%   coordinates, GX, GY and GZ, arrays of M gravity values in the x-axis,
%   y-axis and z-axis of the Planet-Centered Planet-Fixed coordinates
%   are calculated for planet using fourth order zonal coefficients by
%   default. 
%
%   Alternate formats for calling zonal harmonic gravity are:
%   [GX, GY, GZ] = GRAVITYZONAL( P, DEGREE )   
%   [GX, GY, GZ] = GRAVITYZONAL( P, MODEL )   
%   [GX, GY, GZ] = GRAVITYZONAL( P, MODEL, DEGREE )   
%   [GX, GY, GZ] = GRAVITYZONAL( P, MODEL, DEGREE, ACTION )   
%   [GX, GY, GZ] = GRAVITYZONAL( P, 'Custom', RE, GM, J, ACTION )   
%
%   Inputs for zonal harmonic gravity are:
%   P      :a M-by-3 array of Planet-Centered Planet-Fixed coordinates in
%           meters where the z-axis is positive towards the North Pole. For
%           Earth this would be ECEF coordinates.
%   MODEL  :a string specifying the planetary model:
%          'Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn',
%          'Uranus', 'Neptune' or 'Custom'.  The default is 'Earth'.
%   DEGREE :a scalar value specifying the degree of the harmonic gravity model.
%          The default maximum degree is two for 'Mercury', 'Venus', 'Moon',
%          'Uranus' and 'Neptune'. The default maximum degree is three for
%          'Mars'. The maximum degree is four for 'Earth', 'Jupiter',
%          'Saturn' and 'Custom'.          
%   RE     :a scalar value specifying the planetary equatorial radius in
%          meters. This option is only available with MODEL 'Custom'.
%   GM     :a scalar value specifying the planetary gravitational parameter
%          in meters cubed per second squared.  This option is only
%          available with MODEL 'Custom'. 
%   J      :an array specifying the zonal harmonic coefficients used to
%          calculate zonal harmonics planetary gravity. This option is only
%          available with MODEL 'Custom'. 
%   ACTION :a string to determine action for out of range input. Specify if
%          out of range input invokes a 'Warning', 'Error', or no action
%          ('None'). The default is 'Warning'.
%
%   Output calculated for the zonal harmonic gravity includes:
%   GX     :an array of M gravity values in the x-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared.
%   GY     :an array of M gravity values in the y-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared. 
%   GZ     :an array of M gravity values in the z-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared. 
%
%   Limitations:
%
%   This function has the limitations of excluding the centrifugal effects
%   of planetary rotation, and the effects of a precessing reference frame.
%
%   This function will have a non-zero z-axis term at the equator for third
%   and fourth degree models due to the values of Legendre functions.
%
%   Examples:
%
%   Calculate the gravity in the x-axis at the equator on the surface of
%   Earth, using the fourth degree model with warning actions:
%       gx = gravityzonal( [-6378.1363e3 0 0] ) 
%
%   Calculate the gravity at 100 meters over the south pole of Earth with
%   error actions: 
%       [gx, gy, gz] = gravityzonal( [0 0 -6356.851e3], 'Error' )   
%
%   Calculate the gravity at 15000 meters over the equator and 11000 meters
%   over the north pole using a second order Mars model with warning
%   actions:
%       p  = [2412.648e3 -2412.648e3 0; 0 0 3376.2e3]
%       [gx, gy, gz] = gravityzonal( p, 'Mars', 2, 'Warning' )   
%
%   Calculate the gravity at 15000 meters over the equator and 11000 meters
%   over the north pole using a custom planetary model with no actions:  
%       p       = [2412.648e3 -2412.648e3 0; 0 0 3376e3]
%       GM      = 42828.371901e9  % m^3/s^2
%       Re      = 3397e3          % m
%       Jvalues = [1.95545367944545e-3 3.14498094262035e-5 -1.53773961526397e-5]
%       [gx, gy, gz] = gravityzonal( p, 'custom', Re, GM, Jvalues, 'None' )
%
%   See also GRAVITYWGS84, GRAVITYCENTRIFUGAL, GRAVITYSPHERICALHARMONIC

%   Copyright 2009-2018 The MathWorks, Inc.

%   References:  
%   [1] Vallado, D. A., "Fundamentals of Astrodynamics and Applications",
%       McGraw-Hill, New York, 1997.  
%   [2] Fortescue, P., J. Stark, G. Swinerd, (Eds.). "Spacecraft Systems
%       Engineering, Third Edition", Wiley & Sons, West Sussex, 2003.  
%   [3] Tewari, A., "Atmospheric and Space Flight Dynamics Modeling and
%       Simulation with MATLAB and Simulink", Birkhauser, Boston, 2007.  


narginchk(1, 6);

checkinputs();

% set default values
model = 'earth';
action = 'warning';

switch nargin
    case 2
        % set degree or model
        if ~(isnumeric( varargin{1} ) && isreal( varargin{1} ))
            if ~ischar( varargin{1} ) && ~isstring( varargin{1} ) 
                error(message('aero:gravityzonal:inputTypeVar3'));
            else
                if strcmpi( varargin{1}, 'custom')
                    narginchk(6, 6);
                else
                    % assign model or action
                    modeloraction( varargin{1} );
                end
            end
        else
            model = 'earth';
            maxdeg = varargin{1};
        end
    case 3
        if ~ischar( varargin{2} ) && ~isstring( varargin{2} )
            % not setting action
            if (~ischar( varargin{1} ) && ~isstring( varargin{1} )) || ...
              ~(isnumeric( varargin{2} ) && isreal( varargin{2} ))
                error(message('aero:gravityzonal:inputTypeNoAction4'));
            end
            if strcmpi( varargin{1}, 'custom')
                narginchk(6, 6);
            end
            % set model and degree
            model = varargin{1};
            maxdeg = varargin{2};    
        else
            if ~(isnumeric( varargin{1} ) && isreal( varargin{1} ))
                if ~ischar( varargin{1} ) && ~isstring( varargin{1} )
                    error(message('aero:gravityzonal:inputTypeVar4'));
                end
                if strcmpi( varargin{1}, 'custom')
                    narginchk(6, 6);
                end
                % set model and action
                model = varargin{1};
                defineaction( varargin{2} );              
            else
            % set degree and action
            maxdeg = varargin{1};
            defineaction( varargin{2} );    
            end
        end
    case 4
        % set model, degree and action
        if (ischar( varargin{1} ) || isstring( varargin{1} )) && ...
                strcmpi( varargin{1}, 'custom')
            narginchk(6, 6);
        end
        if (~ischar( varargin{1} ) && ~isstring( varargin{1} )) || ...
          ~(isnumeric( varargin{2} ) && isreal( varargin{2} )) || ...
           (~ischar( varargin{3} ) && ~isstring( varargin{3} ))
            error(message('aero:gravityzonal:inputTypeVar5'));
        end
        model = varargin{1};
        maxdeg = varargin{2};
        defineaction( varargin{3} );
    case 5
        if (~ischar( varargin{1} ) && ~isstring( varargin{1} )) || ...
           ~(isnumeric( varargin{2} ) && isreal( varargin{2} )) || ...
           ~(isnumeric( varargin{3} ) && isreal( varargin{3} )) || ...
           ~(isnumeric( varargin{4} ) && isreal( varargin{4} ))
            error(message('aero:gravityzonal:inputTypeVar6'));
        end
      % set model, Re, GM, J values and degree
        model = varargin{1};
        maxdeg = length(varargin{4}) + 1;
    case 6
        if (~ischar( varargin{1} ) && ~isstring( varargin{1} )) || ...
           ~(isnumeric( varargin{2} ) && isreal( varargin{2} )) || ...
           ~(isnumeric( varargin{3} ) && isreal( varargin{3} )) || ...
           ~(isnumeric( varargin{4} ) && isreal( varargin{4} )) || ...
           (~ischar( varargin{5} ) && ~isstring( varargin{5} ))
            error(message('aero:gravityzonal:inputTypeVar7'));
        end
     % set model, Re, GM, J values, degree and action
        model = varargin{1};
        maxdeg = length(varargin{4}) + 1;
        defineaction( varargin{5} );
    otherwise
        % default degree, model, and action
        % JGM-2 gravity model
        maxdeg = 4;
end

model = char(model);
switch lower( model )
    case 'custom'
        Re      = varargin{2};
        GM      = varargin{3};
        Jvalues = varargin{4};
        if exist('maxdeg','var') && maxdeg > 4
            switch action
                case 'none'
                    % no passing of maxlegendre messages
                case 'warning'
                    warning(message('aero:gravityzonal:exceedMaxLegendreWarn'));
                case 'error'
                    error(message('aero:gravityzonal:exceedMaxLegendre'));
                otherwise
                    error(message('aero:gravityzonal:unknownActionMaxLegendre'));
            end
            maxdeg = 4;
        end
    otherwise
        [ Re, GM, Jvalues ] = zonalplanetparams( model );
        checkmaxdeg( length(Jvalues)+1 );
end

r = sqrt( sum( p.^2, 2 ));
rlat = asin( p(:,3)./ r );

% preallocate Legendre functions
P_legd = zeros(maxdeg+1,2,length(rlat));

crlat = cos(rlat);
srlat = sin(rlat);
 
% calculate Legendre functions required
P_legd(1,1,:) = ones(1,length(rlat));
P_legd(2,1,:) = srlat;
P_legd(2,2,:) = crlat;
P_legd(3,1,:) = 0.5*(3*srlat.^2 - 1);
P_legd(3,2,:) = 3*srlat.*crlat;

if maxdeg > 2
    P_legd(4,1,:) = 0.5*(5*srlat.^3 - 3*srlat);
    P_legd(4,2,:) = 0.5*crlat.*(15*srlat.^2 - 3);
    if maxdeg > 3
        P_legd(5,1,:) = 0.125*(35.*srlat.^4 - 30*srlat.*srlat + 3);
        P_legd(5,2,:) = 2.5*crlat.*(7*srlat.^3 - 3*srlat);
    end
end

% find radial and transverse gravity
gradial = -GM./(r.*r);
glat = 0;
for n = 2:maxdeg
    gradial   = gradial + GM./(r.*r).*(Re./r).^n.*(n+1).*reshape(P_legd(n+1,1,:),size(r))*Jvalues(n-1);
    glat = glat - GM./(r.*r).*(Re./r).^n.*reshape(P_legd(n+1,2,:),size(r))*Jvalues(n-1);
end

% gravity in ECEF coordinates
gx = ((1./r).*gradial - (p(:,3)./(r.*sqrt(p(:,1).^2 + p(:,2).^2))).*glat).*p(:,1);
gy = ((1./r).*gradial - (p(:,3)./(r.*sqrt(p(:,1).^2 + p(:,2).^2))).*glat).*p(:,2);
gz = (1./r).*gradial.*p(:,3) + ((sqrt(p(:,1).^2 + p(:,2).^2))./r).*glat;

% special case for poles
gx(( p(:,1) == 0 ) & ( p(:,2) == 0 )) = 0;
gy(( p(:,1) == 0 ) & ( p(:,2) == 0 )) = 0;

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkinputs( )
        if ~isnumeric( p )
            error(message('aero:gravityzonal:notNumeric'));
        end
        
        if (size( p, 2) ~= 3)
            error(message('aero:gravityzonal:wrongDimension'));
        end
    end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkmaxdeg( upperlimit )
        if exist('maxdeg','var') && maxdeg > upperlimit
            switch action
                case 'none'
                    % no passing of maxdeg messages
                case 'warning'
                    warning(message('aero:gravityzonal:exceedMaxDegWarn'));
                case 'error'
                    error(message('aero:gravityzonal:exceedMaxDeg'));
                otherwise
                    error(message('aero:gravityzonal:unknownActionMaxDeg'));
            end
            maxdeg = upperlimit;
        else
            if ~exist('maxdeg','var')
            % maxdeg was not set
            maxdeg = upperlimit;            
            end
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function modeloraction( str )
        str = char(str);
        switch lower( str )
            case { 'mercury', 'venus', 'earth', 'moon', 'mars', 'jupiter', ...
                   'saturn', 'uranus', 'neptune' 'custom' }
                model = lower( str );
            case { 'error', 'warning', 'none' }
                action = lower( str );
            otherwise
                error(message('aero:gravityzonal:unknownString'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function defineaction(str)
        str = char(str);
        if any(strcmpi({'error' 'warning' 'none'}, str))
            action = lower( str );
        else
            error(message('aero:gravityzonal:unknownAction', str));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end
