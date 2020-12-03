function [gx, gy, gz] =  gravitycentrifugal( p, varargin )
%  GRAVITYCENTRIFUGAL Implement a centrifugal effect of planetary gravity.
%   [GX, GY, GZ] = GRAVITYCENTRIFUGAL( P ) implements the mathematical
%   representation of centrifugal effect for planetary gravity based on
%   planetary rotation rate. Using P, a M-by-3 array of Planet-Centered
%   Planet-Fixed coordinates, GX, GY and GZ, arrays of M gravity values in
%   the x-axis, y-axis and z-axis of the Planet-Centered Planet-Fixed
%   coordinates are calculated for planet. 
%
%   Alternate formats for calling gravity centrifugal effect are:
%   [GX, GY, GZ] = GRAVITYCENTRIFUGAL( P, MODEL )   
%   [GX, GY, GZ] = GRAVITYCENTRIFUGAL( P, 'Custom', OMEGA )   
%
%   Inputs for gravity centrifugal effect are:
%   P      :a M-by-3 array of Planet-Centered Planet-Fixed coordinates in
%           meters where the z-axis is positive towards the north pole. For
%           Earth this would be ECEF coordinates.
%   MODEL  :a string specifying the planetary model:
%          'Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn',
%          'Uranus', 'Neptune' or 'Custom'.  The default is 'Earth'.
%   OMEGA  :a scalar value specifying the rotational rate for the planet in
%          radians per second. This option is only available with MODEL
%          'Custom'. 
%
%   Output calculated for the gravity centrifugal effect includes:
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
%   Examples:
%
%   Calculate the centrifugal effect of Earth gravity in the x-axis at the
%   equator on the surface of Earth:
%       gx = gravitycentrifugal( [-6378.1363e3 0 0] ) 
%
%   Calculate the centrifugal effect of Mars gravity at 15000 meters over
%   the equator and 11000 meters over the north pole:
%       p  = [2412.648e3 -2412.648e3 0; 0 0 3376.2e3]
%       [gx, gy, gz] = gravitycentrifugal( p, 'Mars' )   
%
%   Calculate the precessing centrifugal effect of gravity for Earth at
%   15000 meters over the equator and 11000 meters over the north pole
%   using a custom planetary model at Julian date 2451545:   
%       p       = [2412.648e3 -2412.648e3 0; 0 0 3376e3]
%       % Set Julian date to January 1, 2000 at noon GMT
%       JD      = 2451545
%       % Calculate precession rate in right ascension in meters
%       pres_RA = 7.086e-12 + 4.3e-15*(JD - 2451545)/36525
%       % Calculate the rotational rate in a precessing reference frame
%       Omega   = 7.2921151467e-5 + pres_RA
%       [gx, gy, gz] = gravitycentrifugal( p, 'custom', Omega )
%
%   See also GRAVITYWGS84, GRAVITYZONAL, GRAVITYSPHERICALHARMONIC

%   Copyright 2009-2018 The MathWorks, Inc.

%   References:  
%   [1] Vallado, D. A., "Fundamentals of Astrodynamics and Applications",
%       McGraw-Hill, New York, 1997.  
%
%   [2] NIMA TR8350.2: "Department of Defense World Geodetic System
%   1984, Its Definition and Relationship with Local Geodetic Systems."

narginchk(1, 3);

checkinputs();

% set default values
model = 'earth';

switch nargin
    case 2
        % set model
            if ~ischar( varargin{1} ) && ~isstring( varargin{1} )
                error(message('aero:gravitycentrifugal:inputTypeModel'));
            else
                if strcmpi( varargin{1}, 'custom')
                    narginchk(3, 3);
                else
                    % assign model
                    definemodel( varargin{1} );
                end
            end
    case 3
        if (~ischar( varargin{2} ) && ~isstring( varargin{2} ) && ...
                isreal( varargin{2} ))
            % Set rotation rate
            if (~ischar( varargin{1} ) && ~isstring( varargin{1} )) || ...
              ~(isnumeric( varargin{2} ) && isreal( varargin{2} ))
                error(message('aero:gravitycentrifugal:inputTypeCustom'));
            end
            % set model, and rotation rate
            if ~strcmpi( varargin{1}, 'custom' )
                error(message('aero:gravitycentrifugal:noInputOmega'));
            end
            model = lower( varargin{1} ); % can only be custom at this point
            omega = varargin{2};
        else
            error(message('aero:gravitycentrifugal:inputTypeOmega'));
        end
    otherwise
        % default model
end

switch lower( model )
    case 'custom'
        % rotation rate already defined
    otherwise
        [ ~, ~, ~, omega] = zonalplanetparams( model );
end

gx = omega .* omega .* p(:,1);
gy = omega .* omega .* p(:,2);
gz = zeros(size(p(:,3)));

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkinputs( )
        if ~isnumeric( p )
            error(message('aero:gravitycentrifugal:notNumeric'));
        end
        
        if (size( p, 2) ~= 3)
            error(message('aero:gravitycentrifugal:wrongDimension'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function definemodel(str)
    str = char(str);
        if any(strcmpi({'mercury', 'venus', 'earth', 'moon', 'mars',  ...
                  'jupiter', 'saturn', 'uranus', 'neptune' 'custom'}, str))
            model = lower( str );
        else
            error(message('aero:gravitycentrifugal:unknownModel', str));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end
