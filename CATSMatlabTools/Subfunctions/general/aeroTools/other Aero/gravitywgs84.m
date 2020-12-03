function varargout = gravitywgs84( h, lat, varargin )
%  GRAVITYWGS84 Implement the 1984 World Geodetic System (WGS84)
%  representation of Earth's gravity.
%   G  = GRAVITYWGS84( H, LAT ) implements the mathematical representation
%   of the geocentric equipotential ellipsoid of WGS84. Using H, an array
%   of M altitudes in meters, and LAT, an array of M geodetic latitudes in
%   degrees, G,  an array of M gravity values in the direction normal to
%   the Earth's surface at a specific location is calculated using the
%   Taylor Series method by default. Gravity precision is controlled via
%   the METHOD parameter.
%
%   Alternate formats for calling WGS84 are:
%   G  = GRAVITYWGS84( H, LAT, LON, METHOD, [ NOATM, NOCENT, PREC, JD ], ACTION )   
%   GN = GRAVITYWGS84( H, LAT, LON, 'Exact', [ NOATM, NOCENT, PREC, JD ], ACTION )   
%   [GN, GT] = GRAVITYWGS84( H, LAT, LON, 'Exact', [ NOATM, NOCENT, PREC, JD ], ACTION )   
%
%   Inputs for WGS84 are:
%   H      :an array of M altitudes in meters with respect to the WGS84
%           ellipsoid.
%   LAT    :an array of M geodetic latitudes in degrees where north
%          latitude is positive, and south latitude is negative.
%   LON    :an array of M geodetic longitude in degrees where east
%          longitude is positive, west is negative. This input is only
%          available with METHOD 'CloseApprox' or 'Exact'. 
%   METHOD :a string specifying the method to calculate gravity:
%          'TaylorSeries' 'CloseApprox' or 'Exact'.  The default is
%          'TaylorSeries'.
%   NOATM  :a logical value specifying the exclusion of Earth's atmosphere.
%          Set to true for the Earth's gravitational field to exclude the
%          mass of the atmosphere. Set to false for the value for the
%          Earth's gravitational field to include the mass of the
%          atmosphere. This option is only available with METHOD
%          'CloseApprox' or 'Exact'. The default is false.
%   NOCENT :a logical value specifying the removal of centrifugal effects.
%          Set to true to calculate gravity based on pure attraction
%          resulting from the normal gravitational potential. Set to false
%          to calculate gravity including the centrifugal force resulting
%          from the Earth's angular velocity. This option is only available
%          with METHOD 'CloseApprox' or 'Exact'. The default is false.
%   PREC   :a logical value specifying the presence of a precessing
%          reference frame. Set to true for the angular velocity of the
%          Earth to be calculated using the International Astronomical
%          Union (IAU) value of the Earth's angular velocity and the
%          precession rate in right ascension. In order to obtain the
%          precession rate in right ascension, Julian Centuries from Epoch
%          J2000.0 is calculated using the Julian date, JD. Set to false
%          for the angular velocity of the Earth used is the value of the
%          standard Earth rotating at a constant angular velocity. This
%          option is only available with METHOD 'CloseApprox' or 'Exact'.
%          The default is false.
%   JD     :a scalar value specifying Julian date used to calculate Julian
%          Centuries from Epoch J2000.0. This option is only available with
%          METHOD  'CloseApprox' or 'Exact'.
%   ACTION :a string to determine action for out of range input. Specify if
%          out of range input invokes a 'Warning', 'Error', or no action
%          ('None'). The default is 'Warning'.
%
%   Output calculated for the Earth's gravity includes:
%   G      :an array of M gravity values in the direction normal to the
%          Earth's surface at a specific LAT LON location.  A positive
%          value indicates a downward direction.
%   GN     :an array of M total gravity values in the direction normal to
%          the Earth's surface at a specific LAT LON location. A positive
%          value indicates a downward direction.  This option is only
%          available with METHOD 'Exact'. 
%   GT     :an array of M gravity values in the direction tangential to the
%          Earth's surface at a specific LAT LON location. A positive value
%          indicates a northward direction.  This option is only available
%          with METHOD 'Exact'. 
%
%   Limitations:
%
%   This function has the limitations of the 1984 World Geodetic System
%   (WGS84). For more information see the documentation.  
%
%   Examples:
%
%   Calculate the normal gravity at 5000 meters, and 55 degrees latitude
%   using the taylor series approximation method with error actions:
%       g = gravitywgs84( 5000, 55, 'TaylorSeries', 'Error' ) 
%
%   Calculate the normal gravity at 15000 meters, 45 degrees latitude and
%   120 degrees longitude using the close approximation method with
%   warning actions, atmosphere, centrifugal effects, and no precessing:
%       g = gravitywgs84( 15000, 45, 120, 'CloseApprox' )   
%
%   Calculate the normal and tangential gravity at 1000 meters, 0 degrees 
%   latitude and 20 degrees longitude using the exact method with warning
%   actions, atmosphere, centrifugal effects, and no precessing:
%       [gn, gt] = gravitywgs84( 1000, 0, 20, 'Exact' )   
%
%   Calculate the normal and tangential gravity at 1000 meters, 0 degrees 
%   latitude and 20 degrees longitude and 11000 meters, 30 degrees latitude
%   and 50 degrees longitude using the exact method with no actions,
%   atmosphere, centrifugal effects, and no precessing: 
%       h = [1000; 11000]
%       lat = [0; 30]
%       lon = [20; 50]
%       [gn, gt] = gravitywgs84( h, lat, lon, 'Exact', 'None' )
%
%   Calculate the normal gravity at 15000 meters, 45 degrees latitude and
%   120 degrees longitude and 5000 meters, 55 degrees latitude and 100
%   degrees longitude using the close approximation method with 
%   warning actions, atmosphere, no centrifugal effects, and no precessing:
%       h = [15000 5000]
%       lat = [45 55]
%       lon = [120 100]
%       g = gravitywgs84( h, lat, lon, 'CloseApprox', [false true false 0] )   
%
%   Calculate the normal and tangential gravity at 1000 meters, 0 degrees 
%   latitude and 20 degrees longitude using the exact method with warning
%   actions, atmosphere, centrifugal effects, and precessing at Julian date
%   2451545: 
%       [gn, gt] = gravitywgs84( 1000, 0, 20, 'Exact', ...
%                                 [ false, false, true, 2451545 ], 'Warning' )   
%
%   Calculate the normal gravity at 15000 meters, 45 degrees latitude and
%   120 degrees longitude using the close approximation method with error
%   actions, no atmosphere, centrifugal effects, and precessing at Julian
%   date 2451545:
%       g = gravitywgs84( 15000, 45, 120, 'CloseApprox', ...
%                                   [ true false true 2451545 ], 'Error'  )   
%
%   Calculate the total normal gravity at 15000 meters, 45 degrees latitude
%   and 120 degrees longitude using the exact method with error actions, no
%   atmosphere, centrifugal effects, and precessing at Julian date 2451545:
%       gn = gravitywgs84( 15000, 45, 120, 'Exact', ...
%                                   [ true false true 2451545 ], 'Error'  )   
%
%   See also GEOIDEGM96, GRAVITYZONAL, GRAVITYCENTRIFUGAL, GRAVITYSPHERICALHARMONIC

%   Copyright 2000-2018 The MathWorks, Inc.

%   Reference:  NIMA TR8350.2: "Department of Defense World Geodetic System
%   1984, Its Definition and Relationship with Local Geodetic Systems."
 
%   Assumptions and Limitations:  The WGS84 gravity calculations are based
%   on the assumption of a geocentric equipotential ellipsoid of
%   revolution. Since the gravity potential is assumed to be the same
%   everywhere on the ellipsoid, there must be a specific theoretical
%   gravity potential that can be uniquely determined from the four
%   independent constants defining the ellipsoid.  Use of the WGS84 Taylor
%   Series model should be limited to low geodetic heights. It is
%   sufficient near the surface when submicrogal precision is not
%   necessary. At medium and high geodetic heights, it is less accurate.
%   Use of the WGS84 Close Approximation model should be limited to a
%   geodetic height of 20000.0 m (approximately 65620.0 feet). Below this
%   height, it gives results with submicrogal precision.  

narginchk(2, 6);
nargoutchk(0, 2);

checkinputs();

rlat = convang( lat, 'deg', 'rad');

rlon = zeros(size(lat));
method = 'taylorseries';
options = [ false false false 0 ];
action = 'warning';

switch nargin
    case 2
        % Taylor Series calculation
        % inputs are preset
    case 3
         % Taylor Series calculation
         if ~ischar( varargin{1} ) && ~isstring( varargin{1} )
             error(message('aero:gravitywgs84:inputType'));
         end

         % assign method or action
         taylororaction( varargin{1} );
    case 4
         if ~((ischar( varargin{1} ) || isstring( varargin{1} )) && ...
             (ischar( varargin{2} ) || isstring( varargin{2} )))
             if (isnumeric( varargin{1}) && ...
                (ischar( varargin{2} ) || isstring( varargin{2} )))
                 % Close Approximation or Exact calculation
                 methodexactorclose( varargin{1}, varargin{2} );
             else
                 error(message('aero:gravitywgs84:inputOrder4'));
             end
         else
             % Taylor Series calculation
             for k = 1:2
                 % assign method or action
                 taylororaction( varargin{k} );
             end
         end
    case 5
        % Close Approximation or Exact calculation
        if (isnumeric( varargin{1}) && (ischar( varargin{2} ) || isstring( varargin{2} )))
            methodexactorclose( varargin{1}, varargin{2} );
            if ischar( varargin{3} ) || isstring( varargin{3} )
                % Action defined
                defineaction( varargin{3} );
            elseif isnumeric( varargin{3} )
                % Options defined
                options = varargin{3};
            else 
                error(message('aero:gravitywgs84:inputOrder'));
            end
        else
            error(message('aero:gravitywgs84:inputOrder5'));
        end
    case 6
        % Close Approximation or Exact calculation
        if (isnumeric( varargin{1}) && (ischar( varargin{2} ) || isstring( varargin{2} )) ...
                && isnumeric( varargin{3} ) && ...
                (ischar(varargin{4}) || isstring( varargin{4} )))
            methodexactorclose( varargin{1}, varargin{2} );
            % Options defined
            options = varargin{3};
            % Action defined
            defineaction( varargin{4} );
        else
            error(message('aero:gravitywgs84:inputOrder6'));
        end
  otherwise
        % Can't get to this code due to nargin check    
        error(message('aero:gravitywgs84:inputNumber')); 
end

% wrap latitude and longitude if needed
[lat_wrapped, rlat, rlon] = wraplatitude( rlat, rlon, 'rad' );

% check and fix angle wrapping in longitude
[lon_wrapped, rlon] = wraplongitude( rlon, 'rad', 'pi' );

% Handle action for height checking, lat wrapping and lon wrapping
handleaction();

udata = struct('noatmos', options(1), 'precessing', options(3), ...
               'JD', options(4), 'centrifugal', options(2) );
         
% Enumerate the method
type = {'taylorseries' 'closeapprox' 'exact'};           
emethod = find(strcmpi(type,method));

% Call gravity calculations
if (emethod == 3)
    if nargout == 2
        [~, varargout{1}, varargout{2}] = wgs84methods( h, rlat, rlon, emethod, udata );
    else
        [varargout{1}, ~, ~] = wgs84methods( h, rlat, rlon, emethod, udata ); 
    end
else
    if nargout < 2
        varargout{1} = wgs84methods( h, rlat, rlon, emethod, udata );
    else
        error(message('aero:gravitywgs84:numOutputs'));
    end
end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkinputs( )
        if ~isnumeric( h )
            % Altitude should be a numeric array.  Otherwise error.
            error(message('aero:gravitywgs84:notNumericAltitude'));
        end

        if ~isnumeric( lat )
            % Lat should be a numeric array.  Otherwise error.
            error(message('aero:gravitywgs84:notNumericLatitude'));
        end

        if ~all(size(h) == size(lat))
            error(message('aero:gravitywgs84:arraySizeLatitude'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function taylororaction( str )
        str = char(str);
        switch lower( str )
            case 'taylorseries'
                method = lower( str );
            case { 'error', 'warning', 'none' }
                action = lower( str );
            otherwise
                error(message('aero:gravitywgs84:unknownString'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function methodexactorclose( num, str )
        rlon = convang( num, 'deg', 'rad' );
        str = char(str);
        if ~all(size(h) == size(rlon))
            error(message('aero:gravitywgs84:arraySize'));
        end

        if (strcmpi('exact', str) || ...
                strcmpi('closeapprox',str ))
            method = lower( str );
        else
            error(message('aero:gravitywgs84:unknownStringMethod', str));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function defineaction(str)
        str = char(str);
        if any(strcmpi({'error' 'warning' 'none'}, str))
            action = lower( str );
        else
            error(message('aero:gravitywgs84:unknownStringAction', str));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function handleaction()
        switch action
            case 'none'
                % no checking of height nor wrapping messages
            case 'warning'
                if (any(h(:) < 0.0))
                    warning(message('aero:gravitywgs84:tooLowWarn'));
                end
                if ((strcmpi( method,'closeapprox')|| ...
                        strcmpi( method,'taylorseries')) ...
                        && (any(h(:) > 20000.0)))
                    warning(message('aero:gravitywgs84:tooHighWarn'));
                end
                if lat_wrapped
                    warning(message('aero:gravitywgs84:latitudeWrap'));
                    warning(message('aero:gravitywgs84:longitudeAdjust'));
                end
                if lon_wrapped
                    warning(message('aero:gravitywgs84:longitudeWrap'));
                end
            case 'error'
                if (any(h(:) < 0.0))
                    error(message('aero:gravitywgs84:tooLowError'));
                end
                if ((strcmpi( method,'closeapprox')|| ...
                        strcmpi( method,'taylorseries')) ...
                        && (any(h(:) > 20000.0)))
                    error(message('aero:gravitywgs84:tooHighError'));
                end
                if lat_wrapped
                    error(message('aero:gravitywgs84:latitude90'));
                end
                if lon_wrapped
                    error(message('aero:gravitywgs84:longitude180'));
                end
            otherwise
                error(message('aero:gravitywgs84:unknownAction'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end
