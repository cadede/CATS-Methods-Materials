function [T, a, P, rho] = atmoscoesa( h, varargin )
%  ATMOSCOESA Use 1976 COESA atmosphere.
%   [T, A, P, RHO] = ATMOSCOESA(H, ACTION) implements the mathematical
%   representation of the 1976 Committee on Extension to the Standard
%   Atmosphere (COESA) United States standard lower atmospheric values for
%   absolute temperature, pressure, density, and speed of sound for the
%   input geopotential altitude. 
%
%   Below 32000 meters (approximately 104987 feet), the U.S. Standard
%   Atmosphere is identical with the Standard Atmosphere of the
%   International Civil Aviation Organization (ICAO). 
%
%   Inputs for ATMOSCOESA is:
%   H      :an array of M geopotential height in meters. 
%   ACTION :a string to determine action for out of range input. Specify if
%          out of range input invokes a 'Warning', 'Error', or no action
%          ('None'). The default is 'Warning'.
%
%   Outputs calculated for the COESA model are: 
%   T      :an array of M temperature in kelvin.
%   a      :an array of M speed of sound in meters per second.
%   P      :an array of M air pressure in pascal.
%   rho    :an array of M air density in kilograms per meter cubed.
%
%   Limitation:
%
%   This function has the limitations of the COESA model. For more
%   information see the documentation. 
%
%   Examples:
%
%   Calculate the COESA model at 1000 meters with warnings for out of range
%   inputs: 
%      [T, a, P, rho] = atmoscoesa(1000)
%
%   Calculate the COESA model at 1000, 11000 and 20000
%   meters with errors for out of range inputs:
%      [T, a, P, rho] = atmoscoesa([1000 11000 20000], 'Error' )
%
%   See also ATMOSCIRA, ATMOSISA, ATMOSLAPSE, ATMOSNONSTD, ATMOSNRLMSISE00, ATMOSPALT.

%   Copyright 2000-2016 The MathWorks, Inc.

%   Limitation: Below the geopotential altitude of 0 m (0 feet) and above
%   the geopotential altitude of 84852 m (approximately 278386 feet),
%   temperature values are extrapolated linearly and pressure values are
%   extrapolated logarithmically. Density and speed of sound are calculated
%   using a perfect gas relationship.  
%
%   Reference:  U.S. Standard Atmosphere, 1976, U.S. Government Printing 
%   Office, Washington, D.C.

narginchk(1, 2);

if ~isnumeric( h )
    % Altitude should be a numeric array.  Otherwise error.
    error(message('aero:atmoscoesa:notNumeric'));
end

action = 'warning';

if nargin > 1
    if ~ischar( varargin{1} ) && ~isstring( varargin{1} )
        error(message('aero:atmoscoesa:inputType'));
    end
    varargin{1} = char(varargin{1});
    if any(strcmpi({'error' 'warning' 'none'},varargin{1}))
        action = lower( varargin{1} );
    else
        error(message('aero:atmoscoesa:unknownString',varargin{ 1 }));
   end
end

% Handle action for height checking
handleaction();

[T, a, P, rho] = coesamethod( h );

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function handleaction()
        switch action
            case 'none'
                % no checking of height
            case 'warning'
                if (any( h(:) < 0.0 ) || any( h(:) > 84528.0 ))
                    warning(message('aero:atmoscoesa:tooLowWarn'));
                end
            case 'error'
                if (any( h(:) < 0.0 ) || any( h(:) > 84528.0 ))
                    error(message('aero:atmoscoesa:tooLowError'));
                end
            otherwise
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end
