function h = atmospalt( P, varargin )
%  ATMOSPALT Calculate pressure altitude based on ambient pressure.
%   H = ATMOSPALT(P, ACTION) computes the pressure altitude based on
%   ambient pressure. Pressure altitude is the altitude with specified
%   ambient pressure in the 1976 Committee on the Extension of the Standard
%   Atmosphere (COESA) United States standard. Pressure altitude is also
%   known as the mean sea level (MSL) altitude. 
%
%   Inputs for ATMOSPALT is:
%   P      :an array of M ambient pressures in pascal.
%   ACTION :a string to determine action for out of range input. Specify if
%          out of range input invokes a 'Warning', 'Error', or no action
%          ('None'). The default is 'Warning'.
%
%   Outputs calculated for the COESA model are: 
%   H      :an array of M pressure altitudes or MSL altitudes in meters. 
%
%   Limitation: 
%
%   Below the pressure of 0.3961 Pa (approximately 0.00006 psi) and above
%   the pressure of 101325 Pa (approximately 14.7 psi), altitude values are
%   extrapolated logarithmically. Air is assumed to be dry and an ideal gas.
%
%   Examples:
%
%   Calculate the pressure altitude at a static pressure of 101325 pascals
%   with warnings for out of range inputs: 
%      h = atmospalt(101325)
%
%   Calculate the pressure altitude at a static pressure of 101325 and
%   26436 pascals with errors for out of range inputs: 
%      h = atmospalt([101325 26436], 'Error' )
%
%   See also ATMOSCIRA, ATMOSCOESA

%   Copyright 2000-2016 The MathWorks, Inc.

%   Reference:  U.S. Standard Atmosphere, 1976, U.S. Government Printing 
%   Office, Washington, D.C.

narginchk(1, 2);

if ~isnumeric( P )
    % Altitude should be a numeric array.  Otherwise error.
    error(message('aero:atmospalt:notNumeric'));
end

action = 'warning';

if nargin > 1
    if ~ischar( varargin{1} ) && ~isstring( varargin{1} )
        error(message('aero:atmospalt:inputType'));
    end
    varargin{1} = char(varargin{1});
    if any(strcmpi({'error' 'warning' 'none'},varargin{1}))
        action = lower( varargin{1} );
    else
        error(message('aero:atmospalt:unknownString', varargin{ 1 }));
    end
end

% Handle action for pressure altitude checking
handleaction();

h = paltmethod( P );

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function handleaction()
        switch action
            case 'none'
                % no checking of pressure
            case 'warning'
                if (any( P(:) < 0.3961 ) || any( P(:) > 101325 ))
                    warning(message('aero:atmospalt:tooLowWarn'));
                end
            case 'error'
                if (any( P(:) < 0.3961 ) || any( P(:) > 101325 ))
                    error(message('aero:atmospalt:tooLowError'));
                end
            otherwise
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end
