function [mach, nu, mu] = flowprandtlmeyer(gamma, varargin)
%FLOWPRANDTLMEYER   Calculate Prandtl-Meyer functions for expansion waves
%   [MACH, NU, MU] = FLOWPRANDTLMEYER(GAMMA, VAR, MTYPE) computes an array
%   of Mach numbers, MACH, Prandtl-Meyer angles, NU in degrees, and Mach
%   angles, MU in degrees. FLOWPRANDTLMEYER calculates these arrays for a
%   given set of specific heat ratios, GAMMA, and any one of the
%   Prandtl-Meyer variables.  The Prandtl-Meyer flow variable is selected
%   by the string, MTYPE.
%
%   See documentation for more information about assumptions and
%   limitations.
%
%   Inputs for FLOWPRANDTLMEYER are:
%
%   GAMMA   :  An array of N specific heat ratios.  GAMMA must be a scalar
%              or an array of N real numbers greater than 1.  For subsonic
%              area ratio input mode and supersonic area ratio input mode,
%              GAMMA must be a real, finite scalar greater than 1.
%
%   VAR     :  An array of real numerical values for one of the
%              Prandtl-Meyer variables.
%
%              MACH :  An array of Mach numbers.  MACH must be a scalar or
%                      an array of N of real numbers greater than or equal
%                      to 0.  If MACH and GAMMA are both arrays, they must
%                      be the same size.  MACH is used with the MTYPE
%                      variable 'mach'.  Since 'mach' is the default MTYPE
%                      value, the MTYPE variable is optional for the Mach
%                      number input mode.
%
%              NU   :  A scalar value for Prandtl-Meyer angle in degrees.
%                      The Prandtl-Meyer angle is the angle change required
%                      for a Mach 1 flow to achieve a given Mach number
%                      after expansion.  NU must be a real scalar greater
%                      than or equal to 0 (at MACH = 1) and less than or
%                      equal to 90 * (SQRT((GAMMA+1)/(GAMMA-1)) - 1) (as
%                      MACH -> inf).  NU is used with MTYPE variable 'nu'.
%
%              MU   :  An array of Mach angles in degrees.  The Mach angle
%                      is the angle between the flow direction and the
%                      lines of pressure disturbance caused by supersonic
%                      motion.  The Mach angle is a function of Mach number
%                      only.  MU must be a scalar or an array of N real
%                      numbers greater than or equal to 0 (as MACH -> inf)
%                      and less than or equal to 90 (at MACH = 1).  MU is
%                      used with MTYPE variable 'mu'.
%
%   MTYPE  :  A character vector for selecting the isentropic flow variable
%             represented by VAR.
%
%              'mach'  :  Default value.  Indicates that the function is in
%                         Mach number input mode.
%
%              'nu'    :  Indicates that the function is in Prandtl-Meyer
%                         angle input mode.
%
%              'mu'    :  Indicates that the function is in Mach angle
%                         input mode.
%
%   Outputs calculated for FLOWPRANDTLMEYER are:
%   (All outputs are the same size as the array input or array inputs.
%   If there are no array inputs, all outputs are scalars.)
%
%   MACH    :  An array of Mach numbers.  In Prandtl-Meyer angle input
%              mode, MACH outputs are the same size as the array input or
%              array inputs. If there are no array inputs then MACH is a
%              scalar.
%
%   NU      :  An array of Prandtl-Meyer angles.  The Prandtl-Meyer angle
%              is the angle change required for a Mach 1 flow to achieve a
%              given Mach number after expansion.
%
%   MU      :  An array of Mach angles.  The Mach angle is between the flow
%              direction and the lines of pressure disturbance caused by
%              supersonic motion.
%
%   Examples:
%
%   Calculate the Prandtl-Meyer relations for air (gamma = 1.4) for 
%   Prandtl-Meyer angle 61 degrees.  The following returns a scalar for
%   MACH, NU, and MU.
%
%      [MACH, NU, MU] = flowprandtlmeyer(1.4, 61, 'nu')
%
%   Calculate the Prandtl-Meyer functions for gases with specific heat
%   ratios given following yields a 1 x 4 array for NU, but only a scalar
%   for MACH and MU.  Since there is only one Mach number input and MU
%   depends on MACH only, this makes sense.
%
%      gamma = [1.3, 1.33, 1.4, 1.67];
%      [MACH, NU, MU] = flowprandtlmeyer(gamma, 1.5)
%
%   Calculate the Prandtl-Meyer angles for a specific heat ratio of 1.4 and
%   range of Mach angles from 40 degrees to 70 degrees in increments of 10
%   degrees.  The following returns a 4 x 1 column array for MACH, NU, and
%   MU.
%
%      [MACH, NU, MU] = flowprandtlmeyer(1.4, (40:10:70)', 'mu')
%      
%   Calculate the Prandtl-Meyer relations for gases with specific heat
%   ratio and Mach number combinations as shown.  The following returns a
%   1 x 2 array for NU and MU each, where the elements of each vector
%   correspond to the inputs element-wise.
%
%       gamma = [1.3, 1.4]
%       mach = [1.13, 9]
%       [MACH, NU, MU] = flowprandtlmeyer(gamma,mach)
%
%   See also FLOWISENTROPIC, FLOWNORMALSHOCK, FLOWFANNO, FLOWRAYLEIGH.


%   Copyright 2009-2016 The MathWorks, Inc.

%   References:  James, J. E. A., Gas Dynamics, Second Edition, Allyn and 
%   Bacon, Inc, Boston, 1984.
%   
%   NACA Technical Report 1135, 1953, National Advisory Committee on
%   Aeronautics, Ames Research Staff, Moffett Field, Calif.

% Checks

% Error in case of no input or when there are more than three inputs
narginchk(2,3);

% Error if third input exists and is not an acceptable string
if nargin == 3
    if ~ischar(varargin{2}) && ~isstring(varargin{2})
        error(message('aero:flowprandtlmeyer:paramSelectString'));
    end
end

% Check gamma

% Error if specific heat ratio input is not numeric
if ~isnumeric(gamma)
    error(message('aero:flowprandtlmeyer:notNumericGamma'));
end

% gamma > 1 check (specific heat ratio must be greater than one)
if any(gamma<=1)
    error(message('aero:flowprandtlmeyer:gammaOneOrLess'));
end

% gamma real number check
if ~isreal(gamma)
    error(message('aero:flowprandtlmeyer:imaginaryGamma'));
end

% Error if gamma input is a matrix or if second input variable is a matrix
if (~isvector(gamma) || ~isvector(varargin{1}))
    error(message('aero:flowprandtlmeyer:multiDimensional'))
end

% Second input general checks (varargin{1})

% Error if second input variable is not numeric
if ~isnumeric(varargin{1})
    error(message('aero:flowprandtlmeyer:notNumeric'));
end

% Error if inputs are not the same size and if neither are scalars
if ~((isscalar(gamma) || isscalar(varargin{1})) || ...
        (isequal(size(gamma),size(varargin{1}))))
    error(message('aero:flowprandtlmeyer:size'))
end

% Second input variable real number check
if ~isreal(varargin{1})
    error(message('aero:flowprandtlmeyer:imaginary'));
end

% Initialize the string variable as the default (mach)
param = 'mach';

if nargin == 3 % If there is a selector string for the second input
    
    % Making the string inputs more flexible
    
    % Typing the string 'mach' to input Mach number is optional
    if strcmpi(varargin{2},'mach')
        param = 'mach';
        
    % All the user really has to type in is 'nu' to input Prandtl-Meyer
    % angle
    elseif strncmpi(varargin{2}, 'nu', 2)
        param = 'nu';
    
    % All the user really has to type in is 'mu' to input Mach angle
    elseif strncmpi(varargin{2}, 'mu', 2)
        param = 'mu';
        
    % If the user has a third variable and does not use the third variable
    % as an acceptable string then provide an error message
    else
        error(message('aero:flowprandtlmeyer:paramSelectWrongInput'))
    end
end

% Specific cases

switch lower(param)
    
    case 'mach'
        
        % Error if upstream Mach number input is less than unity.
        if any(varargin{1}<1)
            error(message('aero:flowprandtlmeyer:lessThanUnity'));
        end
        
        % The second input (first variable input) is Mach number
        mach = varargin{1};
        
    case 'nu'
        
        % Save the user's Prandtl-Meyer angle
        nu_in = varargin{1};
                
        % Error if Prandtl-Meyer angle or specific heat ratio inputs are
        % not scalars for Prandtl-Meyer angle input mode
        if ~(isscalar(gamma) && isscalar(varargin{1}))
            error(message('aero:flowprandtlmeyer:prandtlMeyerAngleScalar'))
        end

        upperLimit = atand(inf) * (sqrt((gamma+1)/(gamma-1)) - 1);
        
        obLoTurn = upperLimit(varargin{1}<0);
        obLoGamma = gamma(varargin{1}<0);
        obHiTurn = upperLimit(varargin{1}>upperLimit);
        obHiGamma = gamma(varargin{1}>upperLimit);

        % Error if the Prandtl-Meyer angle is less than 0
        if varargin{1}<0
            error(message('aero:flowprandtlmeyer:prandtlMeyerAngle', sprintf( '%.17f', obLoTurn ), sprintf( '%.15f', obLoGamma )))
        end
        
        % Error if the Prandtl-Meyer angle is greater than
        % 90*(SQRT((GAMMA+1)/(GAMMA-1))
        if varargin{1}>upperLimit
            error(message('aero:flowprandtlmeyer:prandtlMeyerAngle', sprintf( '%.17f', obHiTurn ), sprintf( '%.15f', obHiGamma )))
        end
        
        % Error if either input is not finite in Prandtl-Meyer angle input
        % mode
        if ~(isfinite(gamma) && isfinite(varargin{1}))
            error(message('aero:flowprandtlmeyer:prandtlMeyerAngleFinite'))
        end
        
        % The second input (first variable input) is Prandtl-Meyer angle
        nu = varargin{1};
        
        % Find Mach number through this function (below)
        if any(nu==upperLimit)
            mach(nu==upperLimit) = inf;
        else
            mach = invertPrandtlMeyerAngle(gamma, nu);
        end
        
        
    case 'mu'
        
        % Error if the Mach angle is less than 0
        if varargin{1}<0
            error(message('aero:flowprandtlmeyer:machAngle'))
        end
        
        % Error if the Mach angle is greater than 90 degrees
        if varargin{1}>90
            error(message('aero:flowprandtlmeyer:machAngle'))
        end
        
        % The second input (first variable input) is Mach angle
        mu = varargin{1};
        
        % Find Mach number from Mach angle with this closed form equation.
        mach = 1./sind(mu);
        
end

% Function Body

% Prandtl-Meyer angle, nu
nu=((gamma+1)./(gamma-1)).^(1/2).*atand(((gamma-1)./(gamma+1).*(mach.^2 ...
        -1)).^(1/2))-atand((mach.^2-1).^(1/2));
    
% Mach angle, mu
mu=asind(1./mach);

% For inputs other than mach, give the user back exactly what they put in.
switch param
    case 'mach'
        mach = varargin{1}.*ones(size(gamma));
        mu = mu.*ones(size(gamma));
    case 'nu'
        nu = nu_in.*ones(size(mach));
    case 'mu'
        mu = varargin{1}.*ones(size(gamma));
        mach = mach.*ones(size(gamma));
end

%--------------------------------------------------------------------------
function mach = invertPrandtlMeyerAngle(gamma, nu)
% This function uses a residual method of solving for the Mach number when
% given Prandtl-Meyer angle because a closed form solution for Mach number
% of the Prandtl-Meyer angle equation cannot be found.  This method
% involves solving the following residual equation for zeros: 
% fcn(prandtlMeyerAngle) = (known) value of prandtlMeyerAngle - 
%                      equation for prandtlMeyerAngle as a fcn(GAMMA,MACH)
%
% The Mach number is being sought so it must remain anonymous until its
% value has been found by the fzero function.  The quality of the solution
% is critically dependent on the quality of the initial guess in the fzero
% function.

% This equation is fcn(prandtlMeyerAngle) as described above
equation = @(gamma,nu,mach) nu-(((gamma+1)./(gamma-1)).^(1/2).*...
    atand(((gamma-1)./(gamma+1).*(mach.^2-1)).^(1/2))-atand((mach.^2-1).^(1/2)));

% Initial guess of Mach number
guessMach = [1 realmax];

machEval = fzero(@(mach) equation(gamma,nu,mach), guessMach);

% Let the Mach number to be used be the same as evaluated Mach number
mach = machEval;
