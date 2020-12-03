function [mach, T, P, rho, V, P0, F] = flowfanno(gamma, varargin)
%FLOWFANNO      Calculate Fanno line flow relations
%   [MACH, T, P, RHO, V, P0, F] = FLOWFANNO(GAMMA, VAR, MTYPE) returns an
%   array for each of the Fanno line flow relations.  FLOWFANNO calculates
%   these arrays for a given set of specific heat ratios, GAMMA, and any
%   one of the Fanno flow variables.  The Fanno flow variable is selected
%   by the string, MTYPE.  The outputs are the Mach number, MACH, 
%   temperature ratio, T, the pressure ratio, P, the density ratio, RHO,
%   the velocity ratio, V, and the stagnation (total) pressure ratio, P0.
%   All of these ratios are local conditions over the sonic conditions.  F
%   is the Fanno parameter given by F = f*L/D, where f is the friction
%   coefficient, L is the length of constant area duct required to achieve
%   sonic flow, and D is the hydraulic diameter of the duct.
%
%   See documentation for more information about assumptions and
%   limitations.
%
%   Inputs for FLOWFANNO are:
%
%   GAMMA   :  An array of N specific heat ratios.  GAMMA must be a scalar
%              or an array of N real numbers greater than 1.  For subsonic
%              total pressure ratio, supersonic total pressure ratio,
%              subsonic Fanno parameter, and supersonic Fanno parameter
%              input modes, GAMMA must be a real, finite scalar greater
%              than 1.
%
%   VAR     :  An array of real numerical values for one of the Fanno flow
%              variables.
%
%              MACH :  An array of Mach numbers.  Must be a scalar or an
%                      array of N of real numbers greater than or equal to
%                      0.  If MACH and GAMMA are both arrays, they must be
%                      the same size.  MACH is used with MTYPE of 'mach'.
%                      Since 'mach' is the default MTYPE value, the MTYPE
%                      variable is optional for Mach number input mode.
%
%              T    :  An array of temperature ratios.  The temperature
%                      ratio is defined as the local static temperature
%                      over the reference static temperature for sonic
%                      flow.  T must be a scalar or an array of N of real
%                      numbers greater than or equal to 0 (as MACH -> inf)
%                      and less than or equal to (GAMMA+1)/2 (at MACH = 0).
%                      T is used with MTYPE variable 'temp'.
%
%              P    :  An array of pressure ratios.  The pressure ratio is
%                      defined as the local static pressure over the
%                      reference static pressure for sonic flow.  P must be
%                      a scalar or an array of N real numbers greater than
%                      or equal to 0.  If P and GAMMA are both arrays, they
%                      must be the same size.  P is used with MTYPE
%                      variable 'pres'.
%
%              RHO  :  An array of density ratios.  The density ratio is
%                      defined as the local density over the reference
%                      density for sonic flow.  RHO must be a scalar or an
%                      array of N of real numbers greater than or equal to
%                      SQRT((GAMMA-1)./(GAMMA+1)) (as MACH -> inf).  If RHO
%                      and GAMMA are both arrays, they must be the same
%                      size.  RHO is used with MTYPE variable 'dens'.
%
%              V    :  An array of velocity ratios.  The velocity ratio is
%                      defined as the local velocity over the reference
%                      velocity for sonic flow.  V must be a scalar or an
%                      array of N of real numbers greater than or equal to
%                      0 and less than or equal to
%                      SQRT((GAMMA+1)./(GAMMA-1)) (as MACH -> inf).  If V
%                      and GAMMA are both arrays, they must be the same
%                      size.  V is used with MTYPE variable 'velo'.
%
%              P0   :  A scalar value for total pressure ratio.  The total
%                      pressure ratio is defined as the local total
%                      pressure over the reference total pressure for sonic
%                      flow.  P0 must be greater than or equal to 1.  P0 is
%                      used with MTYPE variables 'totalpsub' and
%                      'totalpsup'.
%
%              F    :  A scalar value for Fanno parameter.  The Fanno
%                      parameter is defined as F = f*L/D, where Fanno
%                      parameter F = f*L/D, f is the friction coefficient,
%                      L is the length of constant area duct required to
%                      achieve sonic flow, and D is the hydraulic diameter
%                      of the duct.  In subsonic mode, F must be greater
%                      than or equal to 0.  In supersonic mode, F must be
%                      greater than or equal to 0 (at MACH = 1) and less 
%                      than or equal to
%                      (GAMMA+1)/(2*GAMMA)*log((GAMMA+1)/(GAMMA-1))-1/GAMMA
%                      (as MACH -> inf).  F is used with MTYPE variables
%                      'fannosub' and 'fannosup'.
%
%   MTYPE   :  A character vector for selecting the Fanno flow variable
%              represented by VAR.
%
%              'mach'       :  Default value.  Indicates that the function
%                              is in Mach number input mode.
%
%              'temp'       :  Indicates that the function is in
%                              temperature ratio input mode.
%
%              'pres'       :  Indicates that the function is in pressure
%                              ratio input mode.
%
%              'dens'       :  Indicates that the function is in density
%                              ratio input mode.
%
%              'velo'       :  Indicates that the function is in velocity
%                              ratio input mode.
%
%              'totalpsub'  :  Indicates that the function is in subsonic
%                              total pressure ratio input mode.
%
%              'totalpsup'  :  Indicates that the function is in supersonic
%                              total pressure ratio input mode.
%
%              'fannosub'   :  Indicates that the function is in subsonic
%                              Fanno parameter input mode.
%
%              'fannosup'   :  Indicates that the function is in supersonic
%                              Fanno parameter input mode.
%
%   Outputs calculated for FLOWFANNO are:
%   (All outputs are the same size as the array input or array inputs.
%   If there are no array inputs, all outputs are scalars.)
%
%   MACH    :  An array of Mach numbers.
%
%   T       :  An array of temperature ratios.  The temperature ratio is
%              defined as the local static temperature over the reference
%              static temperature for sonic flow.
%
%   P       :  An array of pressure ratios.  The pressure ratio is defined
%              as the local static pressure over the reference static
%              pressure for sonic flow.
%
%   RHO     :  An array of density ratios.  The density ratio is defined as
%              the local density over the reference density for sonic flow.
%
%   V       :  An array of velocity ratios.  The velocity ratio is defined
%              as the local velocity over the reference velocity for sonic
%              flow.
%
%   P0      :  An array of total pressure ratios.  The total pressure ratio
%              is defined as the local total pressure over the reference
%              total pressure for sonic flow.
%
%   F       :  An array of Fanno parameters.  The Fanno parameter is
%              defined as F = f*L/D, where Fanno parameter F = f*L/D, f is
%              the friction coefficient, L is the length of constant area
%              duct required to achieve sonic flow, and D is the hydraulic
%              diameter of the duct.
%
%   Examples:
%
%   Calculate the Fanno line flow relations for air (gamma = 1.4) for
%   subsonic Fanno parameter 1.2.  The following returns scalar values for
%   MACH, T, P, RHO, V, P0, and F.
%
%      [MACH, T, P, RHO, V, P0, F] = flowfanno(1.4, 1.2, 'fannosub')
%
%   Calculate the Fanno line flow relations for gases with specific heat
%   ratios given in the following 1 x 4 row array for the Mach number 0.5.
%   The following yields a 1 x 4 row array for MACH, T, P, RHO, V, P0, and
%   F.
%
%      gamma = [1.3, 1.33, 1.4, 1.67]
%      [MACH, T, P, RHO, V, P0, F] = flowfanno(gamma, 0.5)
%
%   Calculate the Fanno line flow relations for a specific heat ratio of
%   1.4 and range of temperature ratios from 0.40 to 0.70 in increments of
%   0.10.  The  following returns a 4 x 1 column array for MACH, T, P, RHO,
%   V, P0, and F.
%
%      [MACH, T, P, RHO, V, P0, F] = flowfanno(1.4, [1.1 1.2], 'temp')
%
%   Calculate the Fanno line flow relations for gases with specific heat
%   ratio and velocity ratio combinations as shown.  The following returns
%   a 1 x 2 array for MACH, T, P, RHO, V, P0, and F each, where the
%   elements of each array corresponds to the inputs element-wise.
%
%       gamma = [1.3, 1.4]
%       V = [0.53, 0.49]
%       [MACH, T, P, RHO, V, P0, F] = flowfanno(gamma, V, 'velo')
%
%   See also FLOWISENTROPIC, FLOWNORMALSHOCK, FLOWPRANDTLMEYER,
%   FLOWRAYLEIGH.


%   Copyright 2009-2017 The MathWorks, Inc.

%   References:  James, J. E. A., Gas Dynamics, Second Edition, Allyn and 
%   Bacon, Inc, Boston, 1984.
%
%   NACA Technical Report 1135, 1953, National Advisory Committee on
%   Aeronautics, Ames Research Staff, Moffett Field, Calif.

% Checks

% General checks

% Error in case of less than 2 inputs or when there are more than three inputs
narginchk(2,3);

% Error if third input exists and is not an acceptable string
if nargin > 2
    if ~ischar(varargin{2}) && ~isstring(varargin{2})
        error(message('aero:flowfanno:paramSelectString'))
    end
end

% Check gamma

% Error if specific heat ratio input is not numeric
if ~isnumeric(gamma)
    error(message('aero:flowfanno:notNumericGamma'));
end

% gamma >1 check (specific heat ratio must be greater than one)
if any(gamma<=1)
    error(message('aero:flowfanno:gammaOneOrLess'));
end

% gamma real number check
if ~isreal(gamma)
    error(message('aero:flowfanno:imaginaryGamma'));
end

% Error if gamma input is a matrix or if second input variable input is a matrix
if (~isvector(gamma) || ~isvector(varargin{1}))
    error(message('aero:flowfanno:multiDimensional'))
end

% Second input general checks (varargin{1})

% Error if second input variable is not numeric input is not numeric
if ~isnumeric(varargin{1})
    error(message('aero:flowfanno:notNumeric'));
end

% Error if inputs are not the same size and if neither are scalars
if ~((isscalar(gamma) || isscalar(varargin{1})) || ...
        (isequal(size(gamma),size(varargin{1}))))
    error(message('aero:flowfanno:size'))
end

% Second input variable real number check
if ~isreal(varargin{1})
    error(message('aero:flowfanno:imaginary'));
end

% Initialize the string variable as the default (mach)
param = 'mach';

if nargin > 2 % If there are selector strings for the second input
    
    % Making string inputs more flexible
    
    % Typing the string 'mach' to input Mach number is optional
    if strcmpi(varargin{2},'mach')
        param = 'mach';
    
    % All the user really has to type in is "temp' to input
    % temperatureRatio.
    elseif strncmpi(varargin{2},'tempratio',4)
        param = 'tempratio';
        
    % All the user really has to type in is 'pres' to input pressureRatio.
    elseif strncmpi(varargin{2},'pressratio',4)
        param = 'pressratio';
        
    % All the user really has to type in is 'dens' to input densityRatio.
    elseif strncmpi(varargin{2},'densratio',4)
        param = 'densratio';
        
    % All the user really has to type in is 'velo' to input velocityRatio.
    elseif strncmpi(varargin{2},'veloratio',4)
        param = 'veloratio';
        
    % All the user really has to type in is 'totalpsub' to input
    % subsonicTotalPressureRatio.
    elseif strncmpi(varargin{2},'totalpsub',9)
        param = 'subsonictotalpressureratio';
        
    % All the user really has to type in is 'totalpsup' to input
    % supersonicTotalPressureRatio.
    elseif strncmpi(varargin{2},'totalpsup',9)
        param = 'supersonictotalpressureratio';
        
    % All the user really has to type in is 'fannosub' to input
    % subsonicFannoParameter.
    elseif strncmpi(varargin{2},'fannosub',8)
        param = 'subsonicfannoparameter';
        
    % All the user really has to type in is 'fannosup' to input
    % supersonicFannoParameter.
    elseif strncmpi(varargin{2},'fannosup',8)
        param = 'supersonicfannoparameter';
        
    % If the user has a third variable and does not use the third variable
    % as an acceptable string then provide an error message
    else
        error(message('aero:flowfanno:paramSelectWrongInput'))
    end
end

% Specific cases

switch lower(param)
    
    case 'mach'
        
        % MachNo non-negative check (Mach number must be greater than or equal to 0
        if any(varargin{1}<0)
            error(message('aero:flowfanno:negative'));
        end
        
        % The second input (first variable input) is Mach number.
        mach = varargin{1};
        
    case 'tempratio'
        
        upperLimit = (gamma+1)/2;
        
        if ~isscalar(gamma)
            obLoTemp = upperLimit(varargin{1}<0);
            obLoGamma = gamma(varargin{1}<1);
            obHiTemp = upperLimit(varargin{1}>upperLimit);
            obHiGamma = gamma(varargin{1}>upperLimit);
        else
            obLoTemp = upperLimit;
            obLoGamma = gamma;
            obHiTemp = upperLimit;
            obHiGamma = gamma;
        end
        
        % Error if temperature ratio is less than 0
        if any(varargin{1}<0)
            error(message('aero:flowfanno:temperatureRatio', sprintf( '%.17f', obLoTemp( 1 ) ), sprintf( '%.15f', obLoGamma( 1 ) )));
        end
        
        % Error if temperature ratio is greater than (GAMMA+1)/2.
        if any(varargin{1}>((gamma+1)/2))
            error(message('aero:flowfanno:temperatureRatio', sprintf( '%.17f', obHiTemp( 1 ) ), sprintf( '%.15f', obHiGamma( 1 ) )));
        end
        
        % The second input (first variable input) is temperature ratio T
        T = varargin{1};
        
        % Find Mach number from temperature ratio with this closed form
        % equation.
        mach = sqrt((gamma+1)./(gamma-1).*T.^(-1)-2./(gamma-1));
        
    case 'pressratio'
        
        % Error if pressure ratio is negative
        if any(varargin{1}<0)
            error(message('aero:flowfanno:pressureRatio'));
        end
        
        % The second input (first variable input) is pressure ratio P
        P = varargin{1};
        
        % Find Mach number from pressure ratio with this closed form equation.
        mach = sqrt(-1./(gamma-1)+sqrt(1./(gamma-1).^2+(gamma+1)./(gamma-1).*P.^-2));
        if any(P==inf)
            mach(P==inf) = 0;
        end
        
    case 'densratio'
        
        % Force values very close to the limit density to be the limit
        lowerLimit = sqrt((gamma-1)./(gamma+1));
        
        if ~isscalar(gamma)
            obLoDens = lowerLimit(varargin{1}<lowerLimit);
            obLoGamma = gamma(varargin{1}<lowerLimit);
        else
            obLoDens = lowerLimit;
            obLoGamma = gamma;
        end
        
        % Error if density ratio is less than SQRT((GAMMA-1)/(GAMMA+1))
        if any(varargin{1}<lowerLimit)
            error(message('aero:flowfanno:densityRatio', sprintf( '%.17f', obLoDens( 1 ) ), sprintf( '%.15f', obLoGamma( 1 ) )));
        end
        
         % The second input (first variable input) is density ratio rho
         rho = varargin{1};
         
         % Find Mach number from density ratio with this closed form equation.
         mach = sqrt((rho).^(-2).*(2./(gamma+1))./(1-((gamma-1)./(gamma+1).*rho.^(-2))));
                                                
    case 'veloratio'
        
        % Force values very close the limit velocity to be the limit
        upperLimit = sqrt((gamma+1)./(gamma-1));
        
        if ~isscalar(gamma)
            obLoVelo = upperLimit(varargin{1}<0);
            obLoGamma = gamma(varargin{1}<0);
            obHiVelo = upperLimit(varargin{1}>upperLimit);
            obHiGamma = gamma(varargin{1}>upperLimit);
        else
            obLoVelo = upperLimit;
            obLoGamma = gamma;
            obHiVelo = upperLimit;
            obHiGamma = gamma;
        end
        
        % Error if the velocity ratio is less than 0
        if any(varargin{1}<0)
            error(message('aero:flowfanno:velocityRatio', sprintf( '%.17f', obLoVelo( 1 ) ), sprintf( '%.15f', obLoGamma( 1 ) )));
        end
        
        % Error if the velocity ratio is greater than
        % SQRT((GAMMA+1)/(GAMMA-1))
        if any(varargin{1}>upperLimit)
            error(message('aero:flowfanno:velocityRatio', sprintf( '%.17f', obHiVelo( 1 ) ), sprintf( '%.17f', obHiGamma( 1 ) )));
        end
        
        % The second input (first variable input) is velocity ratio V
        V = varargin{1};
        
        % Find Mach number from velocity ratio with this closed form equation.
        mach = sqrt(2./(gamma+1).*V.^2./(1-((gamma-1)./(gamma+1)).*V.^2));
        
        if any(V == upperLimit)
            mach(V==upperLimit) = inf;
        end
        
    case 'subsonictotalpressureratio'
        
        % The second input (first variable input) is total pressure ratio
        P0 = varargin{1};                                           % P0
        
        % Find the total pressure ratio errors with this function (below)
        totalPressureRatioSpecificErrors(gamma, P0);
        
        % Find the Mach number through this function (below)
        mach = machFromTotalPressureRatio(gamma, P0, param);
        
    case 'supersonictotalpressureratio'
        
        % The second input (first variable input) is total pressure ratio
        P0 = varargin{1};                                           % P0
        
        % Find the total pressure ratio errors with this function (below)
        totalPressureRatioSpecificErrors(gamma, P0);
        
        % Find the Mach number through this function (below)
        mach = machFromTotalPressureRatio(gamma, P0, param);    
        
    case 'subsonicfannoparameter'
        
        % The second input (first variable input) is Fanno parameter F.
        F = varargin{1};
        
        % Find the Fanno parameter errors with this function (below)
        fannoParameterSpecificErrors(gamma, F, param);
        
        % Find the Mach number through this function (below)
        mach = machFromFannoParameter(gamma, F, param);
        
    case 'supersonicfannoparameter'
        
        % The second input (first variable input) is Fanno parameter F.
        F = varargin{1};
        
        % Find the Fanno parameter errors with this function (below)
        fannoParameterSpecificErrors(gamma, F, param);
        
        % Find the Mach number through this function (below)
        mach = machFromFannoParameter(gamma, F, param);    
        
end

% Function Body (once Mach number has been given or deduced)

% Static to sonic temperature ratio, T = Temp/Temp*
T=1./((2./(gamma+1)).*(1+(gamma-1)./2.*mach.^2));

% Static to sonic pressure ratio, P = p/p*
P=1./(mach.*sqrt((2./(gamma+1)).*(1+(gamma-1)./2.*mach.^2)));

% Static to sonic density ratio, rho = rho_local/rho*
rho=sqrt(2./(gamma+1).*(mach.^-2+(gamma-1)./2));

% Ratio of velocity at a point to sonic velocity, V = V/V*
V=sqrt(1./(2./(mach.^2.*(gamma+1))+(gamma-1)./(gamma+1)));

% Ratio of total pressure at a point over the total pressure at sonic
% conditions, P0 = p_0/p_0*
P0=1./mach.*((2./(gamma+1)).*(1+(gamma-1)./2.*mach.^2))...
    .^((gamma+1)./(2.*(gamma-1)));

% Fanno parameter, F = f*L/D, where f is the friction coefficient, L is the
% length of constant area duct required to achieve supersonic flow, and D
% is the hydraulic diameter of the duct.
F=(mach.^-2-1)./gamma+(gamma+1)./(2.*gamma).*...
    log(1./((mach.^-2.*2./(gamma+1))+(gamma-1)./(gamma+1)));

% Special case: Force Fanno parameter to be infinity if Mach number is zero
% (Otherwise evaluated as NaN)

if isequal(size(F),size(mach))
    F(mach==0)=inf;
else
    if mach==0
        F(:)=inf;
    end
end

% Special cases: Use these equation for total pressure ratio if Mach number
% is infinity (Otherwise evaluated as NaN)

b = (gamma+1)./(2.*(gamma-1));

if isequal(size(gamma),size(mach))
    P0(mach==inf) = mach(mach==inf).^(2.*b(mach==inf)-1).*(2./(mach(mach==inf).^2.*...
        (gamma(mach==inf)+1))+(gamma(mach==inf)-1)./(gamma(mach==inf)+1)).^b(mach==inf);

elseif (~isscalar(gamma) && isscalar(mach))
    if mach==inf
        P0 = mach.^(2.*b-1).*(2./(mach.^2.*(gamma+1))+(gamma-1)./(gamma+1)).^b;
    end
elseif (isscalar(gamma) && ~isscalar(mach))
    P0(mach==inf) = mach(mach==inf).^(2.*b-1).*(2./(mach(mach==inf).^2.*(gamma+1))+(gamma-1)./(gamma+1)).^b;
end

% Return the values input by the user, make the size such that all outputs
% are the same.

switch param
    case 'mach'
        mach = varargin{1}.*ones(size(gamma));
    case 'tempratio'
        T = varargin{1}.*ones(size(mach));
    case 'pressratio'
        P = varargin{1}.*ones(size(mach));
    case 'densratio'
        rho = varargin{1}.*ones(size(mach));
    case 'veloratio'
        V = varargin{1}.*ones(size(mach));
    case {'subsonictotalpressureratio','supersonictotalpressureratio'}
        P0 = varargin{1}.*ones(size(mach));
    case {'subsonicfannoparameter','supersonicfannoparameter'}
        F = varargin{1}.*ones(size(mach));
end


%--------------------------------------------------------------------------
function totalPressureRatioSpecificErrors(gamma, P0)
% This function catches errors specific to area ratio.  These errors are
% in a function here because the error messages and error IDs are not
% unique between subsonic and supersonic total pressure ratios.

% Error if specific heat ratio or total pressure ratio are not scalars for
% total pressure ratio input mode
if ~(isscalar(gamma) && isscalar(P0))
    error(message('aero:flowfanno:totalPressureRatioScalar'))
end

% Error if total pressure ratio is less than 1
if P0<1
    error(message('aero:flowfanno:totalPressureRatio'))
end

% Error if total pressure ratio is Not-a-Number
if isnan(P0)
    error(message('aero:flowfanno:nanTotalPressure'))
end

% Error if specific heat ratio is not a finite value
if ~isfinite(gamma)
    error(message('aero:flowfanno:totalPressureRatioGamma'))
end

%--------------------------------------------------------------------------
function mach = invertTotalPressureRatio(gamma, P0, param)
% This function uses a residual method of solving for the Mach number when
% given the total pressure ratio because a closed form solution of the
% total pressure ratio equation for Mach number cannot be found.  This
% method involves solving the following residual equation for zeros:
% fcn(totalPressureRatio) = (known) value of totalPressureRatio -
%                                           equation for totalPressureRatio
%                                           as a fcn(GAMMA,MACH)
%
% The Mach number is being sought so it must remain anonymous until its
% value has been found by the fzero function.  The quality of the solution
% is critically dependent on the quality of the initial guess in the fzero
% function because for each total pressure ratio there is both a subsonic
% and supersonic solution.  The solution that is found corresponds to the
% user's choice of the subsonic or supersonic value from the main function
% flowfanno above.

% This equation is fcn(totalPressureRatio) as described above
equation = @(gamma,P0,mach) P0-(1./mach.*(2./(gamma+1).*(1+(gamma-1)./2.*mach^2)).^...
                                                ((gamma+1)/(2*(gamma-1))));

                                            
if strcmpi(param,'subsonictotalpressureratio')
    
    % Subsonic guess
    guessMach = [realmin 1];
    
elseif strcmpi(param,'supersonictotalpressureratio')
    
    % Supersonic guess
    if gamma < 1.1
        % mach > 37.685 calculates P0 = Inf
        guessMach = [1 realmax/4.86e306];
    elseif gamma < 1.2
        % mach > 2e+15 calculates P0 = inf
        guessMach = [1 realmax/1e293];
    elseif gamma < 1.3
        % mach > 1e+28 calculates P0 = inf
        guessMach = [1 realmax/1e280];
    elseif gamma < 1.4
        % mach > 1e+40 calculates P0 = inf
        guessMach = [1 realmax/1e268];
    elseif gamma < 1.8
        % mach > 1e+51 calculates P0 = inf
        guessMach = [1 realmax/1e257];
    elseif gamma < 1.9
        % mach > 1e+88 calculates P0 = inf
        guessMach = [1 realmax/1e220];
    elseif gamma < 2.0
        % mach > 1e+95 calculates P0 = inf
        guessMach = [1 realmax/1e213];
    else
        % mach > 1e+102 calculates P0 = inf
        guessMach = [1 realmax/1e206];
    end
    
end

% Evaluate the Mach number which corresponds to the zero crossing of the
% fcn(totalPressureRatio) equation
machEval = fzero(@(mach) equation(gamma,P0,mach), guessMach);

% Let the Mach number to be used be the same as evaluated Mach number
mach = machEval;

%--------------------------------------------------------------------------
function mach = machFromTotalPressureRatio(gamma, P0, param)
% This function finds the Mach number when a total pressure ratio is the
% second input in the function flowfanno above.  The subfunction
% invertTotalPressureRatio has the method for solving for the Mach number.

% Find the Mach number when total pressure ratio is finite
if isfinite(P0)
    mach = invertTotalPressureRatio(gamma, P0, param);
    
elseif ~isnan(P0) % If total pressure ratio is inf
    switch param
        case 'subsonictotalpressureratio' % mach = 0 for subsonic
            mach = 0;
        case 'supersonictotalpressureratio' % mach = inf for supersonic
            mach = inf;
    end
end

%--------------------------------------------------------------------------
function fannoParameterSpecificErrors(gamma, F, param)
% This function catches errors specific to Fanno parameter.  These errors
% are in a function here because many of the error messages and error IDs
% are not unique between subsonic and supersonic Fanno parameters.

% Error if specific heat ratio or Fanno parameter are not scalars for
% Fanno parameter input mode
if ~(isscalar(gamma) && isscalar(F))
    error(message('aero:flowfanno:fannoParameterScalar'))
end

% Error if Fanno parameter is less than 0
if F<0
    error(message('aero:flowfanno:fannoParameterNegative'))
end

upperLimit = (gamma+1)/(2*gamma)*log((gamma+1)/(gamma-1))-1/gamma;

% Error if the supersonic Fanno parameter is greater than
% (GAMMA+1)/(2*GAMMA)*LOG((GAMMA+1)/(GAMMA-1))-1/GAMMA in Fanno parameter
% input mode
if strcmpi(param,'supersonicfannoparameter')
    if F>upperLimit
        error(message('aero:flowfanno:fannoParameterSupersonic', sprintf( '%.17f', upperLimit ), sprintf( '%.15f', gamma )));
    end
end

% Error if Fanno parameter is Not-a-Number
if isnan(F)
    error(message('aero:flowfanno:nanFannoParameter'))
end

% Error is specific heat ratio is not a finite value
if ~isfinite(gamma)
    error(message('aero:flowfanno:fannoParameterGamma'))
end

%--------------------------------------------------------------------------
function mach = invertFannoParameter(gamma, F, param)
% This function uses a residual method of solving of solving for the Mach
% number when given Fanno parameter because a closed form solution of the
% Fanno parameter equation for the Mach number cannot be found.  This 
% method involves solving the following residual equation for zeros:
%
% fcn(fannoParameter) = (known) value of fannoParameter - 
%                                               equation for fannoParameter
%                                               as a fcn(GAMMA,MACH)
%
% The Mach number is being sought so it must remain anonymous until its
% value has been found by the fzero function.  The quality of the solution
% is critically dependent on the quality of the initial guess in the fzero
% function because for each Fanno parameter there may be a subsonic and a
% supersonic solution.  The solution that is found corresponds to the
% user's choice of the subsonic or supersonic value from the main function
% flowfanno above.

% This equation is fcn(fannoParameter) as described above
equation = @(gamma,F,mach) F-((1-mach^2)/(gamma*mach^2)+(gamma+1)/(2*gamma)*...
    log(mach^2/(2/(gamma+1)*(1+(gamma-1)/2*mach^2))));

if strcmpi(param,'subsonicfannoparameter')
    
    % Subsonic guess
    guessMach = [sqrt(realmin) 1];
    
elseif strcmpi(param,'supersonicfannoparameter')
    
    % Supersonic guess
    guessMach = [1 sqrt(realmax)];
    
end

% Evaluate the Mach number which corresponds to the zero crossing of the
% fcn(fannoParameter) equation
machEval = fzero(@(mach) equation(gamma,F,mach), guessMach);

% Let the Mach number to be used be the same as evaluated Mach number
mach = machEval;

%--------------------------------------------------------------------------
function mach = machFromFannoParameter(gamma, F, param)
% This function finds the Mach number when a Fanno parameter is the second
% input in the function flowfanno above.  The subfunction
% invertFannoParameter has the method for solving for the Mach number.

switch param
    
    % Subsonic 
    case 'subsonicfannoparameter' 
        if isfinite(F) % Find the Mach number when Fanno parameter is finite
            mach = invertFannoParameter(gamma, F, param);
        elseif ~isnan(F)
            mach = 0;
        end
        
    case 'supersonicfannoparameter'
        
        % Upper limit of supersonic Fanno parameter
        upperLimit = (gamma+1)/(2*gamma)*log((gamma+1)/(gamma-1))-1/gamma;
        
        % Mach is infinity when the Fanno parameter is the upper limit
        if isequal(F, upperLimit)
            mach = inf;
        else % Find the Mach number when Fanno parameter is under the limit
            mach = invertFannoParameter(gamma, F, param);
        end
end
