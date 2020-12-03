function [mach, T, P, rho, V, T0, P0] = flowrayleigh(gamma, varargin)
%FLOWRAYLEIGH   Calculate Rayleigh line relations
%   [MACH, T, P, RHO, V, T0, P0] = FLOWRAYLEIGH(GAMMA, VAR, MTYPE) computes
%   an array for each of the Rayleigh line flow relations.  FLOWRAYLEIGH
%   calculates these arrays for a given set of specific heat ratios, GAMMA,
%   and any one of the Rayleigh line flow variables.  The Rayleigh flow
%   variable is selected by the string, MTYPE.  The outputs are the Mach
%   number, MACH, static temperature ratio, T, the static pressure ratio,
%   P, the density ratio, RHO, the velocity ratio, V, stagnation (total)
%   temperature ratio, T0, and the stagnation (total) pressure ratio, P0.
%   All of these ratios are static conditions over the sonic conditions.
%
%   See documentation for more information about assumptions and
%   limitations.
%
%   Inputs for FLOWRAYLEIGH are:
%
%   GAMMA   :  An array of specific heat ratios.  GAMMA must be a scalar or
%              an array of N of real numbers greater than 1.  For low speed
%              temperature ratio, high speed temperature ratio, subsonic
%              total temperature, supersonic total temperature, subsonic
%              total pressure, and supersonic total pressure input modes,
%              GAMMA must be a real, finite scalar greater than 1.
%
%   VAR     :  An array of real numerical values for one of the Rayleigh
%              line flow variables.
%
%              MACH :  An array of Mach numbers.  MACH must be a scalar or
%                      an array of real numbers greater than or equal to 0.
%                      If MACH and GAMMA are both arrays, they must be the
%                      same size.  MACH is used with MTYPE variable 'mach'.
%                      Since 'mach' is the default MTYPE value, the MTYPE
%                      variable is optional for Mach number input mode.
%
%              T    :  A scalar value of temperature ratio. The temperature
%                      ratio is defined as the local static temperature
%                      over the reference static temperature for sonic
%                      flow.  T must be a real scalar greater than or equal
%                      to 0 (at MACH = 0 for low speeds or as MACH -> inf
%                      for high speeds) and less than or equal to
%                      1/4*(GAMMA+1/GAMMA)+1/2 (at MACH = 1/SQRT(GAMMA)). 
%                      T is used with MTYPE variables 'templo' and
%                      'temphi'.
%
%              P    :  An array of pressure ratios.  The pressure ratio is
%                      defined as the local static pressure over the
%                      reference static pressure for sonic flow.  P must be
%                      a scalar or an array of real numbers less than or
%                      equal to GAMMA+1 (at MACH = 0).  If P and GAMMA are
%                      both arrays, they must be the same size.  P is used
%                      with MTYPE variable 'pres'.
%
%              RHO  :  An array of density ratios.  The density ratio is
%                      defined as the local density over the reference
%                      density for sonic flow. RHO must be a scalar or an
%                      array of real numbers greater than or equal to
%                      GAMMA/(GAMMA+1) (as MACH -> inf).  If RHO and GAMMA
%                      are both arrays, they must be the same size.  RHO is
%                      used with MTYPE variable 'dens'.
%
%              V    :  An array of velocity ratios.  The velocity ratio is
%                      defined as the local velocity over the reference
%                      velocity for sonic flow.  V must be a scalar or an
%                      array of real numbers greater than or equal to 0 and
%                      less than or equal to (GAMMA+1)/GAMMA (as 
%                      MACH -> inf).  If V and GAMMA are both arrays, they
%                      must be the same size.  V is used with MTYPE
%                      variable 'velo'.
%
%              T0   :  A scalar value of total temperature ratio.  The
%                      total temperature ratio is defined as the local
%                      stagnation temperature over the reference stagnation
%                      temperature for sonic flow.  T0 must be a real
%                      scalar greater than or equal to 0 (at MACH = 0) and
%                      less than or equal to 1 (at MACH = 1) in subsonic
%                      mode.  In supersonic mode, T0 must be a real scalar
%                      greater than or equal to
%                      (GAMMA+1)^2*(GAMMA-1)/2/(GAMMA^2*(1+(GAMMA-1)/2)))
%                      (as MACH -> inf) and less than or equal to 1 (at
%                      MACH = 1).  T0 is used with the MTYPE variables
%                      'totaltsub' and 'totaltsup'.
%
%              P0   :  A scalar value of total pressure ratio.  The total
%                      pressure ratio is defined as the local stagnation
%                      pressure over the reference stagnation pressure for
%                      sonic flow.  In subsonic mode, P0 must be a real
%                      scalar greater than or equal to 1 (at MACH = 1) and
%                      less than or equal to
%                      (1+GAMMA)*(1+(GAMMA-1)/2)^(-GAMMA/(GAMMA-1)) (at
%                      MACH = 0).  In supersonic mode, P0 must be a real
%                      scalar greater than or equal to 1.  P0 is used with
%                      the MTYPE variables 'totalpsub' and 'totalpsup'.
%
%   MTYPE   :  A character vector for selecting the isentropic flow variable
%              represented by VAR.
%
%              'mach'      :  Default value.  Indicates that the function
%                             is in Mach number input mode.  Value is of
%                             type string.
%
%              'templo'    :  Indicates that the function is in low speed
%                             static temperature ratio input mode.  The low
%                             speed temperature ratio is defined as the
%                             local static temperature over the reference
%                             sonic temperature when the Mach number of the
%                             upstream flow is less than the critical Mach
%                             number of 1/SQRT(GAMMA).  Value is of type
%                             string.
%
%               'temphi'   :  Indicates that the function is in high speed
%                             static temperature ratio input mode.  The
%                             high speed temperature ratio is defined as
%                             the local static temperature over the
%                             reference sonic temperature when the Mach
%                             number of the upstream flow is greater than
%                             the critical Mach number of 1/SQRT(GAMMA).
%                             Value is of type string.
%
%              'pres'      :  Indicates that the function is in pressure
%                             ratio input mode.  Value is of type string.
%
%              'dens'      :  Indicates that the function is in density
%                             ratio input mode.  Value is of type string.
%
%              'velo'      :  Indicates that the function is in velocity
%                             ratio input mode.  Value is of type string.
%
%              'totaltsub' :  Indicates that the function is in subsonic
%                             total temperature ratio input mode.  Value is
%                             of type string.
%
%              'totaltsup' :  Indicates that the function is in supersonic
%                             total temperature ratio input mode.  Value is
%                             of type string.
%
%              'totalpsub' :  Indicates that the function is in subsonic
%                             total pressure ratio input mode.  Value is of
%                             type string.
%
%              'totalpsup' :  Indicates that the function is in supersonic
%                             total pressure ratio input mode.  Value is of
%                             type string.
%
%   Outputs calculated for Rayleigh line flow are:
%   (All outputs are the same size as the array input or array inputs.  If
%   there are no array inputs, all outputs are scalars.)
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
%   T0      :  An array of total temperature ratios.  The temperature ratio
%              is defined as the local static temperature over the
%              reference static temperature for sonic flow.
%
%   P0      :  An array of total pressure ratios.  The total pressure ratio
%              is defined as the local stagnation pressure over the
%              reference stagnation pressure for sonic flow.
%
%   Examples:
%
%   Calculate the Rayleigh line flow relations for air (gamma = 1.4) for
%   supersonic total pressure ratio 1.2.  The following returns scalar
%   values for MACH, T, P, RHO, V, T0, and P0.
%
%      [MACH, T, P, RHO, V, T0, P0] = flowrayleigh(1.4, 1.2, 'totalpsup')
%
%   Calculate the Rayleigh line flow relations for gases with specific heat
%   ratios given in the following 1 x 4 row array for the Mach number 0.5.
%   The following yields a 1 x 4 row array for MACH, T, P, RHO, V, T0, and
%   P0.
%
%      gamma = [1.3, 1.33, 1.4, 1.67]
%      [MACH, T, P, RHO, V, T0, P0] = flowrayleigh(gamma, 0.5)
%
%   Calculate the Rayleigh line flow relations for a specific heat ratio of
%   1.4 and high speed temperature ratio 0.70.  The following returns
%   scalar values for MACH, T, P, RHO, V, T0, and P0.
%
%      [MACH, T, P, RHO, V, T0, P0] = flowrayleigh(1.4, 0.70, 'temphi')
%
%   Calculate the Rayleigh line flow relations for gases with specific heat
%   ratio and static pressure ratio combinations as shown.  The following
%   returns a 1 x 2 array for MACH, T, P, RHO, V, T0, and P0 each, where
%   the elements of each array corresponds to the inputs element-wise.
%
%       gamma = [1.3, 1.4]
%       P = [0.13, 1.7778]
%       [MACH, T, P, RHO, V, T0, P0] = flowrayleigh(gamma, P, 'pres')
%       
%   See also FLOWISENTROPIC, FLOWNORMALSHOCK, FLOWPRANDTLMEYER, FLOWFANNO.


%   Copyright 2009-2018 The MathWorks, Inc.

%   References:  James, J. E. A., Gas Dynamics, Second Edition, Allyn and 
%   Bacon, Inc, Boston, 1984.
%
%   NACA Technical Report 1135, 1953, National Advisory Committee on
%   Aeronautics, Ames Research Staff, Moffett Field, Calif.

% Checks

% General checks

% Less than two arguments check or more than three inputs check
narginchk(2,3);

% Error if third input exists and is not an acceptable string
if nargin>2
    if ~ischar(varargin{2}) && ~isstring(varargin{2})
        error(message('aero:flowrayleigh:paramSelectString'))
    end
end

% Check gamma

% gamma non-numeric check
if ~isnumeric(gamma)
    error(message('aero:flowrayleigh:notNumericGamma'));
end

% gamma > 1 check (specific heat ratio must be greater than one)
if any(gamma<=1)
    error(message('aero:flowrayleigh:gammaOneOrLess'));
end

% gamma real number check
if ~isreal(gamma)
    error(message('aero:flowrayleigh:imaginaryGamma'));
end

% Error if gamma input is a matrix or if second input variable is a matrix
if (~isvector(gamma) || ~isvector(varargin{1}))
    error(message('aero:flowrayleigh:multiDimensional'))
end

% Second input general checks (varargin{1})

% Second input variable non-numeric check
if ~isnumeric(varargin{1})
    error(message('aero:flowrayleigh:notNumeric'));
end

% Error if inputs are not the same size and if neither are scalars
if ~((isscalar(gamma) || isscalar(varargin{1})) || ...
        (isequal(size(gamma),size(varargin{1}))))
    error(message('aero:flowrayleigh:size'))
end

% Second input variable real number check
if ~isreal(varargin{1})
    error(message('aero:flowrayleigh:imaginary'));
end

% Initialize the string variable as the default (mach)
param = 'mach';

if nargin > 2 % If there are selector strings for the second input
    
    % Making string inputs more flexible
    
    % Typing the string 'mach' to input Mach number is optional
    if strcmpi(varargin{2},'mach')
        param = 'mach';
        
    % All the user really has to type in is "templo' to input
    % lowSpeedTemperatureRatio.
    elseif strncmpi(varargin{2},'templospeed',6)
        param = 'tempratiolospeed';
        
    % All the user really has to type in is "temphi' to input
    % highSpeedTemperatureRatio.
    elseif strncmpi(varargin{2},'temphispeed',6)
        param = 'tempratiohispeed';
        
    % All the user really has to type in is 'pres' to input pressureRatio.
    elseif strncmpi(varargin{2},'pressratio',4)
        param = 'pressratio';
        
    % All the user really has to type in is 'dens' to input densityRatio.
    elseif strncmpi(varargin{2},'densratio',4)
        param = 'densratio';
        
    % All the user really has to type in is 'velo' to input velocityRatio.
    elseif strncmpi(varargin{2},'veloratio',4)
        param = 'veloratio';
        
    % All the user really has to type in is 'totaltsub' to input
    % subsonicTotalTemperatureRatio.
    elseif strncmpi(varargin{2},'totaltsub',9)
        param = 'subsonictotaltemperatureratio';
        
    % All the user really has to type in is 'totaltsup' to input
    % supersonicTotalTemperatureRatio.
    elseif strncmpi(varargin{2},'totaltsup',9)
        param = 'supersonictotaltemperatureratio';
        
    % All the user really has to type in is 'totalpsub' to input
    % subsonicTotalPressureRatio.
    elseif strncmpi(varargin{2},'totalpsub',9)
        param = 'subsonictotalpressureratio';
        
    % All the user really has to type in is 'totalpsup' to input
    % supersonicTotalPressureRatio.
    elseif strncmpi(varargin{2},'totalpsup',9)
        param = 'supersonictotalpressureratio';
    else
        error(message('aero:flowrayleigh:paramSelectWrongInput'))
    end
end

% Specific cases

switch lower(param)
    
    case 'mach'
        
        % MachNo non-negative check (Mach number must be greater than or equal to 0
        if any(varargin{1}<0)
            error(message('aero:flowrayleigh:negative'));
        end

        % The second input (first variable input) is Mach number.
        mach = varargin{1};
        
    case 'tempratiolospeed'

        % The second input (first variable input) is low speed temperature
        % ratio T
        T = varargin{1};
        
        % Find the temperature ratio errors with this function (below)
        tempRatioSpecificErrors(gamma, T);
        
        % Find the Mach number through this function (below)
        mach = machFromTempRatio(gamma, T, param);
        
    case 'tempratiohispeed'
        
        % The second input (first variable input) is high speed temperature
        % ratio T
        T = varargin{1};
        
        % Find the temperature ratio errors with this function (below)
        tempRatioSpecificErrors(gamma, T);
        
        % Find the Mach number through this function (below)
        if isequal(T,0)
            mach = inf;
        else
            mach = machFromTempRatio(gamma, T, param);
        end
        
    case 'pressratio'
        
        upperLimit = gamma+1;
        
        if ~isscalar(gamma)
            obLoPres = upperLimit(varargin{1}<0);
            obLoGamma = gamma(varargin{1}<0);
            obHiPres = upperLimit(varargin{1}>upperLimit);
            obHiGamma = gamma(varargin{1}>upperLimit);
        else
            obLoPres = upperLimit;
            obLoGamma = gamma;
            obHiPres = upperLimit;
            obHiGamma = gamma;
        end
        
        
        % Error if the pressure ratio is negative
        if any(varargin{1}<0)
            error(message('aero:flowrayleigh:pressureRatio', sprintf( '%.17f', obLoPres( 1 ) ), sprintf( '%.15f', obLoGamma( 1 ) )))
        end
        
        % Error if the pressure ratio is greater than GAMMA+1
        if any(varargin{1}>(gamma+1))
            error(message('aero:flowrayleigh:pressureRatio', sprintf( '%.17f', obHiPres( 1 ) ), sprintf( '%.15f', obHiGamma( 1 ) )))
        end
        
        % The second input (first variable input) is pressure ratio P
        P = varargin{1};
        
        % Find Mach number from pressure ratio with this closed form equation.

            mach = ((1+gamma)./(gamma.*P)-1./gamma).^(1/2);
            mach(P==gamma+1) = 0;
            
    case 'densratio'
        
        lowerLimit = gamma./(gamma+1);
        
        if ~isscalar(gamma)
            obLoDens = lowerLimit(varargin{1}<lowerLimit);
            obLoGamma = gamma(varargin{1}<lowerLimit);
        else
            obLoDens = lowerLimit;
            obLoGamma = gamma;
        end
        
        % Error if the density ratio is less than GAMMA/(GAMMA+1)
        if any(varargin{1}<lowerLimit)
            error(message('aero:flowrayleigh:densityRatio', sprintf( '%.17f', obLoDens( 1 ) ), sprintf( '%.15f', obLoGamma( 1 ) )));
        end
        
        % The second input (first variable input) is density ratio rho
        rho = varargin{1};
        
        % Find Mach number from density ratio with this closed form equation.
        mach = (1./(rho.*(gamma+1)-gamma)).^(1/2);
        
    case 'veloratio'
        
        upperLimit = (gamma+1)./gamma;
        
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
        
        % Error if the velocity ratio is negative
        if any(varargin{1}<0)
            error(message('aero:flowrayleigh:velocityRatio', sprintf( '%.17f', obLoVelo( 1 ) ), sprintf( '%.15f', obLoGamma( 1 ) )));
        end
        
        % Error if the velocity ratio is greater than (GAMMA+1)/GAMMA
        if any(varargin{1}>upperLimit)
            error(message('aero:flowrayleigh:velocityRatio', sprintf( '%.17f', obHiVelo( 1 ) ), sprintf( '%.15f', obHiGamma( 1 ) )));
        end
        
        % The second input (first variable input) is velocity ratio V
        V = varargin{1};
        
        % Find Mach number from velocity ratio with this closed form equation.
        mach = (-V./(gamma.*V-(1+gamma))).^(1/2);
        
     case 'subsonictotaltemperatureratio'
        
        % The second input (first variable input) is total temperature
        % ratio T0
        T0 = varargin{1};                                          
        
        % Find the total temperature ratio errors with this function (below)
        totalTemperatureRatioSpecificErrors(gamma, T0, param);
        
        % Find the Mach number through this function (below)
        
        mach = machFromTotalTemperatureRatio(gamma, T0, param);
        
    case 'supersonictotaltemperatureratio'
                
        % The second input (first variable input) is total temperature
        % ratio T0
        T0 = varargin{1};                                 
        
        % Find the total temperature ratio errors with this function (below)
        totalTemperatureRatioSpecificErrors(gamma, T0, param);
        
        lowerLimit = (gamma+1)^2*(gamma-1)/2/(gamma^2*(1+(gamma-1)/2));
                
        % Find the Mach number through this function (below)
        if T0 == lowerLimit
            mach = inf;
        else
            mach = machFromTotalTemperatureRatio(gamma, T0, param);
        end
        
     case 'subsonictotalpressureratio'
        
        % The second input (first variable input) is total pressure ratio
        P0 = varargin{1};                                           % P0
        
        % Find the total pressure ratio errors with this function (below)
        totalPressureRatioSpecificErrors(gamma, P0, param);
        
        % Find the Mach number through this function (below)
        mach = machFromTotalPressureRatio(gamma, P0, param);
        
    case 'supersonictotalpressureratio'
        
        % The second input (first variable input) is total pressure ratio
        P0 = varargin{1};                                           % P0
        
        % Find the total pressure ratio errors with this function (below)
        totalPressureRatioSpecificErrors(gamma, P0, param);
        
        % Find the Mach number through this function (below)
        if P0 == inf
            mach = inf;
        else
            mach = machFromTotalPressureRatio(gamma, P0, param);
        end
        
end

% Function Body

% Static to sonic temperature ratio, T = Temp/Temp*
T=(1+gamma).^2.*mach.^2./(1+gamma.*mach.^2).^2;

% Static to sonic pressure ratio, P = p/p*
P=(1+gamma)./(1+gamma.*mach.^2);

% Static to sonic density ratio, rho = rho_local/rho*
rho=(1+gamma.*mach.^2)./((1+gamma).*mach.^2);

% Ratio of velocity at a point to sonic velocity, V = V/V*
V=(1+gamma).*mach.^2./(1+gamma.*mach.^2);

% Ratio of total pressure at a point over the total pressure at sonic
% conditions, T0 = Temp_0/Temp_0*
T0=(1+gamma).^2.*mach.^2./(1+gamma.*mach.^2).^2.*...
    (1+(gamma-1)./2.*mach.^2)./(1+(gamma-1)./2);

% Ratio of total pressure at a point over the total pressure at sonic
% conditions, P0 = p_0/p_0*
P0=(1+gamma)./(1+gamma.*mach.^2).*((1+(gamma-1)./2)./...
    (1+(gamma-1)./2.*mach.^2)).^(-gamma./(gamma-1));

% Special cases: Use these rearranged equations if Mach number is infinity
% (Otherwise evaluated as NaN)
if isequal(size(gamma),size(mach))
    
    T(mach==inf)=(1+2.*gamma(mach==inf)+gamma(mach==inf).^2)./...
        (mach(mach==inf).^2.*(mach(mach==inf).^-4+ ...
        2.*gamma(mach==inf).*mach(mach==inf).^-2+gamma(mach==inf).^2));
    
    rho(mach==inf)=((gamma(mach==inf)+1).*mach(mach==inf).^2).^-1 + ...
        (gamma(mach==inf)./(gamma(mach==inf)+1));    
    
    V(mach==inf)=(1+gamma(mach==inf))./(mach(mach==inf).^-2+gamma(mach==inf));    
    
    T0(mach==inf)=(1+gamma(mach==inf)).^2.*(mach(mach==inf).^...
        -2+(gamma(mach==inf)-1)./2)./...
        ((mach(mach==inf).^-4+2.*gamma(mach==inf).*mach(mach==inf).^...
        -2+gamma(mach==inf).^2).*(1+(gamma(mach==inf)-1)./2));
    
    P0(mach==inf)=mach(mach==inf).^(2.*gamma(mach==inf)./...
        (gamma(mach==inf)-1)-2).*(1+gamma(mach==inf))./...
        (mach(mach==inf).^-2+gamma(mach==inf)).*...
        ((1+(gamma(mach==inf)-1)./2)./(mach(mach==inf).^...
        -2+(gamma(mach==inf)-1)./2)).^...
        (-gamma(mach==inf)./(gamma(mach==inf)-1));
    
    
elseif (isscalar(gamma) && ~isscalar(mach))
        
        T(mach==inf)=(1+2.*gamma+gamma.^2)./(mach(mach==inf).^2.*(mach(mach==inf).^-4+ ...
        2.*gamma.*mach(mach==inf).^-2+gamma.^2));
        
        rho(mach==inf)=((gamma+1).*mach(mach==inf).^2).^-1 + (gamma./(gamma+1));
        
        V(mach==inf)=(1+gamma)./(mach(mach==inf).^-2+gamma);
        
        T0(mach==inf)=(1+gamma).^2.*(mach(mach==inf).^-2+(gamma-1)./2)./...
        ((mach(mach==inf).^-4+2.*gamma.*mach(mach==inf).^-2+gamma.^2).*(1+(gamma-1)./2));
        
        P0(mach==inf)=mach(mach==inf).^(2.*gamma./(gamma-1)-2).*(1+gamma)./...
        (mach(mach==inf).^-2+gamma).*((1+(gamma-1)./2)./(mach(mach==inf).^-2+(gamma-1)./2)).^...
        (-gamma./(gamma-1));
        
elseif (~isscalar(gamma) && isscalar(mach))
    if mach==inf
        
    T=(1+2.*gamma+gamma.^2)./(mach.^2.*(mach.^-4+ ...
        2.*gamma.*mach.^-2+gamma.^2));
    
    rho=((gamma+1).*mach.^2).^-1 + (gamma./(gamma+1));
    
    V=(1+gamma)./(mach.^-2+gamma);
    
    T0=(1+gamma).^2.*(mach.^-2+(gamma-1)./2)./...
        ((mach.^-4+2.*gamma.*mach.^-2+gamma.^2).*(1+(gamma-1)./2));
    
    P0=mach.^(2.*gamma./(gamma-1)-2).*(1+gamma)./...
        (mach.^-2+gamma).*((1+(gamma-1)./2)./(mach.^-2+(gamma-1)./2)).^...
        (-gamma./(gamma-1));
    
    end
end

% Return the values input by the user, make the size such that all outputs
% are the same.

switch param
    case 'mach'
        mach = varargin{1}.*ones(size(gamma));
    case {'tempratiolospeed' 'tempratiohispeed'}
        T = varargin{1}.*ones(size(mach));
    case 'pressratio'
        P = varargin{1}.*ones(size(mach));
    case 'densratio'
        rho = varargin{1}.*ones(size(mach));
    case 'veloratio'
        V = varargin{1}.*ones(size(mach));
    case {'subsonictotaltemperatureratio','supersonictotaltemperatureratio'}
        T0 = varargin{1}.*ones(size(mach));
    case {'subsonictotalpressureratio','supersonictotalpressureratio'}
        P0 = varargin{1}.*ones(size(mach));
end

%--------------------------------------------------------------------------
function tempRatioSpecificErrors(gamma, T)
% This function catches errors specific to temperature ratio.  These errors
% are in a function here because the error messages and error IDs are not
% unique between low-speed and high-speed temperature ratios.

% Error if the specific heat ratio or temperature ratio are not scalars for
% total pressure ratio input mode
if ~(isscalar(gamma) && isscalar(T))
    error(message('aero:flowrayleigh:temperatureRatioScalar'))
end

upperLimit = 1/4*(gamma+1/gamma)+1/2;

obLoTemp = upperLimit(T<0);
obLoGamma = gamma(T<0);
obHiTemp = upperLimit(T>upperLimit);
obHiGamma = gamma(T>upperLimit);

% Error if the temperature ratio is negative
if T<0
    error(message('aero:flowrayleigh:temperatureRatioLow', sprintf( '%.17f', obLoTemp ), sprintf( '%.15f', obLoGamma )))
end

% Error if the temperature ratio is greater than: 1/4*(GAMMA+1/GAMMA)+1/2
if T>1/4*(gamma+1/gamma)+1/2
    error(message('aero:flowrayleigh:temperatureRatioHigh', sprintf( '%.17f', obHiTemp ), sprintf( '%.15f', obHiGamma )))
end

% Error if temperature ratio is Not-a-Number
if isnan(T)
    error(message('aero:flowrayleigh:nanTemperatureRatio'))
end

% Error if specific heat ratio is not a finite value in temperature ratio
% input mode
if ~isfinite(gamma)
    error(message('aero:flowrayleigh:temperatureRatioGamma'))
end

%--------------------------------------------------------------------------
function mach = machFromTempRatio(gamma, T, param)
% This function finds the Mach number when a temperature ratio is the
% second input in the function flowrayleigh above.  This function uses a
% residual method of solving for the Mach number because a closed form
% solution of the temperature ratio equation for the Mach number cannot be
% found.  This method involves solving the following residual equation for
% zeros:
% fcn(tempRatio) = (known) value of tempRatio - equation for tempRatio as
%                                               a fcn(GAMMA, MACH)
%
% The Mach number is being sought so it must remain anonymous until its
% value has been found by the fzero function.  The quality of the solution
% is critically dependent on the quality of the initial guess in the fzero
% function because for each total temperature ratio there is both a
% subsonic and supersonic solution.  The solution that is found corresponds
% to the user's choice of the low-speed or high-speed value from the main
% function flowrayleigh above.

maxTempMach = 1/sqrt(gamma);

% This equation is fcn(tempRatio) as described above
equation = @(gamma,T,mach) T-((1+gamma).^2.*mach.^2./(1+gamma.*mach.^2).^2);

if strcmpi(param,'tempratiolospeed')
    
    % Low speed Mach number guess range
    guessMach = [0 maxTempMach];
    
elseif strcmpi(param,'tempratiohispeed')
    
    % High speed Mach number guess range
    if gamma < 6 
        guessMach = [maxTempMach realmax/1e155];
    elseif gamma < 1e100
        guessMach = [maxTempMach realmax/1e255];
    end  
end

% Evaluate the Mach number which corresponds to the zero crossing of the
% fcn(tempRatio) equation
machEval = fzero(@(mach) equation(gamma,T,mach), guessMach);

% Let the Mach number to be used be the same as evaluated Mach number
mach = machEval;

%--------------------------------------------------------------------------
function totalTemperatureRatioSpecificErrors(gamma, T0, param)
% This function catches errors specific to total temperature ratio.  These
% errors are in a function here because many of the error messages and
% error IDs are not unique between subsonic and supersonic total
% temperature ratios.

% Error if specific heat ratio or total temperature ratio are not scalars
% for total temperature ratio input mode
if ~(isscalar(gamma) && isscalar(T0))
    error(message('aero:flowrayleigh:totalTemperatureRatioScalar'))
end

% Error if total temperature ratio is greater than 1
if T0>1
    error(message('aero:flowrayleigh:totalTemperatureRatioLow'))
end

% Error if subsonic total temperature ratio is negative
if T0<0
    error(message('aero:flowrayleigh:totalTemperatureRatioHigh'))
end

lowerLimit = (gamma+1)^2*(gamma-1)/2/(gamma^2*(1+(gamma-1)/2));

obLoTotalTemp = lowerLimit(T0<lowerLimit);
obLoGamma = gamma(T0<lowerLimit);

% Error if supersonic total temperature ratio is less than
% (GAMMA+1)^2*(GAMMA-1)/2/(GAMMA^2*(1+(GAMMA-1)/2))
if strcmpi(param, 'supersonictotaltemperatureratio')
    if T0<lowerLimit
        error(message('aero:flowrayleigh:totalTemperatureRatio', sprintf( '%.17f', obLoTotalTemp ), sprintf( '%.15f', obLoGamma )))
    end
end

% Error if total temperature ratio is Not-a-Number
if isnan(T0)
    error(message('aero:flowrayleigh:nanTotalTemperatureRatio'))
end

% Error if specific heat ratio is not a finite value in total temperature
% input mode
if ~isfinite(gamma)
    error(message('aero:flowrayleigh:totalTemperatureRatioGamma'))
end

%--------------------------------------------------------------------------
function mach = machFromTotalTemperatureRatio(gamma, T0, param)
% This function finds the Mach number when a total temperature ratio is the
% second input in the function flowrayleigh above.  This function uses a
% residual method of solving for the Mach number because a closed form
% solution of the total temperature ratio equation for the Mach number
% cannot  be found.  This method involves solving the following residual
% equation for zeros:
% fcn(totalTemperatureRatio) = (known) value of totalTemperatureRatio -
%                                       equation for totalTemperatureRatio as
%                                       a fcn(GAMMA, MACH)
%
% The Mach number is being sought so it must remain anonymous until its
% value has been found by the fzero function.  The quality of the solution
% is critically dependent on the quality of the initial guess in the fzero
% function because for each total temperature ratio there is both a
% subsonic and supersonic solution.  The solution that is found corresponds
% to the user's choice of the subsonic or supersonic value from the main
% function flowrayleigh above.

% This equation is fcn(totalTemperatureRatio) as described above
equation = @(gamma,T0,mach) T0-((1+gamma).^2.*mach.^2./(1+gamma.*mach.^2).^2.*...
    (1+(gamma-1)./2.*mach.^2)./(1+(gamma-1)./2));

if strcmpi(param,'subsonictotaltemperatureratio')
    
    % Subsonic guess
    guessMach = [0 1];
    
elseif strcmpi(param,'supersonictotaltemperatureratio')

    % Supersonic guess range
    if gamma < 6 
        guessMach = [1 realmax/1e155];
    elseif gamma < 1e100
        guessMach = [1 realmax/1e255];
    end
    
end

% Evaluate the Mach number which corresponds to the zero crossing of
% fcn(totalTemperatureRatio) equation
machEval = fzero(@(mach) equation(gamma,T0,mach), guessMach);

% Let the Mach number to be used be the same as evaluated Mach number
mach = machEval;

%--------------------------------------------------------------------------
function totalPressureRatioSpecificErrors(gamma, P0, param)
% This function catches errors specific to total pressure ratio.  These
% errors are in a function here because many of the error messages and
% error IDs are not unique between subsonic an supersonic total pressure
% ratios. 

% Error if the specific heat ratio or total pressure ratio are not scalars
% for total pressure ratio input mode
if ~(isscalar(gamma) && isscalar(P0))
    error(message('aero:flowrayleigh:totalPressureRatioScalar'));
end

% Error if total pressure ratio is less than 1
if P0<1
    error(message('aero:flowrayleigh:totalPressureRatioLow'))
end

upperLimit = (1+gamma)*(1+(gamma-1)/2)^(-gamma/(gamma-1));

obHiTotalPres = upperLimit(P0>upperLimit);
obHiGamma = gamma(P0>upperLimit);

% Error if subsonic total pressure ratio is greater than
% (1+GAMMA)*(1+(GAMMA-1)/2)^(-GAMMA/(GAMMA-1))
if strcmpi(param,'subsonictotalpressureratio')
    if P0>upperLimit
        error(message('aero:flowrayleigh:totalPressureRatioHigh', sprintf( '%.17f', obHiTotalPres ), sprintf( '%.15f', obHiGamma )))
    end
end

% Error if total pressure ratio is Not-a-Number
if isnan(P0)
    error(message('aero:flowrayleigh:nanTotalPressureRatio'))
end

% Error if specific heat ratio is not a finite value in total pressure
% ratio input mode
if ~isfinite(gamma)
    error(message('aero:flowrayleigh:totalPressureRatioGamma'))
end

%--------------------------------------------------------------------------
function mach = machFromTotalPressureRatio(gamma, P0, param)
% This function finds the Mach number when a total pressure ratio is the
% second input in the function flowrayleigh above.  This function uses a
% residual method of solving for the Mach number because a closed form
% solution of the total pressure ratio equation for the Mach number cannot
% be found.  This method involves solving the following residual equation
% for zeros:
% fcn(totalPressureRatio) = (known) value of totalPressureRatio -
%                                       equation for totalPressureRatio as
%                                       a fcn(GAMMA, MACH)
%
% The Mach number is being sought so it must remain anonymous until its
% value has been found by the fzero function.  The quality of the solution
% is critically dependent on the quality of the initial guess in the fzero
% function because for each total pressure ratio there is both a subsonic
% and supersonic solution.  The solution that is found corresponds to the
% user's choice of the subsonic or supersonic value from the main function
% flowrayleigh above.

% This equation is fcn(totalPressureRatio) as described above
equation = @(gamma,P0,mach) P0-((1+gamma)./(1+gamma.*mach.^2).*((1+(gamma-1)./2)./...
    (1+(gamma-1)./2.*mach.^2)).^(-gamma./(gamma-1)));

if strcmpi(param,'subsonictotalpressureratio')
    
    % Subsonic guess
    guessMach = [0 1];
    
elseif strcmpi(param,'supersonictotalpressureratio')

    % Supersonic guess range
    if gamma < 6 
        guessMach = [1 (sqrt(realmax))^((gamma-1)/gamma)];
    elseif gamma < 1e100
        guessMach = [1 realmax/1e255];
    end
    
end

% Evaluate the Mach number which corresponds to the zero crossing of the
% fcn(totalPressureRatio) as described above
machEval = fzero(@(mach) equation(gamma,P0,mach), guessMach);

% Let the Mach number to be used be the same as evaluated Mach number
mach = machEval;
