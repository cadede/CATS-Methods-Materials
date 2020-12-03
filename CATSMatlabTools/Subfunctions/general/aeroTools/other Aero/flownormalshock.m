function [mach, T, P, rho, M, P0, P1] = flownormalshock(gamma, varargin)
%FLOWNORMALSHOCK    Calculate normal shock relations
%   [MACH, T, P, RHO, M, P0, P1] = FLOWNORMALSHOCK(GAMMA, VAR, MTYPE)
%   produces an array for each of the normal shock relations.
%   FLOWNORMALSHOCK calculates these arrays for a given set of specific
%   heat ratios, GAMMA, and any one of the normal shock relations.  The
%   normal shock relations is selected by the string, MTYPE.  The outputs
%   are the Mach number, MACH, temperature ratio, T, the (static) pressure
%   ratio, P, the density ratio, RHO, downstream Mach number, M, and the
%   total (stagnation) pressure ratio, P0.  All of these ratios are
%   downstream value over upstream value.  P1 is the Rayleigh-Pitot ratio,
%   which is the ratio of upstream static pressure over the downstream
%   stagnation pressure.  Note that upstream is said to be "before" or
%   "ahead" of the shock.  Downstream is said to be "after" or "behind" the
%   shock.
%
%   See documentation for more information about assumptions and
%   limitations.
%
%
%   Inputs for FLOWNORMALSHOCK are:
%
%   GAMMA :  An array of N specific heat ratios.  GAMMA must be a scalar
%            or an array of N real numbers greater than 1.  For temperature
%            ratio, total pressure ratio, and Rayleigh-Pitot ratio input
%            modes, GAMMA must be a real, finite scalar greater than 1.
%
%   VAR   :  An array of real numerical values for one of the normal shock
%            relations.
%
%            MACH :  An array of upstream Mach numbers.  MACH must be a
%                    scalar or an array of N real numbers greater than or
%                    equal to 1.  If MACH and GAMMA are both arrays they
%                    must be the same size.  MACH is used with the MTYPE
%                    variable 'mach'.  Since 'mach' is the default MTYPE
%                    value, the MTYPE variable is optional for Mach number
%                    input mode.
%
%            T    :  A scalar value of temperature ratio.  The temperature
%                    ratio is defined as the static temperature downstream
%                    of the shock over the static temperature upstream of
%                    the shock.  T must be a real scalar greater than or
%                    equal to 1.  T is used with MTYPE variable 'temp'.
%
%            P    :  An array of pressure ratios.  The pressure ratio is
%                    defined as the static pressure downstream of the shock
%                    over the static pressure upstream of the shock.  P
%                    must be a scalar or an array of real numbers greater
%                    than or equal to 1.  If P and GAMMA are both arrays,
%                    they must be the same size.  P is used with MTYPE
%                    variable 'pres'.
%
%            RHO  :  An array of density ratios.  The density ratio is
%                    defined as the density of the fluid downstream of the
%                    shock over the density upstream of the shock.  RHO
%                    must be a scalar or an array of real numbers greater
%                    than or equal to 1 (at MACH = 1) and less than or
%                    equal to (GAMMA+1)./(GAMMA-1) (as MACH -> inf).  If
%                    RHO and GAMMA are both arrays, they must be the same
%                    size.  RHO is used with MTYPE variable 'dens'.
%
%            M    :  An array of downstream Mach numbers.  M must be a
%                    scalar or an array of real numbers greater than or
%                    equal to 0 (as MACH -> inf) and less than or equal to
%                    SQRT((GAMMA-1)./(2.*GAMMA)) (at MACH = 1).  If M and
%                    GAMMA are both arrays, they must be the same size.  M
%                    is used with MTYPE variable 'down'.  
%
%            P0   :  A scalar value of total pressure ratio.  The total
%                    pressure ratio is defined as the total pressure
%                    downstream of the shock over the total pressure
%                    upstream of the shock.  P0 must be a real scalar
%                    greater than or equal to 0 (as MACH -> inf) and less
%                    than or equal to 1 (at MACH = 1).  If P0 and GAMMA are
%                    both arrays, they must be the same size.  P0 is used
%                    with MTYPE variable 'totalp'.
%
%            P1   :  A scalar value of Rayleigh-Pitot ratio.  The
%                    Rayleigh-Pitot ratio is defined as the static pressure
%                    upstream of the shock over the total pressure
%                    downstream of the shock.  P1 must be a real scalar
%                    greater than or equal to 0 (as MACH -> inf) and less
%                    than or equal to ((GAMMA+1)./2).^(-GAMMA/(GAMMA-1))
%                    (at MACH = 1).  If P0 and GAMMA are both arrays, they
%                    must be the same size.  P1 is used with MTYPE variable
%                    'pito'.
%
%   MTYPE  :  A character vector for selecting the isentropic flow variable
%             represented by VAR.
%
%             'mach'  :  Default value.  Indicates that the function is in
%                        Mach number input mode.
%   
%             'temp'  :  Indicates that the function is in temperature
%                        ratio input mode.
%
%             'pres'  :  Indicates that the function is in pressure ratio
%                        input mode.
%
%             'dens'  :  Indicates that the function is in density ratio
%                        input mode.
%
%             'down'  :  Indicates that the function is in downstream Mach
%                        number input mode.
%
%             'totalp':  Indicates that the function is in total pressure
%                        ratio input mode.
%
%             'pito'  :  Indicates that the function is in Rayleigh-Pitot
%                        ratio input mode.
%
%   Outputs calculated for FLOWNORMALSHOCK are:
%   (All outputs are the same size as the array input or array inputs.  If
%   there are no array inputs, all outputs are scalars.)
%
%   MACH    :  An array of upstream Mach numbers.
%
%   T       :  An array of temperature ratios.  The temperature ratio is
%              defined as the static temperature downstream of the shock
%              over the static temperature upstream of the shock.
%
%   P       :  An array of pressure ratios.  The pressure ratio is defined
%              as the static pressure downstream of the shock over the
%              static pressure upstream of the shock.
%
%   RHO     :  An array of density ratios.  The density ratio is defined as
%              the density of the fluid downstream of the shock over the
%              density upstream of the shock.
%
%   M       :  An array of downstream Mach numbers.
%
%   P0      :  An array of total pressure ratios.  The total pressure ratio
%              is defined as the total pressure downstream of the shock
%              over the total pressure upstream of the shock.
%
%   P1      :  An array of Rayleigh-Pitot ratios.  The Rayleigh-Pitot ratio
%              is defined as the static pressure upstream of the shock over
%              the total pressure downstream of the shock.
%   
%   Examples:
%
%   Calculate the normal shock relations for air (gamma = 1.4) for total
%   pressure ratio of 0.61.  The following returns scalar values for MACH,
%   T, P, RHO, M, P0, and P1.
%
%      [MACH, T, P, RHO, M, P0, P1] = flownormalshock(1.4, 0.61, 'totalp')
%
%   Calculate the normal shock relations for gases with specific heat
%   ratios given in the following 1 x 4 row array for upstream Mach number
%   1.5.  The follow yields a 1 x 4 array for MACH, T, P, RHO, M, P0, and
%   P1.
%
%      gamma = [1.3, 1.33, 1.4, 1.67]
%      [MACH, T, P, RHO, M, P0, P1] = flownormalshock(gamma, 1.5)
%
%   Calculate the normal shock relations for a specific heat ratio of 1.4
%   and range of density ratios from 2.40 to 2.70 in increments of 0.10.
%   The following returns a 4 x 1 column array for MACH, T, P, RHO, M, P0,
%   and P1.
%
%      [MACH, T, P, RHO, M, P0, P1] = flownormalshock(1.4, (2.4:.1:2.7)', 'dens')
%
%   Calculate the normal shock relations for gases with specific heat ratio
%   and downstream Mach number combinations as shown.  The following
%   returns a 1 x 2 array for MACH, T, P, RHO, M, P0, and P1 each, where
%   the elements of each vector corresponds to the inputs element-wise.
%
%       gamma = [1.3, 1.4]
%       M = [.34, .49]
%       [MACH, T, P, RHO, M, P0, P1] = flownormalshock(gamma, M, 'down')
%       
%
%   See also FLOWISENTROPIC, FLOWPRANDTLMEYER, FLOWFANNO, FLOWRAYLEIGH.


%   Copyright 2009-2018 The MathWorks, Inc.

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
        error(message('aero:flownormalshock:paramSelectString'));
    end
end

% Check gamma

% Error if specific heat ratio input is not numeric
if ~isnumeric(gamma)
    error(message('aero:flownormalshock:notNumericGamma'));
end

% gamma > 1 check (specific heat ratio must be greater than 1)
if any(gamma<=1)
    error(message('aero:flownormalshock:gammaOneOrLess'));
end

% gamma real number check
if ~isreal(gamma)
    error(message('aero:flownormalshock:imaginaryGamma'));
end

% Error if gamma input is a matrix or if second input variable is a matrix
if (~isvector(gamma) || ~isvector(varargin{1}))
    error(message('aero:flownormalshock:multiDimensional'))
end

% Second input general checks (varargin{1})

% Error if second input variable is not numeric
if ~isnumeric(varargin{1})
    error(message('aero:flownormalshock:notNumeric'));
end

% Error if inputs are not the same size and if neither are scalars
if ~((isscalar(gamma) || isscalar(varargin{1})) || ...
        (isequal(size(gamma),size(varargin{1}))))
    error(message('aero:flownormalshock:size'))
end

% Second input variable real number check
if ~isreal(varargin{1})
    error(message('aero:flownormalshock:imaginary'));
end

% Initialize the string variable as the default (mach)
param = 'mach';

if nargin == 3 % If there is a selector string for the second input
    
    % Making string inputs more flexible
    
    % Typing the string 'mach' to input Mach number is optional
    if strcmpi(varargin{2},'mach')
        param = 'mach';

    % All the user really has to type in is 'temp' to input temperatureRatio    
    elseif strncmpi(varargin{2},'tempratio',4)
        param = 'tempratio';
        
    % All the user really has to type in is 'pres' to input pressureRatio
    elseif strncmpi(varargin{2},'pressratio',4)
        param = 'pressratio';
        
    % All the user really has to type in is 'dens' to input densityRatio
    elseif strncmpi(varargin{2},'densityratio',4)
        param = 'densityratio';
        
    % All the user really has to type in is 'down' to input
    % downstreamMachNumber
    elseif strncmpi(varargin{2},'downstreammach',4)
        param = 'downstreammach';
        
    % All the user really has to type in is 'totalp' to input
    % totalPressureRatio
    elseif strncmpi(varargin{2},'totalpressureratio',6)
        param = 'totalpressureratio';
        
    % All the user really has to type in is 'pito' to input
    % pitotRayleighRatio
    elseif strncmpi(varargin{2},'pitotrayleighratio',4)
        param = 'pitotrayleighratio';
        
    % If the user has a third variable and does not use the third variable
    % as an acceptable string then provide an error message
    else
        error(message('aero:flownormalshock:paramSelectWrongInput'))
    end
end

% Specific cases

switch lower(param)
    
    case 'mach'
        
        % Error if upstream Mach number input is less than unity
        if any(varargin{1}<1)
            error(message('aero:flownormalshock:lessThanUnity'));
        end
        
        % The second input (first variable input) is upstream Mach number
        mach = varargin{1};
        
    case 'tempratio'
        
        % Error if temperature ratio or specific heat ratio inputs are not
        % scalars for temperature ratio input mode
        if ~(isscalar(gamma) && isscalar(varargin{1}))
            error(message('aero:flownormalshock:tempRatioScalar'))
        end
        
        % Error if temperature ratio is less than 1
        if varargin{1}<1
            error(message('aero:flownormalshock:tempRatio'))
        end
        
        % Error if specific heat ratio input is not finite (for fzero fcn)
        if ~isfinite(gamma)
            error(message('aero:flownormalshock:tempRatioGammaFinite'))
        end
        
        % Error if temperature ratio input is Not-a-Number (for fzero fcn)
        if isnan(varargin{1})
            error(message('aero:flownormalshock:tempRatioNan'))
        end
        
        % The second input (first variable input) is temperature ratio T
        T = varargin{1};
        
        % Find Mach number through this function (below)
        if T==inf
            mach=inf;
        else
            mach = invertTempRatio(gamma, T);
        end
        
        
    case 'pressratio'
        
        % Error if pressure ratio is less than 1
        if any(varargin{1}<1)
            error(message('aero:flownormalshock:pressureRatio'));
        end
        
        % The second input (first variable input) is pressure ratio P
        P = varargin{1};
        
        % Find upstream Mach number from pressure ratio with this closed
        % form equation.
        mach = ((gamma+1)./(2.*gamma).*(P+(gamma-1)./(gamma+1))).^(1./2);
        
    case 'densityratio'
        
        upperLimit = (gamma+1)./(gamma-1);
        
        if ~isscalar(gamma)
            obLoDens = upperLimit(varargin{1}<1);
            obLoGamma = gamma(varargin{1}<1);
            obHiDens = upperLimit(varargin{1}>upperLimit);
            obHiGamma = gamma(varargin{1}>upperLimit);
        else
            obLoDens = upperLimit;
            obLoGamma = gamma;
            obHiDens = upperLimit;
            obHiGamma = gamma;
        end
        
        % Error if density ratio is less than 1
        if any(varargin{1}<1)
            error(message('aero:flownormalshock:densityRatio', sprintf( '%.17f', obLoDens( 1 ) ), sprintf( '%.15f', obLoGamma( 1 ) )))
        end
        
        % Error if density ratio is greater than (GAMMA+1)./(GAMMA-1)
        if any(varargin{1}>((gamma+1)./(gamma-1)))
            error(message('aero:flownormalshock:densityRatio', sprintf( '%.17f', obHiDens( 1 ) ), sprintf( '%.15f', obHiGamma( 1 ) )))
        end
        
        % The second input (first variable input) is density ratio
        rho = varargin{1};
        
        % Find upstream Mach number from density ratio with this closed
        % form equation.
        mach = (2.*rho./((gamma+1)-(gamma-1).*rho)).^(1./2);
        
    case 'downstreammach'
        
        lowerLimit = sqrt((gamma-1)./(2.*gamma));
        
        if ~isscalar(gamma)
            obHiDown = lowerLimit(varargin{1}>1);
            obHiGamma = gamma(varargin{1}>1);
            obLoDown = lowerLimit(varargin{1}<lowerLimit);
            obLoGamma = gamma(varargin{1}<lowerLimit);
        else
            obHiDown = lowerLimit;
            obHiGamma = gamma;
            obLoDown = lowerLimit;
            obLoGamma = gamma;
        end
        
        % Error if downstream Mach number is greater than 1
        if any(varargin{1}>1)
            error(message('aero:flownormalshock:downstreamMach', sprintf( '%.17f', obHiDown( 1 ) ), sprintf( '%.15f', obHiGamma( 1 ) )))
        end
        
        % Error if downstream Mach number is less than
        % SQRT((GAMMA-1)./(2.*GAMMA))
        if any(varargin{1}<(sqrt((gamma-1)./(2.*gamma))))
            error(message('aero:flownormalshock:downstreamMach', sprintf( '%.17f', obLoDown( 1 ) ), sprintf( '%.15f', obLoGamma( 1 ) )))
        end
        
        % The second input (first variable input) is downstream Mach number
        M = varargin{1};
        
        % Find upstream Mach number from downstream Mach number with this
        % closed form equation., or force mach number to be infinity if abs of
        % equation denominator is eps or smaller
        mach = ((M.^2+2./(gamma-1))./(2.*gamma./(gamma-1).*M.^2-1)).^(1./2);
        
        mach_den = (2.*gamma./(gamma-1).*M.^2-1);
        mach(abs(mach_den)<eps) = inf;
        
        
    case 'totalpressureratio'
        
        % Error if total pressure ratio or specific heat ratio inputs are
        % not scalars for total pressure ratio input mode
        if ~(isscalar(gamma) && isscalar(varargin{1}))
            error(message('aero:flownormalshock:totalPressureRatioScalar'))
        end
        
        % Error if the total pressure ratio is negative
        if varargin{1}<0
            error(message('aero:flownormalshock:totalPressureRatio'))
        end
        
        % Error if the total pressure ratio is greater than 1
        if varargin{1}>1
            error(message('aero:flownormalshock:totalPressureRatio'))
        end
        
        % Error if either input is not finite in total pressure ratio mode
        % (for fzero fcn)
        if ~(isfinite(gamma) && isfinite(varargin{1}))
            error(message('aero:flownormalshock:totalPressureRatioFinite'))
        end
        
        % The second input (first variable input) is total pressure ratio
        P0 = varargin{1};
        
        % Find Mach number through this function (below)
        if P0 == 0
            mach = inf;
        else
            mach = invertTotalPressureRatio(gamma, P0);
        end
        
        
    case 'pitotrayleighratio'
        
        % Error if total pressure ratio or specific heat ratio inputs are
        % not scalars for Rayleigh-Pitot ratio input mode
        if ~(isscalar(gamma) && isscalar(varargin{1}))
            error(message('aero:flownormalshock:pitotRayleighRatioScalar'))
        end
        
        upperLimit = ((gamma+1)/2)^(-gamma/(gamma-1));
        
        obLoPitot = upperLimit(varargin{1}<0);
        obLoGamma = gamma(varargin{1}<0);
        obHiPitot = upperLimit(varargin{1}>upperLimit);
        obHiGamma = gamma(varargin{1}>upperLimit);
        
        % Error if Rayleigh-Pitot ratio is less than 0
        if varargin{1}<0
            error(message('aero:flownormalshock:pitotRayleighRatio', sprintf( '%.17f', obLoPitot ), sprintf( '%.15f', obLoGamma )))
        end
        
        % Error if Rayleigh-Pitot ratio is greater than
        % ((GAMMA+1)./2).^(-GAMMA./(GAMMA-1))
        if varargin{1}>((gamma+1)/2)^(-gamma/(gamma-1))
            error(message('aero:flownormalshock:pitotRayleighRatio', sprintf( '%.17f', obHiPitot ), sprintf( '%.15f', obHiGamma )))
        end    
        
        % Error if either input is not finite in Rayleigh-Pitot ratio mode
        % (for fzero fcn)
        if ~(isfinite(gamma) && isfinite(varargin{1}))
            error(message('aero:flownormalshock:pitotRayleighRatioFinite'))
        end
        
        % The second input (first variable input) is Rayleigh-Pitot ratio
        P1 = varargin{1};
        
        % Find Mach number through this function (below)
        if isequal(P1,0)
            mach=inf;
        elseif isequal(P1,((gamma+1)/2)^(-gamma/(gamma-1)))
            mach=1;
        else
            mach = invertRayleighPitotRatio(gamma, P1);
        end
        
end

% Function main calculations (after Mach number is given or deduced)

% Static temperature ration downstream over upstream, T = Temp_2/Temp_1
T=mach.^2.*(mach.^-2+(gamma-1)./2).*(2.*gamma./(gamma-1)-mach.^-2)./...
    (2*gamma./(gamma-1)+(gamma-1)./2);

% Static pressure ratio downstream over upstream, P = p_2/p_1
P=2.*gamma.*mach.^2./(gamma+1)-(gamma-1)./(gamma+1);

% Static density ratio downstream over upstream, rho = rho_2/rho_1
rho=(gamma+1)./((gamma-1)+2.*mach.^-2);

% Downstream Mach number M
M=((gamma-1+2.*mach.^-2)./(2.*gamma-mach.^-2.*(gamma-1))).^(1/2);

% Total pressure ratio downstream over upstream, P0 = p_02/p_01
P0=((gamma+1)./(2.*mach.^-2+(gamma-1))).^...
    (gamma./(gamma-1)).*(1./(2.*gamma./(gamma+1).*mach.^2-(gamma-1)./...
    (gamma+1))).^(1./(gamma-1));

% Ratio of upstream static pressure over downstream total pressure,
% P1 = p1/p02 (also called the Rayleigh-Pitot relation)
P1=(2.*mach.^-2./(gamma+1)).*((4.*gamma-2.*(gamma-1).*mach.^-2)./...
    ((gamma+1).^2)).^(1./(gamma-1));

% Return the values input by the user, make the size such that all outputs
% are the same.

switch param
    case 'mach'
        mach = varargin{1}.*ones(size(gamma));
    case 'tempratio'
        T = varargin{1}.*ones(size(mach));
    case 'pressratio'
        P = varargin{1}.*ones(size(mach));
    case 'densityratio'
        rho = varargin{1}.*ones(size(mach));
    case 'downstreammach'
        M = varargin{1}.*ones(size(mach));
    case 'totalpressureratio'
        P0 = varargin{1}.*ones(size(mach));
    case 'pitotrayleighratio'
        P1 = varargin{1}.*ones(size(mach));
end


%--------------------------------------------------------------------------
function mach = invertTempRatio(gamma, T)
% This function uses a residual method of solving for the Mach number when
% given temperature ratio because a closed form solution for Mach number of
% the temperature ratio equation cannot be found.  This method involves
% solving the following residual equation for zeros:
% fcn(tempRatio) = (known) value of tempRatio - equation for tempRatio
%                                               as a fcn(GAMMA,MACH)
%
% The Mach number is being sought so it must remain anonymous until its
% value has been found by the fzero function.  The quality of the solution
% is critically dependent on the quality of the initial guess in the fzero
% function.


% This equation is fcn(tempRatio) as described above
equation = @(gamma,T,mach) T-(1+(gamma-1)./2.*mach.^2).*(2*gamma./(gamma-1).*...
    mach.^2-1)./(mach.^2.*(2*gamma./(gamma-1)+(gamma-1)./2));

% Initial guess of Mach number for fzero fcn
    guessMach = [1 realmax^(1/5)];

% Estimate zeros of the residual function
machEval = fzero(@(mach) equation(gamma,T,mach), guessMach);

% Let the Mach number to be used be the same as evaluated Mach number
mach = machEval;

%--------------------------------------------------------------------------
function mach = invertTotalPressureRatio(gamma, P0)
% This function uses a residual method of solving for the Mach number when
% given total pressure ratio because a closed form solution for Mach number
% of the total pressure ratio equation cannot be found.  This method
% involves solving the following residual equation for zeros:
% fcn(totalPressureRatio) = (known) value of totalPressureRatio - 
%                      equation for totalPressureRatio as a fcn(GAMMA,MACH)
%
% The Mach number is being sought so it must remain anonymous until its
% value has been found by the fzero function.  The quality of the solution
% is critically dependent on the quality of the initial guess in the fzero
% function.


% This equation is fcn(totalPressureRatio) as described above
equation = @(gamma,P0,mach) P0-((gamma+1)./2.*mach.^2./...
    (1+(gamma-1)./2.*mach.^2)).^(gamma./(gamma-1)).*...
    (1./(2*gamma./(gamma+1).*mach.^2-(gamma-1)./(gamma+1))).^(1./(gamma-1));

% Initial guess of Mach number for fzero fcn

% ONLY GOOD FROM GAMMA = 1.01 TO GAMMA = 18
    guessMach = [1 sqrt(realmax)^((gamma-1)/gamma)];


% Estimate zeros of the residual function
machEval = fzero(@(mach) equation(gamma,P0,mach), guessMach);

% Let the Mach number to be used be the same as evaluated Mach number
mach = machEval;

%--------------------------------------------------------------------------
function mach = invertRayleighPitotRatio(gamma, P1)
% This function uses a residual method of solving for the Mach number when
% given Rayleigh-Pitot ratio because a closed form solution for Mach number
% of the Rayleigh-Pitot ratio equation cannot be found.  This method
% involves solving the following residual equation for zeros:
% fcn(pitotRayleighRatio) = (known) value of pitotRayleighRatio - 
%                      equation for pitotRayleighRatio as a fcn(GAMMA,MACH)
%
% The Mach number is being sought so it must remain anonymous until its
% value has been found by the fzero function.  The quality of the solution
% is critically dependent on the quality of the initial guess in the fzero
% function.


% This equation is fcn(pitotRayleighRatio) as described above
equation = @(gamma,P1,mach) P1-(1+(gamma-1)./2.*mach.^2).^(-gamma./(gamma-1)).*...
    (((gamma+1)./2.*mach.^2./(1+(gamma-1)./2.*mach.^2)).^(gamma./(gamma-1)).*...
    (1./(2*gamma./(gamma+1).*mach.^2-(gamma-1)./(gamma+1))).^(1./(gamma-1))).^-1;

% Initial guess of Mach number for fzero fcn

% ONLY GOOD FROM GAMMA = 1.03 TO GAMMA = 100
    guessMach = [1 sqrt(realmax)^((gamma-1)/gamma)];


% Estimate zeros of the residual function
machEval = fzero(@(mach) equation(gamma,P1,mach), guessMach);

% Let the Mach number to be used be the same as evaluated Mach number
mach = machEval;
