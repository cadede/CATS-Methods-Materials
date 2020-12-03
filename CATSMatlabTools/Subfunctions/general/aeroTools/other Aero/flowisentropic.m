function [mach, T, P, rho, A] = flowisentropic(gamma, varargin)
%FLOWISENTROPIC Calculate isentropic flow ratios
%   [MACH, T, P, RHO, A] = FLOWISENTROPIC(GAMMA, VAR, MTYPE) returns an
%   array of isentropic flow Mach number, MACH, temperature ratio, T,
%   pressure ratio, P, density ratio, RHO, and area ratio, A.
%   FLOWISENTROPIC calculates these arrays given a set of specific heat
%   ratios, GAMMA, and any one of the isentropic flow variables.  The
%   isentropic flow variable is selected by the string, MTYPE. The
%   temperature, pressure, and density ratios are a comparison of the local
%   static conditions over the stagnation (or total) conditions.  The area
%   ratio is a comparison of the instantaneous streamtube area for a
%   quasi-one-dimensional flow over the reference area (throat area or
%   minimum area) where the Mach number of the flow becomes unity.
%
%   See documentation for more information about assumptions and
%   limitations.
%
%   Inputs for FLOWISENTROPIC are:
%
%   GAMMA :  An array of N specific heat ratios.  GAMMA Must be a scalar
%            or an array of N of real numbers greater than 1.  For
%            subsonic area ratio input mode and supersonic area ratio
%            input mode, GAMMA must be a real, finite scalar greater than
%            1.
%
%   VAR   :  An array of real numerical values for one of the isentropic
%            flow relations.
%
%            MACH :  An array of Mach numbers.  Must be a scalar or array
%                    of N real numbers greater than or equal to 0.  If MACH
%                    and GAMMA are both arrays, they must be the same size.
%                    MACH is used with MTYPE variable 'mach'.  Since 'mach'
%                    is the default MTYPE value, the MTYPE variable is
%                    optional for Mach number input mode.
%
%            T    :  An array of temperature ratios.  The temperature ratio
%                    is defined as the local static temperature over the
%                    stagnation temperature.  T must be a scalar or an
%                    array of real numbers greater than or equal to 0 (as
%                    MACH -> inf) and less than or equal to 1 (at 
%                    MACH = 0).  If T and GAMMA are both arrays, they must
%                    be the same size.  T is used with MTYPE variable
%                    'temp'.
%
%            P    :  An array of pressure ratios.  The pressure ratio is
%                    defined as the local static pressure over the
%                    stagnation pressure.  P must be a scalar or an array
%                    of real numbers greater than or equal to 0 (as 
%                    MACH -> inf) and less than or equal to 1 (at 
%                    MACH = 0).  If P and GAMMA are both arrays, they must
%                    be the same size.  P is used with MTYPE variable
%                    'pres'.
%
%            RHO  :  An array of density ratios.  The density ratio is
%                    defined as the local density over the stagnation
%                    density.  RHO must be a scalar or an array of real
%                    numbers greater than or equal to 0 (as MACH -> inf)
%                    and less than or equal to 1 (at MACH = 0).  If RHO and
%                    GAMMA are both arrays, they must be the same size.
%                    RHO is used with MTYPE variable 'dens'.
%
%            A    :  A scalar value of area ratio.  Must be a real value
%                    A >= 1.  A is used with MTYPE variables 'sub' or
%                    'sup'.
%
%   MTYPE :  A character vector for selecting the isentropic flow variable
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
%             'sub'   :  Indicates that the function is in subsonic area
%                        ratio input mode.  The subsonic area ratio is
%                        defined as the local subsonic streamtube area over
%                        the reference streamtube area for sonic
%                        conditions.
%
%             'sup'   :  Indicates that the function is in supersonic area
%                        ratio input mode.  The supersonic area ratio is
%                        defined as the local supersonic streamtube area
%                        over the reference streamtube area for sonic
%                        conditions.
%
%   Outputs calculated for FLOWISENTROPIC are:
%   (All outputs are the same size as the array input or array inputs.  
%   If there are no array inputs, all outputs are scalars)
%
%   MACH    :  An array of Mach numbers.
%
%   T       :  An array of temperature ratios.  The temperature ratio is
%              defined as the local static temperature over the stagnation
%              temperature.
%
%   P       :  An array of pressure ratios.  The pressure ratio is defined
%              as the local static pressure over the stagnation pressure.
%
%   RHO     :  An array of density ratios.  The density ratio is defined as
%              the local density over the stagnation density.
%
%   A       :  An array of area ratios.  The area ratio is defined as the
%              local streamtube area over the reference streamtube area for
%              sonic conditions.
%
%   Examples:
%
%   Calculate the isentropic flow relations for air (gamma = 1.4) for a
%   design subsonic area ratio of 1.255.  The following returns scalar
%   values for MACH, T, P, RHO, and A.
%
%      [MACH, T, P, RHO, A] = flowisentropic(1.4, 1.255, 'sub')
%
%   Calculate the isentropic flow relations for gases with specific heat
%   ratios given in the following 1 x 4 row array for the Mach number 0.5.
%   The following returns a 1 x 4 row array for MACH, T, P, RHO, and A.
%
%      gamma = [1.3, 1.33, 1.4, 1.67]
%      [MACH, T, P, RHO, A] = flowisentropic(gamma, 0.5)
%
%   Calculate the isentropic flow relations for a specific heat ratio of
%   1.4  and range of temperature ratios from 0.40 to 0.70 in increments of
%   0.10 The following returns a 4 x 1 column array for MACH, T, P, RHO,
%   and A.
%
%      [MACH, T, P, RHO, A] = flowisentropic(1.4, (0.40:0.10:0.70)', 'temp')
%
%   Calculate the isentropic flow relations for gases with specific heat
%   ratio and density ratio combinations as shown.  The following returns a
%   1 x 2 array for MACH, T, P, RHO, and A each, where the elements of each
%   vector correspond to the inputs element-wise.
%
%       gamma = [1.3, 1.4]
%       RHO = [0.13, 0.9]
%       [MACH, T, P, RHO, A] = flowisentropic(gamma, RHO , 'dens')
%
%
%   See also FLOWNORMALSHOCK, FLOWPRANDTLMEYER, FLOWFANNO, FLOWRAYLEIGH.

%   Copyright 2009-2016 The MathWorks, Inc.

%   References:  James, J. E. A., Gas Dynamics, Second Edition, Allyn and
%   Bacon, Inc, Boston, 1984.
%
%   NACA Technical Report 1135, 1953, National Advisory Committee on
%   Aeronautics, Ames Research Staff, Moffett Field, Calif.

% Checks

% General checks

% Error in case of less than 2 inputs or when there are too many inputs
narginchk(2,3);

% Error if third input exists and is not a string
if nargin == 3
    if ~ischar(varargin{2}) && ~isstring(varargin{2})
        error(message('aero:flowisentropic:paramSelectString'))
    end
end

% Check gamma

% Error if specific heat ratio input is not numeric
if ~isnumeric(gamma)
    error(message('aero:flowisentropic:notNumericGamma'));
end

% gamma > 1 check (specific heat ratio must be greater than 1)
if any(gamma<=1)
    error(message('aero:flowisentropic:gammaOneOrLess'));
end

% gamma real number check
if ~isreal(gamma)
    error(message('aero:flowisentropic:imaginaryGamma'));
end

% Error if gamma input is a matrix or if second input variable is a matrix
if (~isvector(gamma) || ~isvector(varargin{1}))
    error(message('aero:flowisentropic:multiDimensional'))
end

% Second input general checks (varargin{1})

% Error if second input variable is not numeric
if ~isnumeric(varargin{1})
    error(message('aero:flowisentropic:notNumeric'));
end


% Error if inputs are not the same size and if neither are scalars
if ~((isscalar(gamma) || isscalar(varargin{1})) || ...
        (isequal(size(gamma),size(varargin{1}))))
    error(message('aero:flowisentropic:size'))
end

% Second input variable real number check
if ~isreal(varargin{1})
    error(message('aero:flowisentropic:imaginary'));
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
        
    % All the user really has to type in is 'sub' to input subSonicAreaRatio
    elseif strncmpi(varargin{2},'subsonicarearatio',3)
        param = 'subsonicarearatio';
        
    % All the user really has to type in is 'sup' to input
    % superSonicAreaRatio
    elseif strncmpi(varargin{2},'supersonicarearatio',3)
        param = 'supersonicarearatio';
        
    % If the user has a third variable and does not use the third variable
    % as an acceptable string then provide an error message
    else
        error(message('aero:flowisentropic:paramSelectWrongInput'))
    end
end

% Specific cases

switch lower(param)
    
    case 'mach'
        
        % Error if Mach number input is negative
        if any(varargin{1}<0)
            error(message('aero:flowisentropic:negative'));
        end
        
        % The second input (first variable input) is Mach number.
        mach = varargin{1};
        
    case 'tempratio'
        
        % Error if temperature ratio is less than 0
        if any(varargin{1}<0)
            error(message('aero:flowisentropic:temperatureRatio'));
        end
        
        % Error if temperature ratio is greater than 1
        if any(varargin{1}>1)
            error(message('aero:flowisentropic:temperatureRatio'));
        end
        
        % The second input (first variable input) is temperature ratio T
        T = varargin{1};
        
        % Find Mach number from temperature ratio with this closed form equation.
        mach = (2./(gamma-1).*(T.^-1-1)).^(1/2);
        
    case 'pressratio'
        
        % Error if pressure ratio is less than 0
        if any(varargin{1}<0)
            error(message('aero:flowisentropic:pressureRatio'));
        end
        
        % Error if pressure ratio is greater than 1
        if any(varargin{1}>1)
            error(message('aero:flowisentropic:pressureRatio'));
        end
        
        % The second input (first variable input) is pressure ratio P
        P = varargin{1};
        
        % Find Mach number from pressure ratio with this closed form equation.
        mach = (2./(gamma-1).*(P.^((1-gamma)./gamma)-1)).^(1/2);
        
    case 'densityratio'
        
        % Error if density ratio is less than 0
        if any(varargin{1}<0)
            error(message('aero:flowisentropic:densityRatio'));
        end
        
        % Error if density ratio is greater than 1
        if any(varargin{1}>1)
            error(message('aero:flowisentropic:densityRatio'));
        end
        
        % The second input (first variable input) is density ratio RHO
        rho = varargin{1};
        
        % Find Mach number from density ratio with this closed form equation.
        mach = (2./(gamma-1).*(rho.^(1-gamma)-1)).^(1/2);
        
    case 'subsonicarearatio'
        
        % The second input (first variable input) is area ratio A
        A=varargin{1};
        
        % Find area ratio errors with this function (below)
        areaRatioSpecificErrors(gamma,A);
        
        % Find Mach number through this function (below)
        if isequal(A,inf)
            mach = 0;
        else
            mach = machFromAreaInput(gamma, A, param, 0);
        end
        
    case 'supersonicarearatio'
        
        % The second input (first variable input) is area ratio A
        A=varargin{1};
        
        % Find area ratio errors with this function (below)
        areaRatioSpecificErrors(gamma,A);
        
        % Find Mach number through this function (below)
        mach = machFromAreaInput(gamma, A, param, inf);
        
end

% Function main calculations (after Mach number is given or deduced)

% Static temperature over stagnation temperature, T = Temp/Temp_0
T=(1+(gamma-1)./2.*mach.^2).^-1;

% Static pressure over stagnation pressure, P = p/p_0
P=((1+(gamma-1)./2.*mach.^2).^(gamma./(gamma-1))).^-1;

% Static density over stagnation density, rho = rho_static/rho_0
rho=(1+(gamma-1)./2.*mach.^2).^(-1./(gamma-1));


% Area ratio at a point over area required for, M = 1, A = Area/Area*
b = (gamma+1)./(2.*(1-gamma)); % power in area ratio calculations

A=((gamma+1)./2).^b./(mach.*(1+(gamma-1)./2.*mach.^2).^b);

% Special case: Use this rearranged equation for Area Ratio if Mach number
% is infinity (otherwise evaluated as NaN)
if isequal(size(gamma),size(mach))
    A(mach==inf)=((gamma(mach==inf)+1)./2).^b(mach==inf).*(mach(mach==inf).^-2+...
    (gamma(mach==inf)-1)./2).^-b(mach==inf).*mach(mach==inf).^(-2.*b(mach==inf)-1);
elseif (isscalar(gamma) && ~isscalar(mach))
    A(mach==inf)=((gamma+1)./2).^b.*(mach(mach==inf).^-2+...
    (gamma-1)./2).^-b.*mach(mach==inf).^(-2.*b-1);
elseif (~isscalar(gamma) && isscalar(mach))
    if mach==inf
    A=((gamma+1)./2).^b.*(mach.^-2+(gamma-1)./2).^-b.*mach.^(-2.*b-1);
    end
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
    case 'densityratio'
        rho = varargin{1}.*ones(size(mach));
    case {'subsonicarearatio','supersonicarearatio'}
        A = varargin{1}.*ones(size(mach));
end

%-------------------------------------------------------------------------
function mach = invertAreaRatio(gamma, A, param)
% This function uses a residual method of solving for the Mach number when
% given area ratio because a closed form solution of the area ratio
% equation for Mach number cannot be found.  This method involves solving
% the following residual equation for zeros:
% fcn(areaRatio) = (known) value of areaRatio - equation for areaRatio
%                                               as a fcn(GAMMA,MACH)
%
% The Mach number is being sought so it must remain anonymous until its
% value has been found by the fzero function.  The quality of the solution
% is critically dependent on the quality of the initial guess in the fzero
% function because for each area ratio there is both a subsonic and
% supersonic solution.  The solution that is found corresponds to the
% user's choice of the subsonic or supersonic value from the main function
% flowisentropic above.

b = (gamma+1)./(2.*(1-gamma)); % power in area ratio calculations

% This equation is fcn(areaRatio) as described above
equation = @(gamma,A,mach) A-((gamma+1)./2).^b./(mach.*(1+(gamma-1)./2.*mach.^2).^b);

if strcmpi(param,'subsonicarearatio')
    
    % Subsonic guess
    guessMach = [realmin 1];
    
elseif strcmpi(param,'supersonicarearatio')
    
    % Supersonic guess
    if gamma < 1.1
        % mach > 37 calculates A* = Inf
        guessMach = [1 realmax/4.76e306];
    elseif gamma <1.2
        % mach > 1.7977e+13 calculates A* = inf
        guessMach = [1 realmax/1e295];
    elseif gamma < 1.3
        % mach > 1.7977e+29 calculates A* = inf
        guessMach = [1 realmax/1e279];
    elseif gamma < 1.4
        % mach > 1.7977e+42 calculates A* = inf
        guessMach = [1 realmax/1e266];
    elseif gamma < 1.8
        % mach > 1.7977e+54 calculates A* = inf
        guessMach = [1 realmax/1e254];
    elseif gamma < 1.9
        % mach > 1.7977e+93 calculates A* = inf
        guessMach = [1 realmax/1e216];
    elseif gamma < 2.0
        % mach > 1.7977e+100 calculates A* = inf
        guessMach = [1 realmax/1e208];
    else
        % mach > 1.7977e+107 calculates A* = inf
        guessMach = [1 realmax/1e201];
    end
end

% Evaluate the Mach number which corresponds to the zero crossing of the
% fcn(areaRatio) equation.
machEval = fzero(@(mach) equation(gamma,A,mach), guessMach);

% Let the Mach number to be used be the same as evaluated Mach number
mach = machEval;

%--------------------------------------------------------------------------
function mach = machFromAreaInput(gamma, A, param, K)
% This function finds the Mach number when an area ratio is the second
% input in the function flowisentropic above.  The subfunction
% invertAreaRatio has the method for solving for the Mach number.

% Find the Mach number when area ratio is finite
if isfinite(A)
    mach = invertAreaRatio(gamma,A,param);
    
elseif ~isnan(A) % If area ratio input is inf
    mach = K; % K = 0 (subsonic) or K = inf (supersonic)
end

%--------------------------------------------------------------------------
function areaRatioSpecificErrors(gamma, A)
% This function catches errors specific to area ratio.  These errors are in
% a function here because the error messages and error IDs are not unique
% between subsonic and supersonic area ratios.

% Error if area ratio or specific heat ratio inputs are not scalars for
% area ratio input mode
if ~(isscalar(gamma) && isscalar(A))
    error(message('aero:flowisentropic:areaRatioScalar'))
end

% Error if area ratio is less than 1
if any(A<1)
    error(message('aero:flowisentropic:areaRatio'))
end

% Error if area ratio is Not-a-Number
if isnan(A)
    error(message('aero:flowisentropic:nan'))
end

% Error if specific heat ratio is not a finite value
if ~isfinite(gamma)
    error(message('aero:flowisentropic:areaRatioGamma'))
end
