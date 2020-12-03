function out = correctairspeed( in, a, P0, ain, aout, varargin )
%  CORRECTAIRSPEED Calculate equivalent airspeed (EAS), calibrated airspeed
%  (CAS) or true airspeed (TAS) from one of the other two airspeeds.
%
%   AS = CORRECTAIRSPEED( V, A, P0, AI, AO ) computes the conversion factor
%   from specified input airspeed, AI, to specified output airspeed, AO,
%   using speed of sound, A, and static pressure P0. The conversion factor
%   is applied to the input airspeed, V, to produce the output, AS, in the
%   desired airspeed. V, AS, A, and P0 are floating-point arrays of size M.
%   All of the values in V must have the same airspeed conversions from AI
%   to AO. AI and AO are strings.
%
%   AS = CORRECTAIRSPEED( V, A, P0, AI, AO, METHOD ) uses METHOD to compute
%   the conversion factor from specified input airspeed, AI, to specified
%   output airspeed, AO, using speed of sound, A, and static pressure P0.
%   METHOD is a string.
%
%   Inputs to CORRECTAIRSPEED are:
%      - V: airspeed in meters per second 
%      - A: speed of sound in meters per second
%      - P0: static air pressure in pascal
%      - AI: input airspeed string
%      - AO: output airspeed string
%      - METHOD (Optional): conversion method string
%
%   Supported strings for AI and AO are:
%      - 'TAS': true airspeed    
%      - 'CAS': calibrated airspeed
%      - 'EAS': equivalent airspeed
%
%   Supported strings for METHOD are:
%      - 'TableLookup' (Default): Generate output AS by looking up or
%        estimating table values based on inputs V, A, and P0.
%      - 'Equation': Compute output AS directly using input values V, A,
%        and P0. Calculations involving supersonic airspeeds require an
%        iterative computation. If a solution is not achieved within 30
%        iterations an error message will be displayed.
%
%   Output calculated is: 
%      - AS: airspeed in meters per second
%
%   Limitations: 
%      - This function assumes compressible dry air with constant specific
%        heat ratio (gamma).
%      - Due to potentially high errors, METHOD 'TableLookup' is not
%        recommended when either of the following is true:
%         - Speed of sound A is less than 200 m/s or greater than 355 m/s.
%         - Static Pressure P0 is less than 1000 Pa or greater than
%           106500 Pa.
%      - Input METHOD is ignored when any of the following criteria is met.
%        In these cases, the function automatically uses METHOD 'Equation':
%         - Conversion from ain 'TAS' to aout 'EAS'.
%         - Conversion from ain 'EAS' to aout 'TAS'.
%         - Conversion where input V is greater than five times the speed
%           of sound at sea level (approximately 1700 m/s).
%
%   Examples:
%
%   Convert three airspeeds from true airspeed to calibrated airspeed at
%   1000 meters using method 'TableLookup':
%      ain = [25.7222; 10.2889; 3.0867];
%      as = correctairspeed(ain, 336.4, 89874.6, 'TAS', 'CAS', 'TableLookup')
%
%   Convert airspeeds from calibrated airspeed to equivalent airspeed at
%   1000 and 0 meters using method 'Equation': 
%      ain = [25.7222; 10.2889; 3.0867];
%      sos = [336.4; 340.3; 340.3];
%      P0 = [ 89874.6; 101325; 101325];
%      as = correctairspeed( ain, sos, P0, 'CAS', 'EAS', 'Equation')
%
%   Convert airspeed from true airspeed to equivalent airspeed at
%   15000 meters: 
%      ain = 376.25;
%      [~, sos, P0, ~] = atmoscoesa(15000);
%      as = correctairspeed( ain, sos, P0, 'TAS', 'EAS')
%
%   See also AIRSPEED, ATMOSCOESA, ATMOSISA, ATMOSLAPSE, ATMOSNONSTD.

%   Copyright 2000-2018 The MathWorks, Inc.

%   References:
%   [1] Lowry, J. T., Performance of Light Aircraft, AIAA Education Series,
%   Washington, DC, 1999.  
%   [2] Aeronautical Vestpocket Handbook, United Technologies Pratt &
%   Whitney, August, 1986. 
%   [3] Gracey, William, "Measurement of Aircraft Speed and Altitude",
%   NASA Reference Publication 1046, 1980.

% Validate Inputs
narginchk(1, 6);
if nargin < 6
    Method = 'TableLookup';
else
    Method = varargin{1};
    validatestring(Method, {'TableLookup', 'Equation'}, '','Method');
end

if ~isfloat( in )
    error(message('aero:correctairspeed:notFloat1'));
end
validateattributes(in,{'numeric'},{'real', 'nonnegative'},'','V');

if ~isfloat( a )
    error(message('aero:correctairspeed:notFloat2'));
end
validateattributes(a,{'numeric'},{'real', 'nonnegative'},'','A');

if ~isfloat( P0 )
    error(message('aero:correctairspeed:notFloat3'));
end
validateattributes(P0,{'numeric'},{'real', 'nonnegative'},'','P0');

if ~ischar( ain ) && ~isstring( ain )
    error(message('aero:correctairspeed:notChar1'));
end

if ~ischar( aout ) && ~isstring( aout )
    error(message('aero:correctairspeed:notChar2'));
end

% Validate Input Sizes
if ~isequal(size(in), size(a), size(P0))
    maxSize = max([size(in); size(a); size(P0)]);
    if numel(in) == 1
        in = ones(maxSize)*in;
    end
    if numel(a) == 1
        a = ones(maxSize)*a;
    end
    if numel(P0) == 1
        P0 = ones(maxSize)*P0;
    end
    if ~isequal(size(in), size(a), size(P0))
        error(message('aero:correctairspeed:wrongDimension'));
    end
end

% If no correction is requested, return input velocity
if strcmpi(ain, aout)
    out = in;
    return;
end

% Reformat inputs to allow for vectorized computation
sizeIn = size(in);
in = in(:);
P0 = P0(:);
a = a(:);

% Standard Environment Parameters
gamma = 1.4;
a_sl = 340.2941;
P_sl = 101325;

% Compute sqrt of the density ratio
sqrtSigma = a_sl./a.*sqrt(rrdelta(P0, 0, gamma));

% Compute Corrected Airspeed
out = nan(size(in));
switch lower(ain)
    case 'tas'
        switch lower(aout)
            case 'cas'
                if strcmpi(Method, 'Equation')
                    q_c = arrayfun(@(v, p, sos) asToQc(v, p, sos, gamma), in, P0, a);
                    out = arrayfun(@(q) qcToAS(q, P_sl, a_sl, gamma), q_c);
                else
                    m_sl = in ./ a_sl;
                    % M_sl < 5
                    for machThres = 0 : 4
                        if ~isempty(m_sl(m_sl >= machThres & m_sl < machThres+1))
                            out(m_sl >= machThres & m_sl < machThres+1) = asLookup(in(m_sl >= machThres & m_sl < machThres+1),...
                                P0(m_sl >= machThres & m_sl < machThres+1), a(m_sl >= machThres & m_sl < machThres+1), machThres, 'c');
                        end
                    end
                    % M_sl > 5
                    if any(m_sl >= 5)
                        q_c = arrayfun(@(v, p, sos) asToQc(v, p, sos, gamma), in(m_sl >= 5), P0(m_sl >= 5), a(m_sl >= 5));
                        out(m_sl >= 5) = arrayfun(@(q) qcToAS(q, P_sl, a_sl, gamma), q_c);
                    end
                end
            case 'eas'
                out = in.*sqrtSigma;
            otherwise
                error(message('aero:correctairspeed:unknownAirspeed'));
        end
    case 'cas'
        tas = nan(size(in));
        if strcmpi(Method, 'Equation')
            q_c = arrayfun(@(v) asToQc(v, P_sl, a_sl, gamma), in);
            tas = arrayfun(@(q, p, sos, guess) qcToAS(q, p, sos, gamma), q_c, P0, a);
        else
            m_sl = in ./ a_sl;
            % M_sl < 5
            for machThres = 0 : 4
                if ~isempty(m_sl(m_sl >= machThres & m_sl < machThres+1))
                    tas(m_sl >= machThres & m_sl < machThres+1) = asLookup(in(m_sl >= machThres & m_sl < machThres+1),...
                        P0(m_sl >= machThres & m_sl < machThres+1), a(m_sl >= machThres & m_sl < machThres+1), machThres, 't');
                end
            end
            % M_sl > 5
            if any(m_sl >= 5)
                q_c = arrayfun(@(v) asToQc(v, P_sl, a_sl, gamma), in(m_sl >= 5));
                tas(m_sl >= 5) = arrayfun(@(q, p, sos, guess) qcToAS(q, p, sos, gamma), q_c, P0(m_sl >= 5), a(m_sl >= 5));
            end
        end
        switch lower(aout)
            case 'tas'
                out = tas;
            case 'eas'
                out = tas.*sqrtSigma;
            otherwise
                error(message('aero:correctairspeed:unknownAirspeed'));
        end
    case 'eas'
        switch lower(aout)
            case 'tas'
                out = in./sqrtSigma;
            case 'cas'
                if strcmpi(Method, 'Equation')
                    q_c = arrayfun(@(v, p, sos) asToQc(v, p, sos, gamma), in./ sqrtSigma, P0, a);
                    out = arrayfun(@(q) qcToAS(q, P_sl, a_sl, gamma), q_c);
                else
                    m_sl = in./sqrtSigma ./ a_sl;
                    % M_sl < 5
                    for machThres = 0 : 4
                        if ~isempty(m_sl(m_sl >= machThres & m_sl < machThres+1))
                            out(m_sl >= machThres & m_sl < machThres+1) =...
                                asLookup(in(m_sl >= machThres & m_sl < machThres+1)./ sqrtSigma(m_sl >= machThres & m_sl < machThres+1),...
                                P0(m_sl >= machThres & m_sl < machThres+1), a(m_sl >= machThres & m_sl < machThres+1), machThres, 'c');
                        end
                    end
                    % M_sl > 5
                    if any(m_sl >= 5)
                        q_c = arrayfun(@(v, p, sos) asToQc(v, p, sos, gamma), in(m_sl >= 5)./ sqrtSigma(m_sl >= 5), P0(m_sl >= 5), a(m_sl >= 5));
                        out(m_sl >= 5) = arrayfun(@(q) qcToAS(q, P_sl, a_sl, gamma), q_c);
                    end
                end
            otherwise
                error(message('aero:correctairspeed:unknownAirspeed'));
        end
    otherwise
        error(message('aero:correctairspeed:unknownAirspeed'));
end
out = reshape(out,sizeIn);
end

%% Helper Functions
function out = asLookup(in, P0, a, Regime, subOut)
velTable = sprintf('aeroasV%s_m%d_table', subOut, Regime);
load('aeroascorrdata.mat', velTable, sprintf('aeroasV_m%d_bp', Regime), 'aeroasP_bp', 'aeroasSOS_bp');
[Vbp, Pbp, SOSbp]= meshgrid(eval(sprintf('aeroasV_m%d_bp', Regime)), aeroasP_bp, aeroasSOS_bp);
out = interp3( Vbp, Pbp, SOSbp, eval(velTable), in, P0, a, 'spline');
end

function q_c = asToQc(vel, P0, a, gamma)
if vel < a
    % Use Subsonic Equation
    q_c = P0 * (((1 + (gamma-1)/2 * (vel/a)^2)).^(gamma/(gamma-1)) - 1);
else
    % Use Supersonic Equation
    q_c = (1+gamma)/2 * vel^2 / a^2 * P0 * ...
        ((gamma+1)^2 / (4*gamma - 2*(gamma-1)*a^2/vel^2))^(1/(gamma-1)) - P0;
end
end

function vel = qcToAS(q_c, P0, a, gamma)
% Use Subsonic Equation
vel = a * sqrt(2/(gamma-1) * ((q_c/P0+1)^((gamma-1)/gamma) - 1));
if vel < a
    return;
else
    % Use Supersonic Equation 
    tolerance = 1e-5;
    maxCount = 30;
    i = 0;
    normVel = 1e2;
    while normVel > tolerance
        i = i + 1;
        if i > maxCount || ~isreal(vel)
            error(message('aero:correctairspeed:noConvergence'));
        end
        velPrev = vel;
        vel = a * sqrt(2/(1+gamma) * ((q_c/P0) + 1) / ...
            ((gamma+1)^2/(4*gamma - 2*(gamma-1)*a^2/vel^2))^(1/(gamma-1)));
        normVel = norm(velPrev - vel);
    end
end
end

% [EOF]