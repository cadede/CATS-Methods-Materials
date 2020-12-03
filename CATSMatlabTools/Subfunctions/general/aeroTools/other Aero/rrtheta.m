function th = rrtheta( temp, mach, gamma )
%  RRTHETA Compute relative temperature ratio.
%   TH  = RRTHETA( T0, MACH, G ) computes M temperature relative ratios,
%   TH, from M static temperatures, T0, M Mach numbers, MACH and M specific
%   heat ratios, G.  T0 must be in kelvin.
%
%   Limitations: 
%
%   For cases in which total temperature ratio is desired (Mach number is
%   nonzero), the total temperature is calculated assuming perfect gas
%   (with constant molecular weight, constant pressure specific heat, and
%   constant specific heat ratio) and dry air. 
%
%   Examples:
%
%   Determine relative temperature ratio for three temperatures:
%      theta = rrtheta([ 273.15  310.9278  373.15 ], 0.5, 1.4 )
%
%   Determine relative temperature ratio for three temperatures at
%   different conditions: 
%      theta = rrtheta([ 273.15  310.9278  373.15 ], 0.5, [1.4 1.35 1.4])
%
%   Determine relative temperature ratio for three temperature at three
%   different conditions:
%      theta = rrtheta([ 273.15  310.9278  373.15 ], [ 0.5 1 2], [1.4 1.35 1.4])
%
%   See also RRDELTA, RRSIGMA.
%

%   Copyright 2000-2018 The MathWorks, Inc.
 
%   Reference: 
%   [1] Aeronautical Vestpocket Handbook, United Technologies Pratt &
%   Whitney, August, 1986.

if ~isnumeric( temp )
    error(message('aero:rrtheta:notNumeric1'));
end

if ~isnumeric( mach )
    error(message('aero:rrtheta:notNumeric2'));
end

if ~isnumeric( gamma )
    error(message('aero:rrtheta:notNumeric3'));
end

try
    th = ( 1 + 0.5.*(gamma - 1).*mach.^2 ).*( temp./288.15 );
catch wrongdim
    newExc = MException('aero:rrtheta:wrongDimension', ...
        getString(message('aero:rrtheta:wrongDimension')));
    newExc = newExc.addCause(wrongdim);
    throw(newExc);
end
