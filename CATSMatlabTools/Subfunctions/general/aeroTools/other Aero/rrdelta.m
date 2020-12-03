function delta = rrdelta( pres, mach, gamma )
%  RRDELTA Compute relative pressure ratio.
%   D  = RRDELTA( P0, MACH, G ) computes M pressure relative ratios,
%   D, from M static pressures, P0, M Mach numbers, MACH and M specific
%   heat ratios, G.  P0 must be in pascals.
%
%   Limitations: 
%
%   For cases in which total pressure ratio is desired (Mach number is
%   nonzero), the total pressures are calculated assuming perfect gas (with
%   constant molecular weight, constant pressure specific heat, and
%   constant specific heat ratio) and dry air.  
% 
%   Examples:
%
%   Determine relative pressure ratio for three pressures:
%      delta = rrdelta( [ 101325  22632.0672  4328.1393 ], 0.5, 1.4 )
%
%   Determine relative pressure ratio for three pressures at
%   different conditions: 
%      delta = rrdelta( [ 101325  22632.0672  4328.1393 ], 0.5, [1.4 1.35 1.4])
%
%   Determine relative pressure ratio for three pressures at three
%   different conditions:
%      delta = rrdelta( [ 101325  22632.0672  4328.1393 ], [ 0.5 1 2], [1.4 1.35 1.4])
%
%   See also RRSIGMA, RRTHETA.

%   Copyright 2000-2018 The MathWorks, Inc.

%   Reference: 
%   [1] Aeronautical Vestpocket Handbook, United Technologies Pratt &
%   Whitney, August, 1986.

if ~isnumeric( pres )
    error(message('aero:rrdelta:notNumeric1'));
end

if ~isnumeric( mach )
    error(message('aero:rrdelta:notNumeric2'));
end

if ~isnumeric( gamma )
    error(message('aero:rrdelta:notNumeric3'));
end

try
    delta = ( pres./101325 ).*( 1 + 0.5.*(gamma - 1).*mach.^2 ).^(gamma./(gamma - 1));
catch wrongdim
    newExc = MException('aero:rrdelta:wrongDimension', ...
        getString(message('aero:rrdelta:wrongDimension')));
    newExc = newExc.addCause(wrongdim);
    throw(newExc);
end
