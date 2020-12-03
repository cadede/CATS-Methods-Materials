function sigma = rrsigma( rho, mach, gamma )
%  RRSIGMA Compute relative density ratio.
%   S  = RRSIGMA( RHO, MACH, G ) computes M density relative ratios,
%   S, from M static densities, RHO, M Mach numbers, MACH and M specific
%   heat ratios, G.  RHO must be in kilograms per meter cubed.
%
%   Limitations: 
%
%   For cases in which total density ratio is desired (Mach number is
%   nonzero), the total density is calculated assuming perfect gas (with
%   constant molecular weight, constant pressure specific heat, and
%   constant specific heat ratio) and dry air.  
% 
%   Examples:
%
%   Determine relative density ratio for three densities:
%      sigma = rrsigma( [ 1.225  0.3639  0.0953 ], 0.5, 1.4 )
%
%   Determine relative density ratio for three densities at
%   different conditions: 
%      sigma = rrsigma( [ 1.225  0.3639  0.0953 ], 0.5, [1.4 1.35 1.4])
%
%   Determine relative density ratio for three densities at three
%   different conditions:
%      sigma = rrsigma( [ 1.225  0.3639  0.0953 ], [ 0.5 1 2], [1.4 1.35 1.4])
%
%   See also RRDELTA, RRTHETA.

%   Copyright 2000-2018 The MathWorks, Inc.

%   Reference: 
%   [1] Aeronautical Vestpocket Handbook, United Technologies Pratt &
%   Whitney, August, 1986.

if ~isnumeric( rho )
    error(message('aero:rrsigma:notNumeric1'));
end

if ~isnumeric( mach )
    error(message('aero:rrsigma:notNumeric2'));
end

if ~isnumeric( gamma )
    error(message('aero:rrsigma:notNumeric3'));
end
    
try
    sigma = ( rho./1.225 ).*( 1 + 0.5.*(gamma - 1).*mach.^2 ).^(1./(gamma - 1));
catch wrongdim
    newExc = MException('aero:rrsigma:wrongDimension', ...
        getString(message('aero:rrsigma:wrongDimension')));
    newExc = newExc.addCause(wrongdim);
    throw(newExc);
end
