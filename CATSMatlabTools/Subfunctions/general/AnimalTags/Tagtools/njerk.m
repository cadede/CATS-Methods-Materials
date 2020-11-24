function    j = njerk(A,fs)

%     j = njerk(A)            % A is a sensor structure
%     or
%     j = njerk(A,fs)         % A is a matrix
%
%     Compute the norm-jerk from triaxial acceleration data. The norm-jerk
%		is ||dA/dt||, where ||x|| is the 2-norm of x, i.e., the square-root of the 
%		sum of the squares of each axis.
%
%		Inputs:
%     A is a sensor structure or a 3 column acceleration matrix with columns 
%      [ax ay az]. Acceleration can be in any consistent unit, e.g., g or m/s^2. 
%      A can be in any frame as the norm-jerk is rotation independent. A must 
%      have at least 2 rows (i.e., n>=2) and be regularly sampled.
%     fs is the sampling rate in Hz of the acceleration signals. This is used to
%		 estimate the differential by a first-order difference. fs is only required
%      if A is not a sensor structure.
%
%	   Result:
%		j is a column vector with the same number of rows as in A. If the unit of A 
%		is m/s^2, the norm-jerk has unit m/s^3. If the unit of A is g, the norm-jerk 
%	   has unit g/s. As j is the norm of the jerk, it is always positive or zero (if
%		the acceleration is constant). The final value in j is always 0 because the last
%		finite difference cannot be calculated.
%
%		Example:
%		njerk([1,2,3;2,2,4;1,-2,4;4,4,4],5)
%		result is: [7.0711;20.6155;33.541;0]
%
%    Valid: Matlab, Octave
%    markjohnson@st-andrews.ac.uk
%    Last modified: 15 aug 2012

if nargin<1,
   help njerk
   return
end

if isstruct(A),
   [A,fs]=sens2var(A,'regular') ;
   if isempty(A), return, end
elseif nargin<2,
   help njerk
   return
end

j = [fs*sqrt(sum(diff(A).^2,2));0] ;
