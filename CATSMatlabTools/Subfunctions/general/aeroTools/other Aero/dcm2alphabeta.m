function [ alpha, beta ] = dcm2alphabeta( dcm, varargin ) 
%  DCM2ALPHABETA Convert direction cosine matrix to angle of attack and sideslip angle.
%   [A, B] = DCM2ALPHABETA( N ) calculates the angle of attack and
%   sideslip angle, A, B, for given direction cosine matrix, N. N is a
%   3-by-3-by-M matrix containing M orthogonal direction cosine matrices. N
%   performs the coordinate transformation of a vector in body-axes into a
%   vector in wind-axes.A is an M array of angles of attack.  B is an M
%   array of sideslip angles. Angles of attack and sideslip angles are
%   output in radians.
%
%   [A, B] = DCM2ALPHABETA( N, ACTION ) uses ACTION for error handling
%   during rotation matrix validation. Specify if an invalid rotation
%   matrix invokes a 'Warning', 'Error', or no action ('None'). The default
%   is 'None'.
%
%   [A, B] = DCM2ALPHABETA( N, ACTION, TOL ) uses TOL as a relative error
%   tolerance for rotation matrix validation. It is a scalar value.
%   Rotation matrix validation confirms that the matrix is orthogonal
%   [transpose(N) * N == 1 +/- TOL] and proper [det(N) == 1 +/- TOL]. The
%   default value is eps(2).
%
%   Examples:
%
%   Determine the angle of attack and sideslip angle from direction cosine matrix:
%      dcm = [ 0.8926    0.1736    0.4162; ...
%             -0.1574    0.9848   -0.0734; ...
%             -0.4226         0    0.9063]; 
%      [alpha, beta] = dcm2alphabeta( dcm )
%
%   Determine the angle of attack and sideslip angle from multiple direction cosine matrices:
%      dcm =        [ 0.8926    0.1736    0.4162; ...
%                    -0.1574    0.9848   -0.0734; ...
%                    -0.4226         0    0.9063]; 
%      dcm(:,:,2) = [ 0.9811    0.0872    0.1730; ...
%                    -0.0859    0.9962   -0.0151; ...
%                    -0.1736         0    0.9848]; 
%      [alpha, beta] = dcm2alphabeta( dcm )
%
%   See also ANGLE2DCM, DCM2ANGLE, DCMBODY2WIND.


%   Copyright 2000-2018 The MathWorks, Inc.

narginchk(1, 3);

%Validate inputs
if nargin < 2
    validateDCM(dcm);
elseif nargin < 3
    validateDCM(dcm, varargin{1});
else
    validateDCM(dcm, varargin{1}, varargin{2});
end

beta = asin(dcm(1,2,:));
alpha = asin(-dcm(3,1,:));

beta = beta(:);
alpha = alpha(:);
