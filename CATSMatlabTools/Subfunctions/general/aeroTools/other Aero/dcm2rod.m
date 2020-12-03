function R = dcm2rod( dcm, varargin)
% DCM2ROD Convert direction cosine matrix to Rodrigues vector.
%   R = DCM2ROD( N ) calculates the Rodrigues vector, R, for a given
%   direction cosine matrix, N. N is a 3-by-3-by-M matrix containing M
%   orthogonal direction cosine matrices.  R returns an M-by-3 matrix
%   containing M Rodrigues vectors.
%
%   R = DCM2ROD( N, ACTION ) uses ACTION for error handling during
%   rotation matrix validation. Specify if an invalid rotation matrix
%   invokes a 'Warning', 'Error', or no action ('None'). The default is
%   'None'.
%
%   R = DCM2ROD( N, ACTION, TOL ) uses TOL as a relative error tolerance
%   for rotation matrix validation. It is a scalar value. Rotation matrix
%   validation confirms that the matrix is orthogonal [transpose(N) * N
%   == 1 +/- TOL] and proper [det(N) == 1 +/- TOL]. The default value is
%   eps(2).
%
%   Example:
%
%   Determine the Rodrigues vector from direction cosine matrix:
%      DCM = [0.433 0.75 0.5;-0.25 -0.433 0.866;0.866 -0.5 0.0];
%      r = dcm2rod( DCM )
%
%   See also ROD2DCM, ROD2ANGLE, ROD2QUAT, QUAT2ROD.

%   Copyright 2016-2017 The MathWorks, Inc.

%   Limitations:
%   Rodrigues transformation is not defined for rotation angles equal to
%   +/- pi radians (+/- 180 deg).
%
%   Reference:
%   [1] Dai, J.S. "Euler-Rodrigues formula variations, quaternion
%   conjugation and intrinsic connections." Mechanism and Machine Theory,
%   92, 144-152. Elsevier.

narginchk(1, 3);

%Validate inputs
if nargin < 2
    validateDCM(dcm);
elseif nargin < 3
    validateDCM(dcm, varargin{1});
else
    validateDCM(dcm, varargin{1}, varargin{2});
end

len = size(dcm,3);
% Initialize output
R = zeros(len,3);
for k=1:len
    % Calculate angle rotation
    th = acos((trace(dcm(:,:,k))-1)/2);
    if th ~= 0
        % Calculate unit vector rotation
        sy = (dcm(3,1,k)-dcm(1,3,k))/(2*sin(th));
        sx = (dcm(2,3,k)-dcm(3,2,k))/(2*sin(th));
        sz = (dcm(1,2,k)-dcm(2,1,k))/(2*sin(th));
        % Calculate Rodrigues vector
        R(k,:) = tan(th/2)*[sx sy sz];
    end
end