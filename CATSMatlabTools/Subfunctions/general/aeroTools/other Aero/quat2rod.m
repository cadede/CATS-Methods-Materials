function R = quat2rod(qin)
% QUAT2ROD Convert quaternion to Rodrigues vector
%   R = QUAT2ROD( Q ) calculates the Rodrigues vector, R, for given
%   quaternion Q. Q is an M array of quaternions and R returns an M-by-3 
%   matrix containing M Rodrigues vectors. Q has its scalar number as the
%   first column.
%
%   Example:
%
%   Determine the Rodrigues vector from quaternion:
%      q = [-0.7071 0 0.7071 0];
%      r = quat2rod( q )
%
%   See also ROD2DCM, ROD2ANGLE, ROD2QUAT.

% Copyright 2016 The MathWorks, Inc.

%   Limitations:
%   Rodrigues transformation is not defined for rotation angles equal to
%   +/- 180 deg or +/- pi radians.
%
%   Reference:
%   [1] Dai, J.S. Euler-Rodrigues formula variations, quaternion
%   conjugation and intrinsic connections. Mechanism and Machine Theory, 92
%   144-152. Elsevier.

% Validate input
validateattributes(qin,{'single','double'},{'real','nonempty','ncols',4});
% Normalize quaternions
qinn = quatnormalize(qin);
len = size(qinn,1);
%Initialie output
R = zeros(len,3);
for k=1:len
    % Calculate rotation angle
    thalf = acos(qinn(k,1));
    if thalf ~= 0
        % Calculate Rodrigues vector
        R(k,:) = qinn(k,2:4)/cos(thalf);
    end
end