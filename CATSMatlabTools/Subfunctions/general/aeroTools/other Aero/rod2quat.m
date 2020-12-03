function qout = rod2quat(R)
% ROD2QUAT Convert Rodrigues vector to quaternion.
%   Q = ROD2QUAT( R ) calculates the quaternion, Q, for given Rodrigues
%   vector R.
%
%   Input for rod2quat is:
%   R:  M array of Rodrigues vectors.
%
%   Output for rod2quat is:
%   Q: M-by-4 matrix of M quaternions. Q has its scalar number as the
%   first column.
%
%
%   Examples:
%   Determine the quaternion from Rodrigues vector:
%      r = [.1 .2 -.1]
%      q = rod2quat(r)
%
%   See also QUAT2ROD, DCM2ROD, DCM2QUAT, ANGLE2QUAT, QUAT2DCM. 

% Copyright 2016 The MathWorks, Inc.

%   Limitations:
%   Rodrigues transformation is not defined for rotation angles equal to
%   +/- pi radians (+/- 180 deg).
%   
%   Reference:
%   [1] Dai, J.S. "Euler-Rodrigues formula variations, quaternion
%   conjugation and intrinsic connections." Mechanism and Machine Theory,
%   92, 144-152. Elsevier.

%Validate inputs
validateattributes(R,{'single','double'},{'real','nonempty','ncols',3});
len = size(R,1);
% Initialize outputs
qout = zeros(len,4);
for k = 1:len
    n = norm(R(k,:));
    if n~=0
        % Calculate rotation angle
        thalf = atan(n);
        % Calculate equivalent quaternion
        qout(k,1)= cos(thalf);
        qout(k,2:4) = R(k,:)*sin(thalf)/n;
    end
end