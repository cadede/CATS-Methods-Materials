function DCM = rod2dcm(R)
% ROD2DCM Convert Rodrigues vector to direction cosine matrix.
%   DCM = ROD2DCM( R ) calculates the direction cosine matrix, for given
%   Rodrigues vector R.
%
%   Input for rod2dcm is:
%   R: M-by-3 matrix containing M
%   Rodrigues vectors.
%
%   Output for rod2dcm is:
%   DCM: 3-by-3-by-M matrix containing M
%   direction cosine matrices.
%
%
%   Examples:
%   Determine the direction cosine matrix from Rodrigues vector:
%      r = [.1 .2 -.1];
%      DCM = rod2dcm(r)
%
%   See also DCM2ROD, ANGLE2ROD, QUAT2ROD, ROD2ANGLE.

% Copyright 2016 The MathWorks, Inc.

%   Limitations:
%   Rodrigues transformation is not defined for rotation angles equal to
%   +/- pi radians (+/- 180 deg).
%
%   Reference:
%   [1] Dai, J.S. "Euler-Rodrigues formula variations, quaternion
%   conjugation and intrinsic connections." Mechanism and Machine Theory,
%   92, 144-152. Elsevier.

% Validate inputs
validateattributes(R,{'single','double'},{'real','ncols',3});
len = size(R,1);
% Initalize outputs
DCM = zeros(3,3,len);
for k = 1:len
    % Calculte norm
    n = norm(R(k,:));
    if n == 0
        DCM(:,:,k) = eye(3);
    else
        th = 2*atan(n);
        s = R(k,:)/n;
        cth = cos(th);
        sth = sin(th);
        DCM(:,:,k) = [  s(1)^2+(1-s(1)^2)*cth        s(1)*s(2)*(1-cth)+s(3)*sth  s(1)*s(3)*(1-cth)-s(2)*sth;...
                        s(1)*s(2)*(1-cth)-s(3)*sth   s(2)^2+(1-s(2)^2)*cth       s(2)*s(3)*(1-cth)+s(1)*sth;...
                        s(1)*s(3)*(1-cth)+s(2)*sth   s(2)*s(3)*(1-cth)-s(1)*sth  s(3)^2+(1-s(3)^2)*cth];
    end
end