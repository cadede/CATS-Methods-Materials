function ROD = angle2rod( r1, r2, r3, varargin )
%  ANGLE2ROD Convert rotation angles to Rodrigues vector.
%   ROD = ANGLE2ROD( R1, R2, R3 ) calculates the Rodrigues vector, ROD,
%   for a given set of rotation angles R1, R2, R3.
%
%   ROD = ANGLE2ROD( R1, R2, R3, S ) calculates the Rodrigues vector, ROD,
%   for a given set of rotation angles, R1, R2, R3, and a specified
%   rotation sequence, S.  
%
%   Inputs for angle2rod are:
%   R1: M array of first rotation angles, in radians.
%
%   R2: M array of second rotation angles, in radians.
%
%   R3: M array of third rotation angles, in radians.
%
%   S: Rotation sequence. The default rotation sequence is 'ZYX', where the
%   rotation angle order is:
%   R1  z-axis rotation
%   R2  y-axis rotation
%   R3  x-axis rotation 
%   
%   Possible sequences: 'ZYX', 'ZYZ', 'ZXY', 'ZXZ',
%   'YXZ', 'YXY', 'YZX', 'YZY', 'XYZ', 'XYX', 'XZY', and 'XZX'.
%
%   Output for angle2rod is:
%   ROD: M-by-3 matrix containing M Rodrigues vector. 
%
%   Examples:
%
%   Determine the Rodrigues vector from rotation angles:
%      yaw = 0.7854; 
%      pitch = 0.1; 
%      roll = 0;
%      r = angle2rod( yaw, pitch, roll )
%
%   Determine the Rodrigues vectors from multiple rotation angles:
%      yaw = [0.7854 0.5]; 
%      pitch = [0.1 0.3]; 
%      roll = [0 0.1];
%      r = angle2rod( pitch, roll, yaw, 'YXZ' )
%
%   See also DCM2ROD, ROD2DCM, ROD2ANGLE.

%   Copyright 2016 The MathWorks, Inc.

%   Limitations:
%   Rodrigues transformation is not defined for rotation angles equal to
%   pi radians (+/- 180 deg).
%
%   Reference:
%   [1] Dai, J.S. "Euler-Rodrigues formula variations, quaternion
%   conjugation and intrinsic connections." Mechanism and Machine Theory,
%   92, 144-152. Elsevier.

narginchk(3,4);
len = length(r1);
validateattributes(r1,{'single','double'},{'real','nonempty','vector'});
validateattributes(r2,{'single','double'},{'real','nonempty','vector','numel',len});
validateattributes(r3,{'single','double'},{'real','nonempty','vector','numel',len});

if nargin == 3
    type = 'zyx';
else
    validateattributes(varargin{1},{'char','string'},{'nonempty'});
    type = lower(validatestring(varargin{1},{'ZYX', 'ZYZ', 'ZXY', 'ZXZ','YXZ', ...
        'YXY', 'YZX', 'YZY', 'XYZ', 'XYX', 'XZY','XZX'}));
end

len = size(r1,1);
% Initialize outputs
ROD = zeros(len,3);
% Calculate sines and cosines 
cps = cos(r1(:));
sps = sin(r1(:));
cth = cos(r2(:));
sth = sin(r2(:));
cph = cos(r3(:));
sph = sin(r3(:));
% Calculate Rodrigues vector depending on rotation sequence
switch type
    case 'zyx'
        tr = cth.*cps + sph.*sth.*sps+cph.*cps+cth.*cph;
        r31 = cph.*sth.*cps+sph.*sps+sth;
        r23 = cth.*sph-cph.*sth.*sps+sph.*cps; 
        r12 = cth.*sps-sph.*sth.*cps+cph.*sps;
    case 'zyz'
        tr = cph.*cth.*cps-sph.*sps-sph.*cth.*sps+cph.*cps+cth;
        r31 = sth.*cps+cph.*sth;
        r23 = sph.*sth-sth.*sps;
        r12 = cph.*cth.*sps+sph.*cps+sph.*cth.*cps+cph.*sps;
    case 'zxy'
        tr = cph.*cps-sph.*sth.*sps+cth.*cps+cph.*cth;
        r31 = sph.*cps+cph.*sth.*sps+cth.*sph;
        r23 = sth-sph.*sps+cph.*sth.*cps;
        r12 = cph.*sps+sph.*sth.*cps+cth.*sps;
    case 'zxz'
        tr = cph.*cps-sph.*cth.*sps-sph.*sps+cph.*cth.*cps+cth;
        r31 = sth.*sps-sph.*sth;
        r23 = cph.*sth+sth.*cps;
        r12 = cph.*sps+sph.*cth.*cps+sph.*cps+cph.*cth.*sps;
    case 'yxz'
        tr = cph.*cps+sph.*sth.*sps+cth.*cph+cth.*cps;
        r31 = cth.*sps+cph.*sps-sph.*sth.*cps;
        r23 = sph.*sps+cph.*sth.*cps+sth;
        r12 = sph.*cth+sph.*cps-cph.*sth.*sps;
    case 'yxy'
        tr = cph.*cps-sph.*cth.*sps+cth-sph.*sps+cph.*cth.*cps;
        r31 = sph.*cps+cph.*cth.*sps+cph.*sps+sph.*cth.*cps;
        r23 = sth.*cps+sth.*cph;
        r12 = sth.*sph-sth.*sps;
    case 'yzx'
        tr = cth.*cps+cph.*cth-sph.*sth.*sps+cph.*cps;
        r31 = sph.*sth.*cps+cph.*sps+cth.*sps;
        r23 = cph.*sth.*sps+sph.*cps+sph.*cth;
        r12 = sth+cph.*sth.*cps-sph.*sps;
    case 'yzy'
        tr = cph.*cth.*cps-sph.*sps+cth-sph.*cth.*sps+cph.*cps;
        r31 = sph.*cth.*cps+cph.*sps+cph.*cth.*sps+sph.*cps;
        r23 = sth.*sps-sph.*sth;
        r12 = cph.*sth+sth.*cps;
    case 'xyz'
        tr = cph.*cth+cph.*cps-sph.*sth.*sps+cps.*cth;
        r31 = sth-sph.*sps+cph.*sth.*cps;
        r23 = cph.*sps+sph.*sth.*cps+sps.*cth;
        r12 = sph.*cps+cph.*sth.*sps+sph.*cth;
    case 'xyx'
        tr = cth+cph.*cps-sph.*cth.*sps-sph.*sps+cph.*cth.*cps;
        r31 = cph.*sth+sth.*cps;
        r23 = cph.*sps+sph.*cth.*cps+sph.*cps+cph.*cth.*sps;
        r12 = sth.*sps-sph.*sth;
    case 'xzy'
        tr = cph.*cth+cth.*cps+sph.*sth.*sps+cph.*cps;
        r31 = sph.*cth-cph.*sth.*sps+sph.*cps;
        r23 = cth.*sps-sph.*sth.*cps+cph.*sps;
        r12 = cph.*sth.*cps+sph.*sps+sth;
    case 'xzx'
        tr = cth+cph.*cth.*cps-sph.*sps-sph.*cth.*sps+cph.*cps;
        r31 = sph.*sth-sps.*sth;
        r23 = cph.*cth.*sps+sph.*cps+sph.*cth.*cps+cph.*sps;
        r12 = cps.*sth+cph.*sth;
end
% Calculate rotation angle
th = acos((tr-1)/2);
t0 = th ~= 0;
% Calculate unit rotation vector
sy(t0) = r31(t0)./(2*sin(th(t0)));
sx(t0) = r23(t0)./(2*sin(th(t0)));
sz(t0) = r12(t0)./(2*sin(th(t0)));
% Assemble Rodrigues vector
ROD(t0,:) = tan(th(t0)/2).*[sx(t0)' sy(t0)' sz(t0)'];

