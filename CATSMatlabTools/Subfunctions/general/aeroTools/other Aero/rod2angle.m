function [r1,r2,r3] = rod2angle( R , varargin )
%  ROD2ANGLE Convert Rodrigues vector to rotation angles.
%   [R1, R2, R3] = ROD2ANGLE( ROD ) calculates the set of
%   rotation angles, R1, R2, R3, for a given Rodrigues vector, ROD.  
%
%   [R1, R2, R3] = ROD2ANGLE( ROD, S ) calculates the set of rotation
%   angles, R1, R2, R3, for a given Rodrigues vector, ROD, and a
%   specified rotation sequence, S. 
%
%   Inputs for rod2angle are:
%   ROD: M-by-3 matrix containing M Rodrigues vector. Each element
%   of ROD must be a real number.
%
%   S: Rotation sequence. The default rotation sequence is 'ZYX', where the
%   rotation angle order is:
%   R1  z-axis rotation
%   R2  y-axis rotation
%   R3  x-axis rotation 
%   
%   Possible sequences: 'ZYX', 'ZYZ', 'ZXY', 'ZXZ',
%   'YXZ', 'YXY', 'YZX', 'YZY', 'XYZ', 'XYX', 'XZY', and 'XZX'
%
%   Outputs for rod2angle are:
%   R1: M array of first rotation angles, in radians.
%
%   R2: M array of second rotation angles, in radians.
%
%   R3: M array of third rotation angles, in radians.
%
%
%   Examples:
%
%   Determine the rotation angles from Rodrigues vector:
%      r = [.1 .2 -.1]:
%      [yaw, pitch, roll] = rod2angle( r )
%
%   Determine the rotation angles from multiple Rodrigues vectors:
%      r = [.1 .2 -.1;.76 -.4 .36];
%      [pitch, roll, yaw] = rod2angle( r, 'YXZ' )
%
%   See also ANGLE2ROD, DCM2ROD, ROD2DCM, ROD2QUAT.

%   Copyright 2016-2018 The MathWorks, Inc.

%   Limitations: 
%   The limitations for the 'ZYX', 'ZXY', 'YXZ', 'YZX', 'XYZ', and 'XZY'
%   sequences generate an R2 angle that lies between +/- pi/2 radians (+/-
%   90 degrees), and R1 and R3 angles that lie between +/- pi radians (+/-
%   180 degrees).
%
%   The limitations for the 'ZYZ', 'ZXZ', 'YXY', 'YZY', 'XYX', and 'XZX'
%   sequences generate an R2 angle that lies between 0 and pi radians (180
%   degrees), and R1 and R3 angles that lie between +/- pi (+/- 180
%   degrees).
%   
%   Rodrigues transformation is not defined for rotation angles equal to
%   +/- pi radians (+/- 180 deg).
%
%   Reference:
%   [1] Dai, J.S. "Euler-Rodrigues formula variations, quaternion
%   conjugation and intrinsic connections." Mechanism and Machine Theory,
%   92, 144-152. Elsevier.

narginchk(1, 2);

%Validate inputs
validateattributes(R,{'single','double'},{'real','finite','nonempty','ncols',3});
% Set rotation order
if nargin == 1
    type = 'zyx';
else
    validateattributes(varargin{1},{'char','string'},{'nonempty'});
    type = lower(validatestring(varargin{1},{'ZYX', 'ZYZ', 'ZXY', 'ZXZ', ...
            'YXZ', 'YXY', 'YZX', 'YZY', 'XYZ', 'XYX', 'XZY','XZX'}));
end
len = size(R,1);
% Initalize outputs
r1 = zeros(len,1);
r2 = r1;
r3 = r1;
for k=1:len
    n = norm(R(k,:));
    if n ~= 0
        th = 2*atan(n);
        s = R(k,:)/n;
        cth = cos(th);
        sth = sin(th);
        sx = s(1);
        sy = s(2);
        sz = s(3);
        % Calculate rotation angles based on Rodrigues vector and rotation
        % order
        switch type
            case 'zyx'
                [r1(k),r2(k),r3(k)] = threeaxisrot(sx*sy*(1-cth)+sz*sth,...
                    sx^2+(1-sx^2)*cth,sy*sth-sx*sz*(1-cth),...
                    sy*sz*(1-cth)+sx*sth,sz^2+(1-sz^2)*cth);
            case 'zyz'
                [r1(k),r2(k),r3(k)] = twoaxisrot(sy*sz*(1-cth)-sx*sth,...
                    sx*sz*(1-cth)+sy*sth,sz^2+(1-sz^2)*cth,...
                    sy*sz*(1-cth)+sx*sth,sy*sth-sx*sz*(1-cth));

            case 'zxy'
                [r1(k),r2(k),r3(k)] = threeaxisrot(sz*sth-sx*sy*(1-cth),...
                    sy^2+(1-sy^2)*cth,sy*sz*(1-cth)+sx*sth,...
                    sy*sth-sx*sz*(1-cth),sz^2+(1-sz^2)*cth);

            case 'zxz'
                [r1(k),r2(k),r3(k)] = twoaxisrot(sx*sz*(1-cth)+sy*sth,...
                    sx*sth-sy*sz*(1-cth),sz^2+(1-sz^2)*cth,...
                    sx*sz*(1-cth)-sy*sth,sy*sz*(1-cth)+sx*sth);

            case 'yxz'
                [r1(k),r2(k),r3(k)] = threeaxisrot(sx*sz*(1-cth)+sy*sth,...
                    sz^2+(1-sz^2)*cth,sx*sth-sy*sz*(1-cth),...
                    sx*sy*(1-cth)+sz*sth,sy^2+(1-sy^2)*cth);

            case 'yxy'
                [r1(k),r2(k),r3(k)] = twoaxisrot(sx*sy*(1-cth)-sz*sth,...
                    sy*sz*(1-cth)+sx*sth,sy^2+(1-sy^2)*cth,...
                    sx*sy*(1-cth)+sz*sth,sx*sth-sy*sz*(1-cth));

            case 'yzx'       
                [r1(k),r2(k),r3(k)] = threeaxisrot(sy*sth-sx*sz*(1-cth),...
                    sx^2+(1-sx^2)*cth,sx*sy*(1-cth)+sz*sth,...
                    sx*sth-sy*sz*(1-cth),sy^2+(1-sy^2)*cth);

            case 'yzy'
                [r1(k),r2(k),r3(k)] = twoaxisrot(sy*sz*(1-cth)+sx*sth,...
                    sz*sth-sx*sy*(1-cth),sy^2+(1-sy^2)*cth,...
                    sy*sz*(1-cth)-sx*sth,sx*sy*(1-cth)+sz*sth);

            case 'xyz'
                [r1(k),r2(k),r3(k)] = threeaxisrot(sx*sth-sy*sz*(1-cth),...
                    sz^2+(1-sz^2)*cth,sx*sz*(1-cth)+sy*sth,...
                    sz*sth-sx*sy*(1-cth),sx^2+(1-sx^2)*cth);

            case 'xyx'
                [r1(k),r2(k),r3(k)] = twoaxisrot(sx*sy*(1-cth)+sz*sth,...
                    sy*sth-sx*sz*(1-cth),sx^2+(1-sx^2)*cth,...
                    sx*sy*(1-cth)-sz*sth,sx*sz*(1-cth)+sy*sth);

            case 'xzy'
                [r1(k),r2(k),r3(k)] = threeaxisrot(sy*sz*(1-cth)+sx*sth,...
                    sy^2+(1-sy^2)*cth,sz*sth-sx*sy*(1-cth),...
                    sx*sz*(1-cth)+sy*sth,sx^2+(1-sx^2)*cth);

            case 'xzx'
                [r1(k),r2(k),r3(k)] = twoaxisrot(sx*sz*(1-cth)-sy*sth,...
                    sx*sy*(1-cth)+sz*sth,sx^2+(1-sx^2)*cth,...
                    sx*sz*(1-cth)+sy*sth,sz*sth-sx*sy*(1-cth));
        end
    end
end
    function [r1,r2,r3] = threeaxisrot(r11, r12, r21, r31, r32)
        % find angles for rotations about X, Y, and Z axes
        r1 = atan2( r11, r12 );
        r21(r21 < -1) = -1;
        r21(r21 > 1) = 1;
        r2 = asin( r21 );
        r3 = atan2( r31, r32 );
    end

    function [r1,r2,r3] = twoaxisrot(r11, r12, r21, r31, r32)
        % Check for singularities
        if (r11==0 && r12==0)||(r31==0 && r32==0)
            [r1,r2,r3] = dcm2angle(rod2dcm(R(k,:)),type,'zeror3');
        else
            r1 = atan2( r11, r12 );
            r21(r21 < -1) = -1;
            r21(r21 > 1) = 1;
            r2 = acos( r21 );
            r3 = atan2( r31, r32 );
        end
    end
end
