function [r1,r2,r3] = quat2angle( q, varargin )
%  QUAT2ANGLE Convert quaternion to rotation angles.
%   [R1, R2, R3] = QUAT2ANGLE( Q ) calculates the calculates the set of
%   rotation angles, R1, R2, R3, for a given quaternion, Q.  Input Q is an
%   M-by-4 matrix containing M quaternions.  R1 returns an M array of
%   first rotation angles.  R2 returns an M array of second rotation
%   angles.  R3 returns an M array of third rotation angles. Each element
%   of Q must be a real number.  Additionally, Q has its scalar number as 
%   the first column. Rotation angles are output in radians.   
%
%   [R1, R2, R3] = QUAT2ANGLE( Q, S ) calculates the set of rotation
%   angles, R1, R2, R3, for a given quaternion, Q, and a
%   specified rotation sequence, S. 
%
%   The default rotation sequence is 'ZYX' where the order of rotation
%   angles for the default rotation are R1 = Z Axis Rotation, R2 = Y Axis
%   Rotation, and R3 = X Axis Rotation. 
%
%   All rotation sequences, S, are supported: 'ZYX', 'ZYZ', 'ZXY', 'ZXZ',
%   'YXZ', 'YXY', 'YZX', 'YZY', 'XYZ', 'XYX', 'XZY', and 'XZX'.
%
%   Examples:
%
%   Determine the rotation angles from q = [1 0 1 0]:
%      [yaw, pitch, roll] = quat2angle( [1 0 1 0] )
%
%   Determine the rotation angles from multiple quaternions:
%      q = [1 0 1 0; 1 0.5 0.3 0.1];
%      [pitch, roll, yaw] = quat2angle( q, 'YXZ' )
%
%   See also ANGLE2DCM, DCM2ANGLE, DCM2QUAT, ANGLE2QUAT, QUAT2DCM. 

%   Copyright 2000-2018 The MathWorks, Inc.

%   Limitations: 
%   The limitations for the 'ZYX', 'ZXY', 'YXZ', 'YZX', 'XYZ', and 'XZY'
%   implementations generate an R2 angle that lies between +/- 90 
%   degrees, and R1 and R3 angles that lie between +/- 180 degrees. 
%
%   The limitations for the 'ZYZ', 'ZXZ', 'YXY', 'YZY', 'XYX', and 'XZX'
%   implementations generate an R2 angle that lies between 0 and 
%   180 degrees, and R1 and R3 angles that lie between +/- 180 degrees. 

narginchk(1, 2);

if nargin == 1
    type = 'zyx';
else
    if ischar( varargin{1} ) || isstring( varargin{1} )
        type = char(varargin{1});
    else
        error(message('aero:quat2angle:notChar'));
    end
end

qin = quatnormalize( q );

switch lower( type )
    case 'zyx'
        [r1,r2,r3] = threeaxisrot( 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), ...
                                   qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2, ...
                                  -2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)), ...
                                   2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)), ...
                                   qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2);
    case 'zyz'
        [r1,r2,r3] = twoaxisrot( 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), ...
                                 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), ...
                                 qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2, ...
                                 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)), ...
                                -2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)));
                
    case 'zxy'
        [r1,r2,r3] = threeaxisrot( -2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), ...
                                    qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, ...
                                    2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)), ...
                                   -2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)), ...
                                    qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2);

    case 'zxz'
        [r1,r2,r3] = twoaxisrot( 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), ...
                                -2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), ...
                                 qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2, ...
                                 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)), ...
                                 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)));

    case 'yxz'
        [r1,r2,r3] = threeaxisrot( 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), ...
                                   qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2, ...
                                  -2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), ...
                                   2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), ...
                                   qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2);
       
    case 'yxy'
        [r1,r2,r3] = twoaxisrot( 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), ...
                                 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)), ...
                                 qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, ...
                                 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), ...
                                -2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)));
      
    case 'yzx'       
        [r1,r2,r3] = threeaxisrot( -2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)), ...
                                    qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2, ...
                                    2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), ...
                                   -2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), ...
                                    qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2);
        
    case 'yzy'
        [r1,r2,r3] = twoaxisrot( 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)), ...
                                -2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), ...
                                 qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, ...
                                 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), ...
                                 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)));

    case 'xyz'
        [r1,r2,r3] = threeaxisrot( -2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), ...
                                    qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2, ...
                                    2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), ...
                                   -2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), ...
                                    qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2);
        
    case 'xyx'
        [r1,r2,r3] = twoaxisrot( 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), ...
                                -2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)), ...
                                 qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2, ...
                                 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), ...
                                 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)));
        
    case 'xzy'
        [r1,r2,r3] = threeaxisrot( 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)), ...
                                   qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, ...
                                  -2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), ...
                                   2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), ...
                                   qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2);
        
    case 'xzx'
        [r1,r2,r3] = twoaxisrot( 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)), ...
                                 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), ...
                                 qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2, ...
                                 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), ...
                                -2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)));
    otherwise
        error(message('aero:quat2angle:unknownRotation', type));
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
        r1comp = and((r11==0),(r12==0));
        r3comp = and((r31==0),(r32==0));
        OrR1R3 = or(r1comp,r3comp);
        NorR1R3 = not(OrR1R3);
        if any(OrR1R3)
            [r1(OrR1R3),r2(OrR1R3),r3(OrR1R3)] = dcm2angle(quat2dcm(qin(OrR1R3,:)),type,'zeror3');      
        end
        if any(NorR1R3)
            r1(NorR1R3) = atan2( r11(NorR1R3), r12(NorR1R3) );
            r21(r21 < -1) = -1;
            r21(r21 > 1) = 1;
            r2(NorR1R3) = acos( r21(NorR1R3) );
            r3(NorR1R3) = atan2( r31(NorR1R3), r32(NorR1R3) );
        end
        r1 = r1';
        r2 = r2';
        r3 = r3';
    end
end
