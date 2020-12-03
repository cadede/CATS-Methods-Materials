function [r1,r2,r3] = dcm2angle( dcm, varargin )
%  DCM2ANGLE Create rotation angles from direction cosine matrix.
%   [R1, R2, R3]  = DCM2ANGLE( N ) calculates the set of rotation angles,
%   R1, R2, R3, for a given direction cosine matrix, N. N is a 3-by-3-by-M
%   matrix containing M orthogonal direction cosine matrices. R1 returns 
%   an M array of  first rotation angles.  R2 returns an M array of second
%   rotation angles.  R3 returns an M array of third rotation angles.
%   Rotation angles are output in radians.  
%
%   [R1, R2, R3] = DCM2ANGLE( N, S ) calculates the set of rotation
%   angles, R1, R2, R3, for a given direction cosine matrix, N, and a
%   specified rotation sequence, S. 
%
%   The default rotation sequence is 'ZYX' where the order of rotation
%   angles for the default rotation are R1 = Z Axis Rotation, R2 = Y Axis
%   Rotation, and R3 = X Axis Rotation. 
%
%   All rotation sequences, S, are supported: 'ZYX', 'ZYZ', 'ZXY', 'ZXZ',
%   'YXZ', 'YXY', 'YZX', 'YZY', 'XYZ', 'XYX', 'XZY', and 'XZX'.
%
%   [R1, R2, R3] = DCM2ANGLE( N, S, LIM ) uses LIM as an angle constraint.
%   LIM is a string specifying either 'Default' or 'ZeroR3'. See the
%   Limitations section for full definitions of angle constraints.
%
%   [R1, R2, R3] = DCM2ANGLE( N, S, LIM, ACTION ) uses ACTION for error
%   handling during rotation matrix validation. Specify if an invalid
%   rotation matrix invokes a 'Warning', 'Error', or no action ('None').
%   The default is 'None'.
%
%   [R1, R2, R3] = DCM2ANGLE( N, S, LIM, ACTION, TOL ) uses TOL as a relative
%   tolerance for rotation matrix validation. TOL is a scalar value.
%   Rotation matrix validation confirms that the matrix is orthogonal
%   [transpose(N) * N == 1 +/- TOL] and proper [det(N) == 1 +/- TOL]. The
%   default value is eps(2).
%
%   Limitations:
%
%   In general, the R1 and R3 angles lie between +/- 180 degrees and R2
%   lies either between +/- 90 degrees or between 0 and 180 degrees.  For
%   more information see the documentation. 
%
%   Examples:
%
%   Determine the rotation angles from direction cosine matrix:
%      dcm = [1 0 0; 0 1 0; 0 0 1];
%      [yaw, pitch, roll] = dcm2angle( dcm )
%
%   Determine the rotation angles from multiple direction cosine matrices:
%      dcm        = [ 1 0 0; 0 1 0; 0 0 1]; 
%      dcm(:,:,2) = [ 0.85253103550038   0.47703040785184  -0.21361840626067; ...
%                    -0.43212157513194   0.87319830445628   0.22537893734811; ...
%                     0.29404383655186  -0.09983341664683   0.95056378592206];
%      [pitch, roll, yaw] = dcm2angle( dcm, 'YXZ' )
%
%   See also ANGLE2DCM, DCM2QUAT, QUAT2DCM, QUAT2ANGLE.

%   Copyright 2000-2018 The MathWorks, Inc.

%   Limitations: 
%   The 'Default' limitations for the 'ZYX', 'ZXY', 'YXZ', 'YZX', 'XYZ',
%   and 'XZY' implementations generate an R2 angle that lies between +/- 90
%   degrees, and R1 and R3 angles that lie between +/- 180 degrees. 
%
%   The 'Default' limitations for the 'ZYZ', 'ZXZ', 'YXY', 'YZY', 'XYX',
%   and 'XZX' implementations generate an R2 angle that lies between 0 and
%   180 degrees, and R1 and R3 angles that lie between +/- 180 degrees. 
%
%   The 'ZeroR3' limitations for the 'ZYX', 'ZXY', 'YXZ', 'YZX', 'XYZ',
%   and 'XZY' implementations generate an R2 angle that lies between +/- 90
%   degrees, and R1 and R3 angles that lie between +/- 180 degrees.
%   However, when R2 is +/- 90 degrees, R3 is set to 0 degrees.
%
%   The 'ZeroR3' limitations for the 'ZYZ', 'ZXZ', 'YXY', 'YZY', 'XYX',
%   and 'XZX' implementations generate an R2 angle that lies between 0 and 180
%   degrees, and R1 and R3 angles that lie between +/- 180 degrees.
%   However, when R2 is 0 or +/- 180 degrees, R3 is set to 0 degrees.

narginchk(1, 5);

lim = 'default';

%Validate inputs
if nargin == 1
    type = 'zyx';
else
    if ischar( varargin{1} ) || isstring( varargin{1} )
        type = char(varargin{1});
    else
        error(message('aero:dcm2angle:notChar'));
    end
    if nargin == 3
        if ischar( varargin{2} ) || isstring( varargin{2} )
            lim = char(varargin{2});
        else
            error(message('aero:dcm2angle:notChar1'));
        end
    end
end
if nargin < 4
    validateDCM(dcm);
elseif nargin < 5
    validateDCM(dcm, varargin{3});
else
    validateDCM(dcm, varargin{3}, varargin{4});
end

switch lower( type )
    case 'zyx'
        %     [          cy*cz,          cy*sz,            -sy]
        %     [ sy*sx*cz-sz*cx, sy*sx*sz+cz*cx,          cy*sx]
        %     [ sy*cx*cz+sz*sx, sy*cx*sz-cz*sx,          cy*cx]

        [r1,r2,r3] = threeaxisrot( dcm(1,2,:), dcm(1,1,:), -dcm(1,3,:), ...
                                   dcm(2,3,:), dcm(3,3,:), ...
                                  -dcm(2,1,:), dcm(2,2,:));

    case 'zyz'
        %     [  cz2*cy*cz-sz2*sz,  cz2*cy*sz+sz2*cz,           -cz2*sy]
        %     [ -sz2*cy*cz-cz2*sz, -sz2*cy*sz+cz2*cz,            sz2*sy]
        %     [             sy*cz,             sy*sz,                cy]
        
        [r1,r2,r3] = twoaxisrot( dcm(3,2,:), dcm(3,1,:), dcm(3,3,:), ...
                                 dcm(2,3,:), -dcm(1,3,:), ...
                                -dcm(2,1,:), dcm(2,2,:));
                
    case 'zxy'
        %     [ cy*cz-sy*sx*sz, cy*sz+sy*sx*cz,         -sy*cx]
        %     [         -sz*cx,          cz*cx,             sx]
        %     [ sy*cz+cy*sx*sz, sy*sz-cy*sx*cz,          cy*cx]

        [r1,r2,r3] = threeaxisrot( -dcm(2,1,:), dcm(2,2,:), dcm(2,3,:), ...
                                   -dcm(1,3,:), dcm(3,3,:), ...
                                    dcm(1,2,:), dcm(1,1,:));

    case 'zxz'
        %     [  cz2*cz-sz2*cx*sz,  cz2*sz+sz2*cx*cz,            sz2*sx]
        %     [ -sz2*cz-cz2*cx*sz, -sz2*sz+cz2*cx*cz,            cz2*sx]
        %     [             sz*sx,            -cz*sx,                cx]

        [r1,r2,r3] = twoaxisrot( dcm(3,1,:), -dcm(3,2,:), dcm(3,3,:), ...
                                 dcm(1,3,:), dcm(2,3,:), ...
                                 dcm(1,2,:), dcm(1,1,:));

    case 'yxz'
        %     [  cy*cz+sy*sx*sz,           sz*cx, -sy*cz+cy*sx*sz]
        %     [ -cy*sz+sy*sx*cz,           cz*cx,  sy*sz+cy*sx*cz]
        %     [           sy*cx,             -sx,           cy*cx]

        [r1,r2,r3] = threeaxisrot( dcm(3,1,:), dcm(3,3,:), -dcm(3,2,:), ...
                                   dcm(1,2,:), dcm(2,2,:), ...
                                  -dcm(1,3,:), dcm(1,1,:));
       
    case 'yxy'
        %     [  cy2*cy-sy2*cx*sy,            sy2*sx, -cy2*sy-sy2*cx*cy]
        %     [             sy*sx,                cx,             cy*sx]
        %     [  sy2*cy+cy2*cx*sy,           -cy2*sx, -sy2*sy+cy2*cx*cy]
     
        [r1,r2,r3] = twoaxisrot( dcm(2,1,:), dcm(2,3,:), dcm(2,2,:), ...
                                 dcm(1,2,:), -dcm(3,2,:), ...
                                -dcm(1,3,:), dcm(1,1,:));
      
    case 'yzx'
        %     [           cy*cz,              sz,          -sy*cz]
        %     [ -sz*cx*cy+sy*sx,           cz*cx,  sy*cx*sz+cy*sx]
        %     [  cy*sx*sz+sy*cx,          -cz*sx, -sy*sx*sz+cy*cx]
        
        [r1,r2,r3] = threeaxisrot( -dcm(1,3,:), dcm(1,1,:), dcm(1,2,:), ...
                                   -dcm(3,2,:), dcm(2,2,:), ...
                                    dcm(3,1,:), dcm(3,3,:));
        
    case 'yzy'
        %     [  cy2*cz*cy-sy2*sy,            cy2*sz, -cy2*cz*sy-sy2*cy]
        %     [            -cy*sz,                cz,             sy*sz]
        %     [  sy2*cz*cy+cy2*sy,            sy2*sz, -sy2*cz*sy+cy2*cy]

        [r1,r2,r3] = twoaxisrot( dcm(2,3,:), -dcm(2,1,:), dcm(2,2,:), ...
                                 dcm(3,2,:), dcm(1,2,:), ...
                                 dcm(3,1,:), dcm(3,3,:));

    case 'xyz'
        %     [          cy*cz, sz*cx+sy*sx*cz, sz*sx-sy*cx*cz]
        %     [         -cy*sz, cz*cx-sy*sx*sz, cz*sx+sy*cx*sz]
        %     [             sy,         -cy*sx,          cy*cx]

        [r1,r2,r3] = threeaxisrot( -dcm(3,2,:), dcm(3,3,:), dcm(3,1,:), ...
                                   -dcm(2,1,:), dcm(1,1,:), ...
                                    dcm(2,3,:), dcm(2,2,:));
        
    case 'xyx'
        %     [                cy,             sy*sx,            -sy*cx]
        %     [            sx2*sy,  cx2*cx-sx2*cy*sx,  cx2*sx+sx2*cy*cx]
        %     [            cx2*sy, -sx2*cx-cx2*cy*sx, -sx2*sx+cx2*cy*cx]

        [r1,r2,r3] = twoaxisrot( dcm(1,2,:), -dcm(1,3,:), dcm(1,1,:), ...
                                 dcm(2,1,:), dcm(3,1,:), ...
                                 dcm(2,3,:), dcm(2,2,:));
        
    case 'xzy'
        %     [          cy*cz, sz*cx*cy+sy*sx, cy*sx*sz-sy*cx]
        %     [            -sz,          cz*cx,          cz*sx]
        %     [          sy*cz, sy*cx*sz-cy*sx, sy*sx*sz+cy*cx]

        [r1,r2,r3] = threeaxisrot( dcm(2,3,:), dcm(2,2,:), -dcm(2,1,:), ...
                                   dcm(3,1,:), dcm(1,1,:), ...
                                  -dcm(3,2,:), dcm(3,3,:));
        
    case 'xzx'
        %     [                cz,             sz*cx,             sz*sx]
        %     [           -cx2*sz,  cx2*cz*cx-sx2*sx,  cx2*cz*sx+sx2*cx]
        %     [            sx2*sz, -sx2*cz*cx-cx2*sx, -sx2*cz*sx+cx2*cx]

        [r1,r2,r3] = twoaxisrot( dcm(1,3,:), dcm(1,2,:), dcm(1,1,:), ...
                                 dcm(3,1,:), -dcm(2,1,:), ...
                                -dcm(3,2,:), dcm(3,3,:));

    otherwise
        error(message('aero:dcm2angle:unknownRotation', type));
end

r1 = r1(:);
r2 = r2(:);
r3 = r3(:);

    function [r1,r2,r3] = threeaxisrot(r11, r12, r21, r31, r32, r11a, r12a)
        % find angles for rotations about X, Y, and Z axes
        r1 = atan2( r11, r12 );
        r21(r21 < -1) = -1;
        r21(r21 > 1) = 1;
        r2 = asin( r21 );
        r3 = atan2( r31, r32 );
        if strcmpi( lim, 'zeror3')
            for i = find(abs( r21 ) == 1.0)
                r1(i) = atan2( r11a(i), r12a(i) );
                r3(i) = 0;
            end
        end
    end

    function [r1,r2,r3] = twoaxisrot(r11, r12, r21, r31, r32, r11a, r12a)
        r1 = atan2( r11, r12 );
        r21(r21 < -1) = -1;
        r21(r21 > 1) = 1;
        r2 = acos( r21 );
        r3 = atan2( r31, r32 );
        if strcmpi( lim, 'zeror3')
            for i = find(abs( r21 ) == 1.0)
                r1(i) = atan2( r11a(i), r12a(i) );
                r3(i) = 0;
            end
        end
    end
end
