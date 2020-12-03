function q = dcm2quat( dcm, varargin )
%  DCM2QUAT Convert direction cosine matrix to quaternion.
%   Q = DCM2QUAT( N ) calculates the quaternion, Q, for given direction
%   cosine matrix, N. N is a 3-by-3-by-M matrix containing M orthogonal
%   direction cosine matrices. N performs the coordinate transformation of
%   a vector in inertial axes to a vector in body axes. Q returns an M-by-4
%   matrix containing M quaternions. Q has its scalar number as the first
%   column.
% 
%   Q = DCM2QUAT( N, ACTION ) uses ACTION for error handling during
%   rotation matrix validation. Specify if an invalid rotation matrix
%   invokes a 'Warning', 'Error', or no action ('None'). The default is
%   'None'.
%
%   Q = DCM2QUAT( N, ACTION, TOL ) uses TOL as a relative error tolerance
%   for rotation matrix validation. It is a scalar value. Rotation matrix
%   validation confirms that the matrix is orthogonal [transpose(N) * N ==
%   1 +/- TOL] and proper [det(N) == 1 +/- TOL]. The default value is
%   eps(2).
%
%   Examples:
%
%   Determine the quaternion from direction cosine matrix:
%      dcm = [0 1 0; 1 0 0; 0 0 -1];
%      q = dcm2quat(dcm)
%
%   Determine the quaternions from multiple direction cosine matrices:
%      dcm        = [ 0 1 0; 1 0 0; 0 0 -1]; 
%      dcm(:,:,2) = [ 0.4330    0.2500   -0.8660; ...
%                     0.1768    0.9186    0.3536; ...
%                     0.8839   -0.3062    0.3536];
%      q = dcm2quat(dcm)
%
%   See also ANGLE2DCM, DCM2ANGLE, QUAT2DCM, QUAT2ANGLE, ANGLE2QUAT.

%   Copyright 2000-2017 The MathWorks, Inc.

narginchk(1, 3);

%Validate inputs
if nargin < 2
    validateDCM(dcm);
elseif nargin < 3
    validateDCM(dcm, varargin{1});
else
    validateDCM(dcm, varargin{1}, varargin{2});
end

for i = size(dcm,3):-1:1

    q(i,4) =  0; 
    
    tr = trace(dcm(:,:,i));

    if (tr > 0)
        sqtrp1 = sqrt( tr + 1.0 );
        
        q(i,1) = 0.5*sqtrp1; 
        q(i,2) = (dcm(2, 3, i) - dcm(3, 2, i))/(2.0*sqtrp1);
        q(i,3) = (dcm(3, 1, i) - dcm(1, 3, i))/(2.0*sqtrp1); 
        q(i,4) = (dcm(1, 2, i) - dcm(2, 1, i))/(2.0*sqtrp1); 
    else
        d = diag(dcm(:,:,i));
        if ((d(2) > d(1)) && (d(2) > d(3)))
            % max value at dcm(2,2,i)
            sqdip1 = sqrt(d(2) - d(1) - d(3) + 1.0 );
            
            q(i,3) = 0.5*sqdip1; 
            
            if ( sqdip1 ~= 0 )
                sqdip1 = 0.5/sqdip1;
            end
            
            q(i,1) = (dcm(3, 1, i) - dcm(1, 3, i))*sqdip1; 
            q(i,2) = (dcm(1, 2, i) + dcm(2, 1, i))*sqdip1; 
            q(i,4) = (dcm(2, 3, i) + dcm(3, 2, i))*sqdip1; 
        elseif (d(3) > d(1))
            % max value at dcm(3,3,i)
            sqdip1 = sqrt(d(3) - d(1) - d(2) + 1.0 );
            
            q(i,4) = 0.5*sqdip1; 
            
            if ( sqdip1 ~= 0 )
                sqdip1 = 0.5/sqdip1;
            end
            
            q(i,1) = (dcm(1, 2, i) - dcm(2, 1, i))*sqdip1;
            q(i,2) = (dcm(3, 1, i) + dcm(1, 3, i))*sqdip1; 
            q(i,3) = (dcm(2, 3, i) + dcm(3, 2, i))*sqdip1; 
        else
            % max value at dcm(1,1,i)
            sqdip1 = sqrt(d(1) - d(2) - d(3) + 1.0 );
            
            q(i,2) = 0.5*sqdip1; 
            
            if ( sqdip1 ~= 0 )
                sqdip1 = 0.5/sqdip1;
            end
            
            q(i,1) = (dcm(2, 3, i) - dcm(3, 2, i))*sqdip1; 
            q(i,3) = (dcm(1, 2, i) + dcm(2, 1, i))*sqdip1; 
            q(i,4) = (dcm(3, 1, i) + dcm(1, 3, i))*sqdip1; 
        end
    end
end
