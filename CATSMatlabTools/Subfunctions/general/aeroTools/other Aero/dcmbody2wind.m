function dcm = dcmbody2wind( alpha, beta )
%  DCMBODY2WIND Convert angle of attack and sideslip angle to direction cosine matrix.
%   N = DCMBODY2WIND( A, B ) calculates the direction cosine matrix, N,
%   for given set of angle of attack and sideslip angle, A, B.   A is an M array of
%   angles of attack.  B is an M array of sideslip angles.  N returns an 3-by-3-by-M
%   matrix containing M direction cosine matrices.  N performs the
%   coordinate transformation of a vector in body-axes into a vector in
%   wind-axes.  Angles of attack and sideslip angles are input in radians.  
%
%   Examples:
%
%   Determine the direction cosine matrix from angle of attack and sideslip angle:
%      alpha = 0.4363; 
%      beta = 0.1745;
%      dcm = dcmbody2wind( alpha, beta )
%
%   Determine the direction cosine matrix from multiple angle of attack and sideslip angle:
%      alpha = [0.4363 0.1745]; 
%      beta = [0.1745 0.0873];
%      dcm = dcmbody2wind( alpha, beta )
%
%   See also ANGLE2DCM, DCM2ALPHABETA, DCM2ANGLE.


%   Copyright 2000-2010 The MathWorks, Inc.

if any(~isreal(alpha) || ~isnumeric(alpha))
    error(message('aero:dcmbody2wind:isNotReal1'));
end

if any(~isreal(beta) || ~isnumeric(beta))
    error(message('aero:dcmbody2wind:isNotReal2'));
end

if (length(alpha) ~= length(beta))
    error(message('aero:dcmbody2wind:wrongDimension'));
end

angles = [alpha(:) beta(:)];

dcm = zeros(3,3,size(angles,1));
cang = cos( angles );
sang = sin( angles );

dcm(1,1,:) = cang(:,2).*cang(:,1);
dcm(1,2,:) = sang(:,2);
dcm(1,3,:) = sang(:,1).*cang(:,2);
dcm(2,1,:) = -sang(:,2).*cang(:,1);
dcm(2,2,:) = cang(:,2);
dcm(2,3,:) = -sang(:,1).*sang(:,2);
dcm(3,1,:) = -sang(:,1);
dcm(3,2,:) = 0.0;
dcm(3,3,:) = cang(:,1);
