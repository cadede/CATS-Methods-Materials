function [ lat, lon ] = dcm2latlon( dcm, varargin ) 
%  DCM2LATLON Convert direction cosine matrix to geodetic latitude and longitude.
%   [LAT, LON] = DCM2LATLON( N ) calculates the geodetic latitude and
%   longitude, LAT, LON, for given direction cosine matrix, N. N is a
%   3-by-3-by-M matrix containing M orthogonal direction cosine matrices. N
%   performs the coordinate transformation of a vector in Earth-centered
%   Earth-fixed (ECEF) axes into a vector in north-east-down (NED) axes.
%   LAT is an M array of geodetic latitudes.  LON is an M array of
%   longitudes. Geodetic latitudes and longitudes are output in degrees.
%
%   [LAT, LON] = DCM2LATLON( N, ACTION ) uses ACTION for error handling
%   during rotation matrix validation. Specify if an invalid rotation
%   matrix invokes a 'Warning', 'Error', or no action ('None'). The default
%   is 'None'. 
%
%   [LAT, LON] = DCM2LATLON( N, ACTION, TOL ) uses TOL as a relative error
%   tolerance for rotation matrix validation. It is a scalar value.
%   Rotation matrix validation confirms that the matrix is orthogonal
%   [transpose(N) * N == 1 +/- TOL] and proper [det(N) == 1 +/- TOL]. The
%   default value is eps(2).
%
%   Examples:
%
%   Determine the geodetic latitude and longitude from direction cosine matrix:
%      dcm = [ 0.3747    0.5997    0.7071; ...
%              0.8480   -0.5299         0; ...
%              0.3747    0.5997   -0.7071]; 
%      [lat, lon] = dcm2latlon( dcm )
%
%   Determine the geodetic latitude and longitude from multiple direction cosine matrices:
%      dcm =        [ 0.3747    0.5997    0.7071; ...
%                     0.8480   -0.5299         0; ...
%                     0.3747    0.5997   -0.7071]; 
%      dcm(:,:,2) = [-0.0531    0.6064    0.7934; ...
%                     0.9962    0.0872         0; ...
%                    -0.0691    0.7903   -0.6088]; 
%      [lat, lon] = dcm2latlon( dcm )
%
%   See also ANGLE2DCM, DCM2ANGLE, DCMECEF2NED.


%   Copyright 2000-2018 The MathWorks, Inc.

narginchk(1, 3);

%Validate inputs
if nargin < 2
    validateDCM(dcm);
elseif nargin < 3
    validateDCM(dcm, varargin{1});
else
    validateDCM(dcm, varargin{1}, varargin{2});
end

lat = convang( asin(-dcm(3,3,:)),'rad','deg');
lon = convang( atan2(-dcm(2,1,:), dcm(2,2,:)),'rad','deg');

lat = lat(:);
lon = lon(:);
