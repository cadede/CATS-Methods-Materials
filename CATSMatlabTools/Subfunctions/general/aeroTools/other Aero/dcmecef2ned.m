function dcm = dcmecef2ned( lat, lon )
%  DCMECEF2NED Convert geodetic latitude and longitude to direction cosine matrix.
%   N = DCMECEF2NED( LAT, LON ) calculates the direction cosine matrix, N,
%   for given set of geodetic latitude and longitude, LAT, LON.   LAT is an M array of
%   geodetic latitudes.  LON is an M array of longitudes.  N returns an 3-by-3-by-M
%   matrix containing M direction cosine matrices.  N performs the
%   coordinate transformation of a vector in Earth-centered Earth-fixed
%   (ECEF) axes into a vector in north-east-down (NED) axes.  Geodetic
%   latitudes and longitudes are input in degrees.  
%
%   Examples:
%
%   Determine the direction cosine matrix from geodetic latitude and longitude:
%      lat = 45; 
%      lon = -122;
%      dcm = dcmecef2ned( lat, lon )
%
%   Determine the direction cosine matrix from multiple geodetic latitude and longitude:
%      lat = [45 37.5]; 
%      lon = [-122 -85];
%      dcm = dcmecef2ned( lat, lon )
%
%   See also ANGLE2DCM, DCM2ANGLE, DCM2LATLON.


%   Copyright 2000-2010 The MathWorks, Inc.

if any(~isreal(lat) || ~isnumeric(lat))
    error(message('aero:dcmecef2ned:isNotReal1'));
end

if any(~isreal(lon) || ~isnumeric(lon))
    error(message('aero:dcmecef2ned:isNotReal2'));
end

if (length(lat) ~= length(lon))
    error(message('aero:dcmecef2ned:wrongDimension'));
end

if (abs(lat)>90)
    error(message('aero:dcmecef2ned:exceed90'));
end

if (abs(lon)>180)
    error(message('aero:dcmecef2ned:exceed180'));
end

angles = convang( [lat(:) lon(:)] ,'deg','rad');

dcm = zeros(3,3,size(angles,1));
cang = cos( angles );
sang = sin( angles );

dcm(1,1,:) = -cang(:,2).*sang(:,1);
dcm(1,2,:) = -sang(:,2).*sang(:,1);
dcm(1,3,:) = cang(:,1);
dcm(2,1,:) = -sang(:,2);
dcm(2,2,:) = cang(:,2);
dcm(2,3,:) = 0.0;
dcm(3,1,:) = -cang(:,2).*cang(:,1);
dcm(3,2,:) = -sang(:,2).*cang(:,1);
dcm(3,3,:) = -sang(:,1);
