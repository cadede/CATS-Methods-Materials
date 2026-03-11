function [xEast, yNorth, zUp] = latlon2local(lat, lon, alt, origin)
%latlon2local Convert from geographic to local Cartesian coordinates
%  [xEast, yNorth, zUp] = latlon2local(lat, lon, alt, origin) 
%  transforms point locations from (lat, lon, alt) in degrees and
%  meters to local Cartesian coordinates (xEast, yNorth, zUp) in meters.
%  The local coordinate system is anchored at origin specified as
%  [latOrigin, lonOrigin, altOrigin] vector. Local x, y, z coordinates 
%  are lined up with east, north and up directions respectively.
%  alt and altOrigin are altitudes as returned by a typical GPS sensor.
%
%  Notes
%  -----
%  - The geographic coordinate system's latitude and longitude assume use
%    of WGS84 standard that is commonly used by the GPS.
%  - Altitude in this function is defined as height in meters above
%    WGS84 reference ellipsoid. 
%  - Some GPS receivers may use standards other than WGS84. Conversions 
%    using other ellipsoids are in the Mapping Toolbox. This function 
%    addresses the most common conversion between geographic locations and
%    Cartesian coordinates used by the vehicle's on-board sensors.
%
%  Example
%  -------
%  % Load GPS route
%  d = load('geoRoute.mat');
%
%  % Convert route to Cartesian coordinates
%  alt = 10;  % 10 meters is an approximate altitude in Boston, MA
%  origin = [d.latitude(1), d.longitude(1), alt];
%  [xEast, yNorth] = latlon2local(d.latitude, d.longitude, alt, origin);
%
%  % Visualize route's shape
%  figure;
%  plot(xEast, yNorth)
%  axis('equal'); % set 1:1 aspect ratio to see real-world shape
%
%  See also local2latlon, geoplayer, geoplot

% Copyright 2019 The MathWorks, Inc.

% Do only nominal error checking to avoid performance penalty associated
% with exhaustive input parsing.
validateattributes(origin, {'double','single'}, {'real', 'finite', 'numel', 3}, 'latlon2local', 'origin', 4);

[xEast, yNorth, zUp] = driving.internal.MapUtils.geodetic2enu(lat, lon, alt,...
    origin(1), origin(2), origin(3));

