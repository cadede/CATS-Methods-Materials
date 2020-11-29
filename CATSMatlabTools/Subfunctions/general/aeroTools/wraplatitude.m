function [lat_wrapped, lat, lon] = wraplatitude( lat, lon, units)
% WRAPLATITUDE internal function to check and fix angle wrapping in
% latitude and longitude if needed.

%   Copyright 2007-2013 The MathWorks, Inc.

% Convert units if necessary
if strcmp(units,'deg')
    pideg = 180;
else
    pideg = pi;
end

% Set lat_wrapped flag to false
lat_wrapped = false;

% Re-set latitudes to values between pi and -pi
flat = abs(lat);

if any(flat>pideg)
    lat(flat>pideg) = mod(lat(flat>pideg)+pideg,2*pideg)-pideg;
    % Set lat_wrapped flag to true if necessary
    flat = abs(lat);
    lat_wrapped = true;
end

% Determine if any latitudes need to be wrapped
idx = flat>pideg/2;

if any(idx)
    
    % Adjustments for -pi to pi
    flat = abs(lat);
    latp2 = flat>pideg/2;
    lon(idx) = lon(idx) + pideg;
    lat(latp2) = sign(lat(latp2)).*(pideg/2-(flat(latp2)-pideg/2));
    
    % Set lat_wrapped flag to true if necessary
    lat_wrapped = true;
    
end