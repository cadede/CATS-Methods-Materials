function [lon_wrapped, lon] = wraplongitude( lon, inputUnits, maxRange )
% WRAPLONGITUDE internal function to check and fix angle wrapping in
% longitude.

%   Copyright 2007-2010 The MathWorks, Inc.

% Convert units if necessary
if strcmp(inputUnits,'deg')
    lon = convang(lon,'deg','rad');
end

lon_wrapped = false;

switch maxRange
    % Range goes between -pi and pi
    case {'180','pi'}
        idx = lon > pi | lon < -pi;
        if any(idx)
            lon(idx) = rem(lon(idx),2*pi)- (2*pi)*fix(rem(lon(idx),2*pi)/pi);
            lon_wrapped = true;
        end
    % Range goes between 0 and 2*pi
    case {'360','2*pi'}
        idx = lon > 2*pi | lon < 0;
        if any(idx)
            lon(idx) = mod(lon(idx),2*pi);
            lon_wrapped = true;
        end
    otherwise
        % do nothing
end

% Convert units if necessary
if strcmp(inputUnits,'deg')
    lon = convang(lon,'rad','deg');
end