function N = geoidegm96( lat, lon, varargin )
%  GEOIDEGM96 Implement the EGM96 Geopotential Model
%   N = GEOIDEGM96( LAT, LON ) calculates the geoid height as determined
%   from the EGM96 Geopotential Model. Geoid heights are interpolated from
%   a 15-minute grid of point values in the tide-free system, using the
%   EGM96 Geopotential Model to degree and order 360. The geoid undulations
%   are with respect to the WGS84 ellipsoid.
%
%   The GEOIDEGM96 function calculates geoid heights to 0.01 meters. The
%   output is of the same data type as the input LAT (see below for LAT
%   data type restrictions).
% 
%   The interpolation scheme wraps over the poles to allow for geoid height
%   calculations at and near these locations.
%
%   Inputs for GEOIDEGM96 are:
%   LAT:    an array of M geocentric latitudes in degrees where north latitude
%           is positive, and south latitude is negative. LAT must be of type
%           single or double. If LAT is not in the range of [-90,90] it is
%           wrapped to be within the range.
%   LON:    an array of M geocentric longitude in degrees where east
%           longitude is positive, west is negative. LON must be of type
%           single or double. If LON is not in the range of [0,360] it is
%           wrapped to be within the range.
%   ACTION: a string to determine action for out of range latitude or
%           longitude. Specify if out of range input invokes a 'Warning',
%           'Error', or no action ('None'). The default is 'Warning'.
%
%   Limitations:
%
%   This function has the limitations of the 1996 Earth Geopotential Model
%   For more information see the documentation.
%
%   The WGS84 EGM96 geoid undulations have an error range of +/- 0.5 to
%   +/- 1.0 meters worldwide.
%
%   Examples:
%
%   Calculate the geoid height at 42.4 degrees N latitude and 71.0 degrees 
%   W longitude:
%       N = geoidegm96( 42.4, -71.0)
%
%   Calculate the geoid height at two different locations.
%       N = geoidegm96( [39.3,33.4], [-77.2, 36.5])
%
%   See also GRAVITYWGS84, GRAVITYSPHERICALHARMONIC

%   Copyright 1990-2011 The MathWorks, Inc.

%   References:
%   
%   NIMA TR8350.2: "Department of Defense World Geodetic System
%   1984, Its Definition and Relationship with Local Geodetic Systems."
%
%   NASA/TP-1998-206861: "The Development of the Joint NASA GSFC and NIMA
%   Geopotential Model EGM96"
%
%   National Geospatial-Intelligence Agency website:
%   http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html

warning(message('aero:geoidegm96:ObsoleteFunction'));

persistent geoidegm96data

narginchk(2, 3);

actionstr = {'none' 'warning' 'error'};

checkinputs();

action = 'warning';

if nargin == 3
    evalaction( varargin{1} );
end

actionidx = find(strcmpi(actionstr, action));

checklatitude();
checklongitude();

if isempty(geoidegm96data)
    geoidegm96data = load('geoidegm96grid.mat');
end

% checkinputs above should ensure lat and lon are of equal dimension
n = numel(lat);

% Preallocate
N  = zeros(1,n,'single');

halfWindowSize = geoidegm96data.windowSize/2;

% Get lat data type to set output later
latType = class(lat);

% Convert to single for interpolation
if ~strcmp(latType,'single')
    lat = single(lat);
    lon = single(lon);
end

% Loop through all of the lon/lat input pairs
for k = 1:n
    % Find indices for window of height data. If there is an exact match,
    % it will be in the middle of the window (i.e. vector will be of size
    % windowSize+1-by-windowsize+1 with matching index in the middle).
    lonbpi = [find(geoidegm96data.lonbp<lon(k),halfWindowSize,'last'), ...
        find(geoidegm96data.lonbp==lon(k)), ...
        find(geoidegm96data.lonbp>lon(k),halfWindowSize,'first')];
    latbpi = [find(geoidegm96data.latbp<lat(k),halfWindowSize,'last'), ...
        find(geoidegm96data.latbp==lat(k)), ...
        find(geoidegm96data.latbp>lat(k),halfWindowSize,'first')];

    % As mentioned above, if there are windowSize+1 elements in lonbpi, the
    % middle value is the index for an exact value of lon.
    if numel(lonbpi) == geoidegm96data.windowSize+1
        % Pull out the lat window values for the exact lon value.
        latsforlon = geoidegm96data.grid(latbpi,lonbpi(geoidegm96data.windowSize/2 + 1));

        % Otherwise perform interpolation to find height values for given lon
        % at each lat in window.
    else
        % This transpose is here because interp1 will interp down each
        % column of the input matrix, and in this instance we want it to
        % interp across each row of grid(latbpi,lonbpi).
        latsforlon = interp1(geoidegm96data.lonbp(lonbpi),geoidegm96data.grid(latbpi,lonbpi)',lon(k),'spline');
    end

    % As mentioned above, if there are windowSize+1 elements in latbpi, the
    % middle value is the index for an exact value of lat.
    if numel(latbpi) == geoidegm96data.windowSize+1
        % Pull out the height value for exact value of lat from range of
        % height possibilities determined above from lon.
        N(k) = latsforlon(halfWindowSize + 1);

        % Otherwise perform interpolation to find height value.
    else
        N(k) = interp1(geoidegm96data.latbp(latbpi),latsforlon,lat(k),'spline');
    end
end

% Expand to double for unit changes and truncating.
N = double(N);

N = fix(N.*100)*0.01;

% Cast to desired output data type
N = cast(N,latType);

%%
    function checkinputs()
        if ~isnumeric(lat)
            % Latitude should be a numeric array.  Otherwise error.
            error(message('aero:geoidegm96:latitudeNotNumeric'));
        end
        if ~isnumeric(lon)
            % Altitude should be a numeric array.  Otherwise error.
            error(message('aero:geoidegm96:longitudeNotNumeric'));
        end
        if ~all(size(lat) == size(lon))
            error(message('aero:geoidegm96:arraySize'));
        end
        if ~isfloat(lat) || ~isfloat(lon)
            error(message('aero:geoidegm96:notFloat'));
        end
        if ~strcmpi(class(lat),class(lon))
            error(message('aero:geoidegm96:differentInputType'));
        end
    end

%%
    function checklatitude()
        % Wrap latitude and longitude if necessary
        [latwrapped,lat,lon] = wraplatitude(lat,lon,'deg');
        if latwrapped
            % Create function handle array to handle messages based on
            % action
            fhlat = {@() disp('')... % do nothing for 'none'
                  @() warning(message('aero:geoidegm96:warnLatitudeWrap')), ...
                  @() error(message('aero:geoidegm96:latitudeWrap'))};
            
           % Call appropriate function
           fhlat{actionidx}()
        end
    end

%%
    function checklongitude()
        [lonwrapped, lon] = wraplongitude(lon,'deg','360');

        if lonwrapped
            % Create function handle array to handle messages based on
            % action
            fhlon = {@() disp('')... % do nothing for 'none'
                @() warning(message('aero:geoidegm96:warnLongitudeWrap')), ...
                @() error(message('aero:geoidegm96:longitudeWrap'))};

            % Call appropriate function
            fhlon{actionidx}()
        end
    end
%%
    function evalaction( input )
        if any(strcmpi(actionstr, input))
            action = lower( input );
        else
            error(message('aero:geoidegm96:unknownStringAction',input));
        end
    end
%%
end
