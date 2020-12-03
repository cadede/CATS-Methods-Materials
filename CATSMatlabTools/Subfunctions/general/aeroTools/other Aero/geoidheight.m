function N = geoidheight( lat, lon, varargin )
%  GEOIDHEIGHT Implement a geopotential model to calculate geoid height
%   N = GEOIDHEIGHT( LAT, LON, MODEL ) calculates the geoid height as
%   determined from a selected geopotential model, MODEL. An array of M
%   geoid heights, N, are interpolated at M geocentric latitude, LAT, and M
%   geocentric longitude, LON, from a grid of point values in the tide-free
%   system, using the selected geopotential model, MODEL.  
%
%   The GEOIDHEIGHT function calculates geoid heights to 0.01 meters for
%   EGM96 and custom.  The GEOIDHEIGHT function calculates geoid heights to
%   0.001 meters for EGM2008.
% 
%   The interpolation scheme wraps over the poles to allow for geoid height
%   calculations at and near these locations.
%
%   Alternate formats for calling geoid height are:
%   N = GEOIDHEIGHT( LAT, LON )
%   N = GEOIDHEIGHT( LAT, LON, ACTION )
%   N = GEOIDHEIGHT( LAT, LON, MODEL, ACTION )
%   N = GEOIDHEIGHT( LAT, LON, 'Custom', DATAFILE )
%   N = GEOIDHEIGHT( LAT, LON, 'Custom', DATAFILE, ACTION )
%
%   Inputs for GEOIDHEIGHT are:
%   LAT      :an array of M geocentric latitudes in degrees where north
%             latitude is positive, and south latitude is negative. LAT
%             must be of type single or double. If LAT is not in the range
%             of [-90,90] it is wrapped to be within the range.
%   LON      :an array of M geocentric longitude in degrees where east
%             longitude is positive, west is negative. LON must be of type
%             single or double. If LON is not in the range of [0,360] it is
%             wrapped to be within the range.
%   MODEL    :a string specifying the geopotential model: 'EGM2008' (Earth),
%             'EGM96' (Earth), or 'Custom'.  The default is 'EGM96'.  
%             'EGM96' uses a 15-minute grid of point values in the
%             tide-free system, using EGM96 Geopotential Model to degree
%             and order 360.  The EGM2008 MODEL uses a 2.5-minute grid of
%             point values in the tide-free system, using the EGM2008
%             Geopotential Model to degree and order 2159. The geoid
%             undulations for EGM96 and EGM2008 are with respect to the
%             WGS84 ellipsoid. 
%   DATAFILE :a mat-file containing an array of geocentric
%             latitude breakpoints, 'latbp', an array of geocentric
%             longitude breakpoints, 'lonbp', a table of geoid height
%             values,'grid' and a even integer scalar greater than 2 for
%             number of interpolation points, 'windowSize'.  This is only
%             needed for a 'Custom' geopotential model.     
%   ACTION   :a string to determine action for out of range latitude or
%             longitude. Specify if out of range input invokes a 'Warning',
%             'Error', or no action ('None'). The default is 'Warning'.
%
%   Output for GEOIDHEIGHT is:
%   N        :an array of M geoid heights in meters with the same data type
%             as the input LAT.
%
%   Limitations:
%
%   This function using the 'EGM96' MODEL has the limitations of the 1996
%   Earth Geopotential Model.  For more information see the documentation.
%
%   The WGS84 EGM96 geoid undulations have an error range of +/- 0.5 to
%   +/- 1.0 meters worldwide.
%
%   This function using the 'EGM2008' MODEL has the limitations of the 2008
%   Earth Geopotential Model.  For more information see the documentation.
%
%   Examples:
%
%   Calculate the EGM96 geoid height at 42.4 degrees N latitude and 71.0 degrees 
%   W longitude with warning actions:
%       N = geoidheight( 42.4, -71.0 )
%
%   Calculate the EGM2008 geoid height at two different locations with
%   error actions.
%       N = geoidheight( [39.3, 33.4], [77.2, 36.5], 'egm2008','error')
%
%   Calculate a custom geoid height at two different locations with
%   no actions.
%       N = geoidheight( [39.3, 33.4], [-77.2, 36.5], 'custom', ...
%           'geoidegm96grid','none')
%
%   Note: This function uses geoid data that can be obtained using the
%   aeroDataPackage command.
%
%   See also GRAVITYWGS84, GRAVITYSPHERICALHARMONIC

%   Copyright 2010-2018 The MathWorks, Inc.

%   References:
%   [1] NIMA TR8350.2: "Department of Defense World Geodetic System
%       1984, Its Definition and Relationship with Local Geodetic Systems."
%   [2] NASA/TP-1998-206861: "The Development of the Joint NASA GSFC and NIMA
%       Geopotential Model EGM96"
%   [3] Pavlis, N.K., S.A. Holmes, S.C. Kenyon, and J.K. Factor, "An Earth
%       Gravitational Model to Degree 2160: EGM2008", presented at the 2008
%       General Assembly of the European Geosciences Union, Vienna,
%       Austria, April 13-18, 2008. 
%   [4] Vallado, D. A., "Fundamentals of Astrodynamics and Applications",
%       McGraw-Hill, New York, 1997.  
%
%   National Geospatial-Intelligence Agency websites:
%   http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
%   http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html

persistent astgeoiddata astgeoidCustomError

narginchk(2, 5);

checkinputs();

% set default values
model  = 'EGM96';
action = 'warning';
datafile = 'geoidegm96grid.mat';

switch nargin                                                              
    case 3
        if ~ischar( varargin{1} ) && ~isstring( varargin{1} )
            error(message('aero:geoidheight:inputTypeVar3'));
        end
        if strcmpi( varargin{1}, 'custom')
            narginchk(4, 5);
        else
            % check that the 2008 model is in the path	 
            checkegm2008( varargin{1} );
            % assign model or action
            modeloraction( varargin{1} );
        end       
    case 4
        if (~ischar( varargin{1} ) && ~isstring( varargin{1} ) || ...
                (~ischar( varargin{2} ) && ~isstring( varargin{2} )))
            error(message('aero:geoidheight:inputTypeVar4'));
        end
        if strcmpi( varargin{1}, 'custom')
            % custom model with datafile
            model = lower( varargin{1} );
            datafile = varargin{2};
        else
            % check that the 2008 model is in the path	 
            checkegm2008( varargin{1} );
            % assign model and action
            checkmodel( varargin{1} );
            checkaction( varargin{2} );
        end
    case 5
        if (~ischar( varargin{1} ) && ~isstring( varargin{1} )) || ...
                (~ischar( varargin{2} ) && ~isstring( varargin{2} )) || ...
                (~ischar( varargin{3} ) && ~isstring( varargin{3} ))
            error(message('aero:geoidheight:inputTypeVar5'));
        end
        % set action for custom model with datafile
        if strcmpi( varargin{1}, 'custom')
            % custom model with datafile
            model = lower( varargin{1} );
            datafile = varargin{2};
            checkaction( varargin{3} );
        else
            % This option is only for 'custom'
            error(message('aero:geoidheight:wrongModel'));
        end
end

actionstr = {'none' 'warning' 'error'};
actionidx = find(strcmp(actionstr, action));

checklatitude();
checklongitude();

if ( isempty(astgeoiddata) || ~strcmp( astgeoiddata.type, model ) || ...
                   ~strcmpi( astgeoiddata.file, datafile ) || ~isempty(astgeoidCustomError))
    % data needs to be initialized
    datafile = char(datafile);
    try
        astgeoiddata = load(datafile);
    catch MECustomFileLoad
        throwAsCaller(MECustomFileLoad)
    end
    astgeoiddata.type = model;
    astgeoiddata.file = datafile;
    if strcmp('custom', astgeoiddata.type)
        % check for the existence of correct variables in mat-file
        fieldsExist = ~isfield(astgeoiddata,{'latbp' 'lonbp' 'grid' 'windowSize'});
        
        fhFields = { @() disp(''), ... % all fields exist
            @() error(message('aero:geoidheight:noLatitudeBP', datafile)), ...
            @() error(message('aero:geoidheight:noLongitudeBP', datafile)), ...
            @() error(message('aero:geoidheight:noGrid', datafile)), ...
            @() error(message('aero:geoidheight:noWindowSize', datafile))};
        
        idx = find(fieldsExist,1);
        if isempty(idx)
            idx = 0;
        end
        
        % allow same custom filename to run through check if had previous error
        astgeoidCustomError = 1;
        
        % Call appropriate function
        fhFields{idx+1}()    
        
        % check data type and sizes of custom data
        if  ~isnumeric( astgeoiddata.latbp ) || ~isnumeric( astgeoiddata.lonbp ) || ...
                ~isnumeric( astgeoiddata.grid ) || ~isnumeric( astgeoiddata.windowSize )
            error(message('aero:geoidheight:notNumeric'))
        end
        if  ~isscalar( astgeoiddata.windowSize )
            error(message('aero:geoidheight:notScalar'))
        end
        if  ~isvector( astgeoiddata.latbp ) || ~isvector( astgeoiddata.lonbp ) 
            error(message('aero:geoidheight:not1DArray'))
        end
         if ~all( size( astgeoiddata.grid ) == [ numel(astgeoiddata.latbp) numel(astgeoiddata.lonbp) ] ) 
            error(message('aero:geoidheight:wrongMatrixSize'))
        end
        if  any(~isfinite( astgeoiddata.latbp )) || any(~isfinite( astgeoiddata.lonbp )) || ...
             any(any(~isfinite( astgeoiddata.grid ))) || any(~isfinite( astgeoiddata.windowSize ))
            error(message('aero:geoidheight:notFinite'))
        end
        if (mod(astgeoiddata.windowSize,2) || (astgeoiddata.windowSize <= 2))
            error(message('aero:geoidheight:notEvenInteger'))
        end
        % successful custom file read
        astgeoidCustomError = [];
     end
end

% checkinputs above should ensure lat and lon are of equal dimension
n = numel(lat);

% Preallocate
N  = zeros(1,n,'single');

halfWindowSize = astgeoiddata.windowSize/2;

% Get lat data type to set output later
latType = class(lat);

% Convert to single for interpolation
if ~isa(latType,'single')
    lat = single(lat);
    lon = single(lon);
end

% Loop through all of the lon/lat input pairs
for k = 1:n
    % Find indices for window of height data. If there is an exact match,
    % it will be in the middle of the window (i.e. vector will be of size
    % windowSize+1-by-windowsize+1 with matching index in the middle).
    lonbpi = [find(astgeoiddata.lonbp<lon(k),halfWindowSize,'last'), ...
        find(astgeoiddata.lonbp==lon(k)), ...
        find(astgeoiddata.lonbp>lon(k),halfWindowSize,'first')];
    latbpi = [find(astgeoiddata.latbp<lat(k),halfWindowSize,'last'), ...
        find(astgeoiddata.latbp==lat(k)), ...
        find(astgeoiddata.latbp>lat(k),halfWindowSize,'first')];

    % As mentioned above, if there are windowSize+1 elements in lonbpi, the
    % middle value is the index for an exact value of lon.
    if numel(lonbpi) == astgeoiddata.windowSize+1
        % Pull out the lat window values for the exact lon value.
        latsforlon = astgeoiddata.grid(latbpi,lonbpi(astgeoiddata.windowSize/2 + 1));

        % Otherwise perform interpolation to find height values for given lon
        % at each lat in window.
    else
        % This transpose is here because interp1 will interp down each
        % column of the input matrix, and in this instance we want it to
        % interp across each row of grid(latbpi,lonbpi).
        latsforlon = interp1(astgeoiddata.lonbp(lonbpi),astgeoiddata.grid(latbpi,lonbpi)',lon(k),'spline');
    end

    % As mentioned above, if there are windowSize+1 elements in latbpi, the
    % middle value is the index for an exact value of lat.
    if numel(latbpi) == astgeoiddata.windowSize+1
        % Pull out the height value for exact value of lat from range of
        % height possibilities determined above from lon.
        N(k) = latsforlon(halfWindowSize + 1);

        % Otherwise perform interpolation to find height value.
    else
        N(k) = interp1(astgeoiddata.latbp(latbpi),latsforlon,lat(k),'spline');
    end
end

% Expand to double for unit changes and truncating.
N = double(N);

if strcmp(astgeoiddata.type,'egm2008')
    N = fix(N.*1000)*0.001;
else
    N = fix(N.*100)*0.01;
end

% Cast to desired output data type
N = cast(N,latType);

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkinputs()
        if ~isnumeric(lat)
            % Latitude should be a numeric array.  Otherwise error.
            error(message('aero:geoidheight:latitudeNotNumeric'));
        end
        if ~isnumeric(lon)
            % Altitude should be a numeric array.  Otherwise error.
            error(message('aero:geoidheight:longitudeNotNumeric'));
        end
        if ~all(size(lat) == size(lon))
            error(message('aero:geoidheight:arraySize'));
        end
        if ~isfloat(lat) || ~isfloat(lon)
            error(message('aero:geoidheight:notFloat'));
        end
        if ~strcmpi(class(lat),class(lon))
            error(message('aero:geoidheight:differentInputType'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checklatitude()
        % Wrap latitude and longitude if necessary
        [latwrapped,lat,lon] = wraplatitude(lat,lon,'deg');
        if latwrapped
            % Create function handle array to handle messages based on
            % action
            fhlat = {@() disp('')... % do nothing for 'none'
                  @() warning(message('aero:geoidheight:warnLatitudeWrap')), ...
                  @() error(message('aero:geoidheight:latitudeWrap'))};
            
           % Call appropriate function
           fhlat{actionidx}()
        end
    end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checklongitude()
        [lonwrapped, lon] = wraplongitude(lon,'deg','360');

        if lonwrapped
            % Create function handle array to handle messages based on
            % action
            fhlon = {@() disp('')... % do nothing for 'none'
                @() warning(message('aero:geoidheight:warnLongitudeWrap')), ...
                @() error(message('aero:geoidheight:longitudeWrap'))};

            % Call appropriate function
            fhlon{actionidx}()
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkmodel( str )
        switch lower( str )
            case 'egm2008'
                model = lower( str );
                datafile = 'geoidegm2008grid.mat';
            case 'egm96'
                model = lower( str );
                datafile = 'geoidegm96grid.mat';
            otherwise
                error(message('aero:geoidheight:unknownModel'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkaction( str )
        switch lower( str )
            case { 'error', 'warning', 'none' }
                action = lower( str );
            otherwise
                error(message('aero:geoidheight:unknownAction'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function modeloraction( str )
        switch lower( str )
            case 'egm2008'
                model = lower( str );
                datafile = 'geoidegm2008grid.mat';
            case 'egm96'
                model = lower( str );
                datafile = 'geoidegm96grid.mat';
            case { 'error', 'warning', 'none' }
                action = lower( str );
            otherwise
                error(message('aero:geoidheight:unknownString'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 function checkegm2008( str )	 
     if strcmp(str,'egm2008')
         if isempty(which('geoidegm2008grid.mat','-ALL'))
             error(message('aero:geoidheight:noDownloadData',...
                 '<a href="matlab:aeroDataPackage">aeroDataPackage</a>'))
         end
     end
 end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end
