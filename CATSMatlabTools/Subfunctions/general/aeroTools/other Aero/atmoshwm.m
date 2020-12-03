function varargout = atmoshwm(latitude,longitude,altitude,varargin)
%  ATMOSHWM Implement Horizontal Wind Model
%  W = ATMOSHWM( LATITUDE, LONGITUDE, ALTITUDE, 'ADDPARAMNAME', ADDPARAMVALUE)
%  implements the U.S. Naval Research Laboratory HWM(TM) routine to
%  calculate the meridional and zonal components of the wind for a given 
%  set of geophysical data, ALTITUDE, LATITUDE, LONGITUDE, APINDEX, DAY,
%  and SECONDS.
%
%  Inputs for ATMOSHWM are:
%  LATITUDE:    M-by-1 array of geodetic latitude values, in degrees.
%  LONGITUDE:   M-by-1 array of geodetic longitude values, in degrees.
%  ALTITUDE:    M-by-1 array of geopotential height values,in meters.
%               Altitude values must be within the range 0 to 500 km.
%
%  Input parameter Name,Value pairs are:
%  'APINDEX':   M-by-1 array of 3 hour geomagnetic disturbance index
%               values. Ap index values must be within the range 0 to
%               400. Default value is an M-by-1 null array.
%  'DAY':       M-by-1 array of day of the year in Universal Coordinated
%               Time. Specify the day as a value between 1 and 366 (for a
%               leap year). Default value is an M-by-1 array set to 1.
%  'SECONDS':   M-by-1 array of elapsed seconds since midnight for the
%               selected day, in Universal Coordinated Time. Default value
%               is an M-by-1 null array.
%  'MODEL':     String to indicate the type of model used. It can be one of the
%               following: 
%               'QUIET':        Calculates the horizontal wind model
%                               without the magnetic disturbances. Default.
%               'DISTURBANCE': 	Calculates the effect of only magnetic
%                               disturbances in the wind. For this model,
%                               specify APINDEX values greater or equal to
%                               0.
%               'TOTAL':        Calculates the combined effect of the quiet
%                               and magnetic disturbances. For this model,
%                               specify APINDEX values greater or equal to
%                               0.
%  'ACTION':    String to determine action if inputs are out of range.
%               Specify if out-of-range input invokes a 'WARNING', 'ERROR',
%               or no action ('NONE'). The default is 'ERROR'.
%  'VERSION':   String to determine the version of the wind model. Specify
%               '14' for Horizontal Wind Model 14 or '07' for Horizontal
%               Wind Model 07. The default value is '14'.
%
%  Output calculated by ATMOSHWM is:
%  W:   M-by-2 array of meridional wind (northward) and the zonal wind
%       (eastward) values, in m/s.
%
%  Examples:
%  Calculate the total horizontal wind model for an altitude of 25000 m
%  above mean sea level (msl) for latitude of 45 degrees south, longitude
%  of 85 degrees west, the 150th day of the year at 11 am UTC, an AP
%  index of 80, and for version 14.
%
%  w = atmoshwm(-45,-85,25000,'day',150,'seconds',39600,'apindex',80,'model','total','version','14')
%
%  Calculate the quiet horizontal wind model 14 for altitudes of 100000 m and
%  150000 m above msl at latitude 50 degrees north, longitude 20 degrees
%  west, at midnight UTC of January 30th.
%
%  w = atmoshwm([50;50],[-20;-20],[100000;150000],'day',[30;30])
%
%
%  See also ATMOSCOESA, ATMOSNRLMSISE00, ATMOSCIRA.

%  Copyright 2016-2018 The MathWorks, Inc.

%  References:
%  [1] Drob, D.P., Emmert, J.T., Crowley, G., Picone, J.M., Shepherd, G.G.,
%  Skinner, W., Hays, P., Niciejewski, R.J., Larsen, M., She, C.Y.,
%  Meriwether, J.W., Hernandez, G., Jarvis, Martin J., Sipler, D.P.,
%  Tepley, C.A., O'Brien, M.S., Bowman, J.R., Wu, Q., Murayama, Y.,
%  Kawamura, S., Reid, I.M., Vincent, R.A. An empirical model of the
%  Earth's horizontal wind fields: HWM07. Journal of Geophysical Research,
%  2008.
%  [2] Drob, D.P.,Emmert, J.T.,Meriwether, J.,Makela, J.J., Doornbos, E., 
%  Conde, M., Hernandez, G.,Noto, J., Zawdie, K.A., McDonald, S.E., Huba,
%  J. D., Klenzing, J. H. An update to the Horizontal Wind Model (HWM): The
%  quiet time thermosphere. Earth Space and Science. 2015.

%% Validate i/o:
% Validate number of inputs
narginchk(3,15);
% Validate number of outputs
nargoutchk(0,2);
% Validate length
len = length(altitude);
% Parse and validate inputs
p = inputParser;
addRequired(p,'latitude',@(x) validateattributes(x,{'numeric'},...
    {'real','nonnan','size',[len,1]}));
addRequired(p,'longitude',@(x) validateattributes(x,{'numeric'},...
    {'real','nonnan','size',[len,1]}));
addRequired(p,'altitude',@(x) validateattributes(x,{'numeric'},{'real',...
    'nonnan'}));
addParameter(p,'day',ones(len,1),@(x) validateattributes(x,{'numeric'},...
    {'real','nonnan','size',[len,1]}));
addParameter(p,'apindex',zeros(len,1),@(x) validateattributes(x,{'numeric'},...
    {'real','nonnan','size',[len,1]}));
addParameter(p,'seconds',zeros(len,1),@(x) validateattributes(x,{'numeric'},...
    {'real','nonnan','size',[len,1]}));
addParameter(p,'model','quiet',@(x) validateattributes(x,{'char','string'},{'nonempty'}));
addParameter(p,'action','error',@(x) validateattributes(x,{'char','string'},{'nonempty'}));
addParameter(p,'version','14',@(x) validateattributes(x,{'char','string'},{'nonempty'}));
parse(p,latitude,longitude,altitude,varargin{:});
windModel = lower(p.Results.model);
validModelList = {'quiet','disturbance','total'};
validModel = validatestring(windModel,validModelList);
actionInput = p.Results.action;
validActionList = {'None','Warning','Error'};
validAction = lower(validatestring(actionInput,validActionList));
versionInput = lower(p.Results.version);
validVersionList = {'07','14'};
validVersion = validatestring(versionInput,validVersionList);
% Out of range actions
switch validAction
    case 'error'
        if any(altitude<0) || any(altitude>500e3)
            error(message('aero:atmoshwm07:altitudeRange'));
        end
        if any(p.Results.day<1) || any(p.Results.day>366)
            error(message('aero:atmoshwm07:dayRange'));
        end
        if (any(p.Results.apindex>400) || any(p.Results.apindex<0)) && ...
            (strcmp(validModel,'total') || strcmp(validModel,'disturbance'))
            error(message('aero:atmoshwm07:apIndexRange'));
        end
        if any(p.Results.seconds<0) || any(p.Results.seconds>86400)
            error(message('aero:atmoshwm07:secondsRange'));
        end
        if ~any(strcmp(p.UsingDefaults,'apindex')) && strcmp(validModel,'quiet')
            error(message('aero:atmoshwm07:quietAPError'));
        end
        apIndex = p.Results.apindex;
        altitude = p.Results.altitude;
    case 'warning'
        if any(altitude<0) || any(altitude>500e3)
            warning(message('aero:atmoshwm07:altitudeRange'));
        end
        if any(p.Results.day<1) || any(p.Results.day>366)
            warning(message('aero:atmoshwm07:dayRange'));
        end
        if (any(p.Results.apindex>400) || any(p.Results.apindex<0)) && ...
            (strcmp(validModel,'total') || strcmp(validModel,'disturbance'))
            warning(message('aero:atmoshwm07:apIndexRange'));
        end
        if any(p.Results.seconds<0) || any(p.Results.seconds>86400)
            warning(message('aero:atmoshwm07:secondsRange'));
        end
        if ~any(strcmp(p.UsingDefaults,'apindex')) && strcmp(validModel,'quiet')
            warning(message('aero:atmoshwm07:quietAPWarning'));
        end
        % If Ap index values are negative set to 0
        apIndex = p.Results.apindex;
        apIndex(p.Results.apindex<0) = 0;
        % If altitude values are negative set to 0
        altitude = p.Results.altitude;
        altitude(p.Results.altitude<0) = 0;
        
    case 'none'
        % If Ap index values are negative set to 0
        apIndex = p.Results.apindex;
        apIndex(p.Results.apindex<0) = 0;
        % If altitude values are negative set to 0
        altitude = p.Results.altitude;
        altitude(p.Results.altitude<0) = 0;
end
% Initialize output matrix
w = zeros(len,2);
% Make sure the days are whole numbers
intDay = floor(p.Results.day);
% Make sure the latitude and longitudes are within -/+ 90 and +/- 180
% respectively
[~, lat, lon] = wraplatitude( latitude, longitude,'deg');
[~, lon] = wraplongitude( lon, 'deg', '180' );
switch validVersion
    case '07'
        %% Call the hwm07 method
        switch validModel
            case 'quiet'
                for k=1:len
                    w(k,:)=hwm07methods(intDay(k),p.Results.seconds(k),altitude(k),...
                        lat(k),lon(k),-1,0);
                end
            case 'total'
                for k=1:len
                    w(k,:)=hwm07methods(intDay(k),p.Results.seconds(k),altitude(k),...
                        lat(k),lon(k),apIndex(k),1);
                end                
            case 'disturbance'
                for k=1:len
                    w(k,:)=hwm07methods(intDay(k),p.Results.seconds(k),altitude(k),...
                        lat(k),lon(k),apIndex(k),2);
                end
        end
    case '14'
        %% Call the hwm14 method
        switch validModel
            case 'quiet'
                for k=1:len
                    w(k,:)=hwm14methods(intDay(k),p.Results.seconds(k),altitude(k),...
                        lat(k),lon(k),-1,0);
                end
            case 'total'
                for k=1:len
                    w(k,:)=hwm14methods(intDay(k),p.Results.seconds(k),altitude(k),...
                        lat(k),lon(k),apIndex(k),1);
                end                
            case 'disturbance'
                for k=1:len
                    w(k,:)=hwm14methods(intDay(k),p.Results.seconds(k),altitude(k),...
                        lat(k),lon(k),apIndex(k),2);
                end
        end
end
%% Outputs
if nargout == 2
    varargout{1} = w(:,1);
    varargout{2} = w(:,2);
else
    varargout{1} = w;
end
