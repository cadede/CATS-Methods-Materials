function varargout = atmoshwm07(latitude,longitude,altitude,varargin)
%  ATMOSHWM07 Implement Horizontal Wind Model 07
%  W = ATMOSHWM07( LATITUDE, LONGITUDE, ALTITUDE, 'ADDPARAMNAME', ADDPARAMVALUE)
%  implements the U.S. Naval Reseach Laboratory HWM(TM) routine to
%  calculate the meridional and zonal components of the wind for a given 
%  set of geophysical data, ALTITUDE, LATITUDE, LONGITUDE, APINDEX, DAY,
%  and SECONDS.
%
%  Inputs for ATMOSHWM07 are:
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
%  'MODEL':     String to indicate the type of model used. It be one of the
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
%
%  Output calculated by ATMOSHWM07 is:
%  W:   M-by-2 array of meridional wind (northward) and the zonal wind
%       (eastward) values, in m/s.
%
%  Examples:
%  Calculate the total horizontal wind model for an altitude of 25000 m
%  above mean sea level (msl) for latitude of 45 degrees south, longitude
%  of 85 degrees west, the 150th day of the year at 11 am UTC, and an AP
%  index of 80.
%
%  w = atmoshwm07(-45,-85,25000,'day',150,'seconds',39600,'apindex',80,'model','total')
%
%  Calculate the quiet horizontal wind model for altitudes of 100000 m and
%  150000 m above msl at latitude 50 degrees north, longitude 20 degrees
%  west, at midnight UTC of January 30th.
%
%  w = atmoshwm07([50;50],[-20;-20],[100000;150000],'day',[30;30])
%
%
%  See also ATMOSCOESA, ATMOSNRLMSISE00, ATMOSCIRA.

%  Copyright 2014-2016 The MathWorks, Inc.

%  References:
%  [1] Drob, D.P., Emmert, J.T., Crowley, G., Picone, J.M., Shepherd, G.G.,
%  Skinner, W., Hays, P., Niciejewski, R.J., Larsen, M., She, C.Y.,
%  Meriwether, J.W., Hernandez, G., Jarvis, Martin J., Sipler, D.P.,
%  Tepley, C.A., O'Brien, M.S., Bowman, J.R., Wu, Q., Murayama, Y.,
%  Kawamura, S., Reid, I.M., Vincent, R.A. An empirical model of the
%  Earth's horizontal wind fields: HWM07. Journal of Geophysical Research,
%  2008.

warning(message('aero:atmoshwm07:updateWarning'));
%% Validate i/o:
% Validate number of inputs
narginchk(3,13);
% Validate number of outputs
nargoutchk(0,2);
%% Call the generalized Horizontal Wind Model Function
w = atmoshwm(latitude,longitude,altitude,'version','07',varargin{:});
%% Outputs
if nargout == 2
    varargout{1} = w(:,1);
    varargout{2} = w(:,2);
else
    varargout{1} = w;
end