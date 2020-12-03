function obj = fganimation
%  FGANIMATION Construct FlightGear animation object.
%   H = FGANIMATION constructs a FlightGear animation object. The
%   FlightGear animation object is returned to H.
%
%   Example:
%
%   Construct a FlightGear animation object, h:
%      h = fganimation
%
%   See Also GENERATERUNSCRIPT, Aero.FlightGearAnimation.Play

% Copyright 1990-2007 The MathWorks, Inc.

%   Object Properties
%       TimeseriesSource           Specify variable that contains
%                                  timeseries data.
%       TimeseriesSourceType       Specify the type of timeseries data
%                                  stored in 'TimeseriesSource'. The
%                                  default value is 'Array6DoF'.
%          'Timeseries'            MATLAB timeseries data with 6 values per
%                                  time: LAT LON ALT PHI THETA PSI. Values
%                                  are linearly interpolated vs. time using
%                                  interp1.
%          'StructureWithTime'     Simulink struct with time (Simulink root
%                                  outport logging 'Structure with time'):
%                                  - signals(1).values:  LAT LON ALT 
%                                  - signals(2).values:  PHI THETA PSI
%                                  - Signals are linearly interpolated vs.
%                                  time using interp1.
%          'Array6DoF'             double precision array in N rows and 7
%                                  columns for 6-DoF data: TIME LAT LON ALT
%                                  PHI THETA PSI.  Values are linearly
%                                  interpolated vs. time using interp1.
%          'Array3DoF'             double precision array in N rows and 4
%                                  columns for 3-DoF data: TIME LAT ALT
%                                  THETA.  Values are linearly
%                                  interpolated vs. time using interp1.
%          'Custom'                Position and angle data is retrieved
%                                  from 'TimeseriesSource' by the currently
%                                  registered 'TimeseriesReadFcn'
%       TimeseriesReadFcn          Specify a function to read a 'Custom'
%                                  'TimeseriesSourceType'.
%       TimeScaling                Specify the seconds of animation data
%                                  per second of wall-clock time. The
%                                  default ratio is 1.
%       FramesPerSecond            Specify the number of frames per second
%                                  used to animate the 'TimeseriesSource'.
%                                  The default value is 12 frames per
%                                  second.
%       FlightGearVersion          Select your FlightGear software version:
%                                  '0.9.3', '0.9.8', '0.9.9' or '0.9.10'. 
%                                  The default version is '0.9.10'.
%       OutputFileName             Specify the name of the output file. The
%                                  file name is the name of the command you
%                                  will use to start FlightGear with these
%                                  initial parameters. The default value is
%                                  'runfg.bat'.
%       FlightGearBaseDirectory    Specify the name of your FlightGear
%                                  installation directory. The default
%                                  value is 'C:\Program Files\FlightGear'.
%       GeometryModelName          Specify the name of the folder
%                                  containing the desired model geometry in
%                                  the FlightGear\data\Aircraft directory.
%                                  The default value is 'HL20'.
%       DestinationIpAddress       Specify your destination IP address. The
%                                  default value is '127.0.0.1'.
%       DestinationPort            Specify your network flight dynamics
%                                  model (fdm) port. This destination port
%                                  should be an unused port that you can
%                                  use when you launch FlightGear. The
%                                  default value is '5502'.
%       AirportId                  Specifies the airport ID. The list of
%                                  supported airports is available in the
%                                  FlightGear interface, under Location.
%                                  The default value is 'KSFO'.
%       RunwayId                   Specify the runway ID. The default value
%                                  is '10L'.
%       InitialAltitude            Specify the initial altitude of the
%                                  aircraft, in feet. The default value is
%                                  7224 feet.
%       InitialHeading             Specify the initial heading of the
%                                  aircraft, in degrees. The default value
%                                  is 113 degrees. 
%       OffsetDistance             Specify the offset distance of the
%                                  aircraft from the airport, in miles. The
%                                  default value is 4.72 miles. 
%       OffsetAzimuth              Specify the offset azimuth of the
%                                  aircraft, in degrees. The default value
%                                  is 0 degrees.

obj = Aero.FlightGearAnimation;
