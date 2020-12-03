function varargout = tdbjuliandate(TT)
%   TDBJULIANDATE Estimate the julian date for Barycentric Dynamical Time
%   [JDTDB,TTTDB] = TDBJULIANDATE(TT) estimates an M-by-1 array of the julian date
%   for the Barycentric Dynamical Time (TDB) based on the Terrestrial Time
%   (TT).
%
%   Input argument for TDBJULIANDATE is:
%   TT:     Array of Terrestrial Time (TT) in year, month, day, hour,
%           minutes and seconds for which the function calculates the
%           julian date for Barycentric Dynamical Time. Define the array as
%           one of the following: 
%           * Array with 1 row and 6 columns.
%           * M by 6 array for M julian dates, one for each TT date. 
%           Values for year, month, day, hour and minutes should be whole
%           numbers.
%
%   Outputs calculated by TDBJULIANDATE are:
%   JDTDB:  Julian date for the Barycentric Dynamical Time. Its size is M
%           by 1 with M rows, one for each Terrestrial Time input.
%   TTTDB:  Difference in seconds between Terrestrial Time and Barycentric
%           Dynamical Time (TT-TDB).Its size is M by 1 with M rows, one for
%           each Terrestrial Time input. 
%
%   Examples:
%
%   Estimate the julian date for the Barycentric Dynamical Time for the
%   Terrestrial Time 2014/10/15 16:22:31.
%   JDTDB = tdbjuliandate([2014,10,15,16,22,31])
%
%   Estimate the julian dates for the Barycentric Dynamical Time and TT-TDB
%   in seconds for the terrestrial time 2014/10/15 16:22:31 and 2010/7/22
%   1:57:17.
%   [JDTDB,TTTTDB] = tdbjuliandate([2014,10,15,16,22,31;2010,7,22,1,57,17])
%
%   Limitations:
%   This estimation method is valid for the year period from 1980 to 2050.
%
%   See also JULIANDATE, MJULIANDATE.

%   Copyright 2014 The MathWorks,Inc. 

%   References:
%   [1]United States Naval Observatory.The Astronomical Almanac. U.S.
%   Government Printing Office. Washington, D.C. 2015.
%   [2]Vallado, D., Seago, J., Seidelmann, P.K., Implementation Issues
%   Surrounding the New IAU Reference Systems for Astrodynamics. American
%   Astronautical Society Printing Office.,Tampa, FL. 2006.

%% Validate i/o

% Validate outputs
nargoutchk(0,2)

% Validate date
validateattributes(TT,{'numeric'},{'ncols',6,'real','finite','nonnan'})

% Assign date vectors
year = TT(:,1);
month = TT(:,2);
day = TT(:,3);
hour = TT(:,4);
min = TT(:,5);
sec = TT(:,6);

% Validate vectors
validateattributes(year,{'numeric'},{'>=',1980,'<=',2050});
validateattributes(month,{'numeric'},{'>=',1,'<=',12});
validateattributes(day,{'numeric'},{'>=',1,'<=',31});
validateattributes(hour,{'numeric'},{'>=',0,'<=',24});
validateattributes(min,{'numeric'},{'>=',0,'<=',60});
validateattributes(sec,{'numeric'},{'>=',0,'<=',60});

%% Calculate the difference in the mean longitudes of Earth and Jupiter.
% Calculate julian date for TT
jdTT = juliandate(TT);
g = 357.53 + 0.98560028*(jdTT-2451545);
deltaL = 246.11 + 0.90251792*(jdTT-2451545);

%% Calculate Barycentric Dynamical Time
TDBmTT = 0.001657*sind(g)+0.000022*sind(deltaL);
TDB = jdTT + TDBmTT/86400;

varargout{1} = TDB;
if nargout==2   
    varargout{2} = -TDBmTT;
end
% [EOF]
