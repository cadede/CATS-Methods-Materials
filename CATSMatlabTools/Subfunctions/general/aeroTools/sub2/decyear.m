function dyear = decyear( varargin )
%  DECYEAR Calculate decimal year.
%   DY = DECYEAR( V ) converts one or more date vectors V into decimal year
%   DY.  Input V can be an M-by-6 or M-by-3 matrix containing M full or
%   partial date vectors respectively.  DECYEAR returns a column vector of
%   M decimal years.
%
%	A date vector contains six elements, specifying year, month, day, hour, 
%	minute, and second. A partial date vector has three elements, specifying 
%	year, month, and day.  Each element of V must be a positive double 
%	precision number.  
%
%	DY = DECYEAR(S,F) converts one or more date strings S to decimal year
%	DY using format string F. S can be a character array where each
%	row corresponds to one date string, or one dimensional cell array of 
%	strings.  DECYEAR returns a column vector of M decimal years, where M is 
%	the number of strings in S. 
%
%	All of the date strings in S must have the same format F, which must be
%	composed of date format symbols according to Table 2 in DATESTR help.
%	Formats with 'Q' are not accepted by DECYEAR.  
%
%	Certain formats may not contain enough information to compute a date
%	number.  In those cases, hours, minutes, and seconds default to 0, days
%	default to 1, months default to January, and years default to the
%	current year. Date strings with two character years are interpreted to
%	be within the 100 years centered around the current year.
%
%	DY = DECYEAR(Y,MO,D) and DY = DECYEAR([Y,MO,D]) return the decimal
%	year for corresponding elements of the Y,MO,D (year,month,day)
%	arrays. Y, MO, and D must be either one dimensional arrays of the same 
%   length or scalar values. 
%
%       DY = DECYEAR(Y,MO,D,H,MI,S) and DY = DECYEAR([Y,MO,D,H,MI,S]) return the
%	decimal years for corresponding elements of the Y,MO,D,H,MI,S
%	(year,month,day,hour,minute,second) arrays.  The six arguments must be
%   either one dimensional arrays of the same length or scalar values.
%   Limitation: 
%
%   The calculation of decimal year does not take into account leap seconds.
%
%   Examples:
%
%   Calculate decimal year for May 24, 2005:
%	   dy = decyear('24-May-2005','dd-mm-yyyy')
%
%   Calculate decimal year for December 19, 2005:
%	   dy = decyear(2006,12,19)
%
%   Calculate decimal year for October 10, 2004 at 12:21:00 pm:
%	   dy = decyear(2004,10,10,12,21,0)  
%
%   See also JULIANDATE, LEAPYEAR, WRLDMAGM, IGRF11MAGM.

%   Copyright 2000-2012 The MathWorks, Inc.

% get serial date number
day = datenum( varargin{:} );

% Check that the input of the function is either a one dimensional vector
% or a scalar.
validateattributes(day,{'numeric'},{'vector'})

year = zeros(size(day));

% get date components from serial date
output = datevec( day );
year(:) = output(:,1);

ndays = 365*ones(size(year));
ndays(leapyear(year)) = 366;

dayofyear = day - datenum([year, 01*ones(size(year)), 01*ones(size(year)) ]);

% calculate the decimal year
dyear = year + dayofyear./ndays;
