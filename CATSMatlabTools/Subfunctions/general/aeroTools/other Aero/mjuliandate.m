function mjd = mjuliandate( varargin ) 
%  MJULIANDATE Calculate modified Julian date.
%   MJD = MJULIANDATE( V ) converts one or more date vectors V into
%   modified Julian date MJD.  Input V can be an M-by-6 or M-by-3 matrix
%   containing M full or partial date vectors, respectively.  MJULIANDATE
%   returns a column vector of M modified Julian dates.  Modified Julian
%   dates begin at midnight rather than noon and have the first two digits of 
%   the corresponding Julian date removed.
%
%   A date vector contains six elements, specifying year, month, day, hour, 
%   minute, and second. A partial date vector has three elements, specifying 
%   year, month, and day.  Each element of V must be a positive double 
%   precision number.  
%   
%   MJD = MJULIANDATE(S,F) converts one or more date strings S to modified
%   Julian date MJD using format string F. S can be a character array where
%   each row corresponds to one date string, or one dimensional cell array
%   of strings.  MJULIANDATE returns a column vector of M Julian dates,
%   where M is the number of strings in S. 
%   
%   All of the date strings in S must have the same format F, which must be
%   composed of date format symbols according to DATESTR help. Formats with
%   'Q' are not accepted by MJULIANDATE. 
%
%   Certain formats may not contain enough information to compute a date
%   number.  In those cases, hours, minutes, and seconds default to 0, days
%   default to 1, months default to January, and years default to the
%   current year. Date strings with two character years are interpreted to
%   be within the 100 years centered around the current year.
%   
%   MJD = MJULIANDATE(Y,MO,D) and MJD = MJULIANDATE([Y,MO,D]) return the
%   modified Julian date for corresponding elements of the Y,MO,D
%   (year,month,day) arrays. Y, MO, and D must be either one dimensional 
%   arrays of the same length or scalar values. 
%   
%   MJD = MJULIANDATE(Y,MO,D,H,MI,S) and MJD = MJULIANDATE([Y,MO,D,H,MI,S])
%   return the modified Julian dates for corresponding elements of the
%   Y,MO,D,H,MI,S (year,month,day,hour,minute,second) arrays.  The six
%   arguments must be either one dimensional arrays of the same length or 
%   scalar values. 
%   
%   Limitations: 
%   This function is valid for all common era (CE) dates in the Gregorian
%   calendar.
%
%   The calculation of Julian date does not take into account leap seconds.
%   
%   Examples:
%   
%   Calculate modified Julian date for May 24, 2005:
%	   mjd = mjuliandate('24-May-2005','dd-mm-yyyy')
%
%   Calculate modified Julian date for December 19, 2006:
%	   mjd = mjuliandate(2006,12,19)
%
%   Calculate modified Julian date for October 10, 2004 at 12:21:00 pm:
%	   mjd = mjuliandate(2004,10,10,12,21,0)  
%
%   See also JULIANDATE.

%   Copyright 2000-2016 The MathWorks, Inc.

% Second correction
if nargin == 1 %Case for a vector with 6 elements
    if size(varargin{:},2)==6
        datemat = varargin{:};
        %Get second value
        sec = datemat(:,6);
        %Remove mantissa
        modSec= mod(sec,1);
        datemat(:,6) = floor(sec);
        %Take care of carryover values
        [year,month,day,hour,min,sec] = datevec(datenum(datemat));
        %Add the mantissa back
        sec = sec + modSec;
    else
        [year,month,day,hour,min,sec] = datevec(datenum(varargin{:}));
    end
elseif nargin == 2 % Entry is char array with date and format
    [year,month,day,hour,min,sec] = datevec(datenum(char(varargin{1}),...
        char(varargin{2})));
elseif nargin ==6 %Case for 6 element input
    %Get second value
    sec = varargin{6};
    %Remove mantissa
    modSec = mod(sec,1);
    varargin{6} = floor(sec);
    %Take care of carryover values
    [year,month,day,hour,min,sec] = datevec(datenum(varargin{:}));
    %Add the mantissa back
    sec = sec + modSec;
else %Any other case
    [year,month,day,hour,min,sec] = datevec(datenum(varargin{:}));
end

% Check that the input of the function is either a one dimensional vector
% or a scalar.
validateattributes(year,{'numeric'},{'vector'})

% Check for non-BCE values
if any(year<1)
    warning(message('aero:aeroephemerides:invalidCalendarBCE'));
end

for k = length(month):-1:1
    if ( month(k) <= 2 ) % January & February
        year(k)  = year(k) - 1.0;
        month(k) = month(k) + 12.0;
    end
end

dayFraction= (hour + min/60 + sec/3600)/24;

day = floor( 365.25*year ) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
    floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day + ...
     - 679006.0;
 
mjd = dayFraction + day;

