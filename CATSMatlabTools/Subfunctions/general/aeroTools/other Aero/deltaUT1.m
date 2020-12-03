function varargout = deltaUT1(MJD,varargin)
%   DELTAUT1 Calculate difference between Coordinated Universal Time (UTC)
%   and Principal Universal Time (UT1)
%
%   DUT1 = DELTAUT1( UTC ) calculates the difference between Coordinated
%   Universal Time (UTC) and Principal Universal Time (UT1) for UTC,
%   specified as a modified Julian date. By default, this function uses a
%   pre-populated list of IAU 2000A Earth orientation (IERS) data. This list
%   contains measured and calculated (predicted) data.
%
%   [DUT1, DUT1Error] = DELTAUT1( UTC ) calculates the difference between
%   UTC and UT1 and the error for the diference between UTC and UT1.
%
%   [DUT1] = DELTAUT1( UTC, Name, Value ) calculates the difference
%   between UTC and UT1 using additional options specified by one or more
%   Name,Value pair arguments.
%
%   Input Arguments:
%   UTC: M-by-1 array of UTC dates represented as modified Julian dates.
%   Use the mjuliandate function to convert the UTC date to a modified
%   Julian date.
%
%   Name-Value Pair Arguments:
%   'SOURCE': MAT-file containing custom list of Earth orientation data,
%   specified as a string. Use the aeroReadIERSData function to create an
%   up-to-date list of Earth orientation data. The folder for this file
%   must be on the MATLAB path.
%   By default, this function uses the pre-populated list contained in
%   matlab\toolbox\aero\aero\aeroiersdata.mat.
%
%   'ACTION': Action to take in case of out-of-range or predicted value
%   dates, specified as a string. Specify if the action should be
%   'Error', 'Warning', or 'None'. Errors and warnings display the cause of
%   the action (out-of-range or predicted value dates). The default is
%   'Warning'.
%
%   Output Arguments:
%   DUT1: M-by-1 array of the difference between UT1 and UTC (UT1-UTC)
%   according to the International Astronomical Union (IAU) 2000A
%   resolutions, in seconds.
%
%   DUT1Error: M-by-1 array of the error for the difference between UT1 and
%   UTC (UT1-UTC) according to the International Astronomical Union (IAU)
%   2000A resolutions, in seconds.
%
%   Examples:
%
%   Calculate the UT1-UTC value for December 28, 2015.
%
%   mjd = mjuliandate(  2015, 12, 28 )
%   dUT1 = deltaUT1( mjd )
%
%   Calculate the UT1-UTC value for December 28, 2015 and January 10, 2016
%   using the file aeroiersdata.mat.
%
%   mjd = mjuliandate( [2015 12 28;2016 1 10] )
%   [dUT1, dUT1Err] = deltaUT1( mjd, 'Source', 'aeroiersdata.mat' )
%
%
%   See also AEROREADIERSDATA, DCMECI2ECEF, LLA2ECI, ECI2LLA, ECI2AER,
%   MJULIANDATE, DELTACIP, POLARMOTION.

% Copyright 2017-2018 The MathWorks, Inc.

%   References: 
%   [1] IERS website. http://www.iers.org/ 
%   [2] United States Naval Observatory website. http://www.usno.navy.mil/

% Validate input
narginchk(1,5);
p=inputParser;
addRequired(p,'MJD',@(x) validateattributes(x,{'numeric'},{'real','nonempty'}));
addParameter(p,'action','warning',@(x) validateattributes(x,{'char',...
    'string'},{'nonempty'}));
addParameter(p,'source','aeroiersdata.mat',@(x) validateattributes(x,{'char',...
    'string'},{'nonempty'}));
parse(p,MJD,varargin{:});
action = p.Results.action;
action = lower(validatestring(action,{'None','Warning','Error'}));
validPath = which(p.Results.source,'-all');
if isempty(validPath) && ~isfile(p.Results.source)
    error(message('aero:aerout1utc:invalidSource',char(p.Results.source)));
end

% Validate output
nargoutchk(0,2);
errFlag= false;
if nargout==2
    errFlag = true;
end

% Load data
dataIERS = load(p.Results.source);

% Initialize data
out = zeros(size(MJD));
if errFlag
    outerr = zeros(size(MJD));
end

% Use final index for UT1-UTC vector in stored data
endut1 = length(dataIERS.ut1utc);

for k=1:length(MJD)
    % Check if the date is within the range of the shipped data
    if MJD(k)>=dataIERS.mjd(endut1)+1
        % This is the case when the selected modified Julian date is
        % later than 1992/1/1 but it is not listed in the data.
        switch action
            case 'error'
                error(message('aero:aerout1utc:latePrediction',num2str(MJD(k))));
            case 'warning'
                warning(message('aero:aerout1utc:latePrediction',num2str(MJD(k))));
        end
        % Calculate Besselian date
        bessDate = 2000+(MJD(k)-51544.03)/365.2422;
        % Calculate UT2-UT1
        dUT2UT1 = 0.022*sin(2*pi*bessDate)-0.012*cos(2*pi*bessDate)...
            -0.006*sin(4*pi*bessDate)+0.007*cos(4*pi*bessDate);
        % Calculate UT1-UTC
        out(k) = 0.5382-0.00124*(MJD(k)-57801)-dUT2UT1;
        % Saturate at +/- 0.9 seconds according to definition of UT1.
        if abs(out(k))>0.9
            out(k) = sign(out(k))*0.9;
        end
        if errFlag
            outerr(k) = NaN;
        end
    elseif MJD(k)<dataIERS.mjd(1)
        % This is the case for dates previous to 1992/1/1.
        switch action
            case 'error'
                error(message('aero:aerout1utc:earlyPrediction',num2str(MJD(k))));
            case 'warning'
                warning(message('aero:aerout1utc:earlyPrediction',num2str(MJD(k))));
        end
        out(k) = dataIERS.ut1utc(1);
        if errFlag
            outerr(k) = dataIERS.ut1utcerr(1);
        end
    else
        % This is the case for the data available in the .MAT file
        if MJD(k)<dataIERS.mjd(end)+1 && MJD(k)>=dataIERS.mjd(end)
            % Special case for the last day of the file.
            out(k) = dataIERS.ut1utc(endut1);
            ipFlag = dataIERS.ut1utcip{endut1};
            if errFlag
                outerr(k) = dataIERS.ut1utcerr(endut1);
            end
        else
            idx = find((MJD(k)<dataIERS.mjd));
            ipFlag = dataIERS.ut1utcip{idx(1)-1};
            out(k) = dataIERS.ut1utc(idx(1)-1);
            if errFlag
                outerr(k) = dataIERS.ut1utcerr(idx(1)-1);
            end
        end
        if strcmp(ipFlag,'P')
            switch action
                case 'error'
                    error(message('aero:aerout1utc:tablePrediction',...
                        char(p.Results.source),num2str(MJD(k))));
                case 'warning'
                    warning(message('aero:aerout1utc:tablePrediction',...
                        char(p.Results.source),num2str(MJD(k))));
            end
        end
    end
end
if errFlag
    varargout{1} = out;
    varargout{2} = outerr;
else
    varargout{1} = out;
end
