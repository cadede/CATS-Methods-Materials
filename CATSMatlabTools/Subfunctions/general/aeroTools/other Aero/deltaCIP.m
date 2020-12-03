function varargout = deltaCIP(MJD,varargin)
%   DELTACIP Calculate Celestial Intermediate Pole (CIP) location
%   adjustment
%
%   DCIP = DELTACIP( UTC ) calculates the adjustment to location of the
%   Celestial Intermediate Pole (CIP) for a specific Universal Coordinated
%   Time (UTC), specified as a modified Julian date. By default, this
%   function uses a pre-populated list of IAU 2000A Earth orientation (IERS)
%   data. This list contains measured and calculated (predicted) data.
%
%   [DCIP, DCIPError] = DELTACIP( UTC ) calculates the adjustment to the
%   location of the CIP and the error for the adjustment to the location of
%   the CIP. 
%
%   [DCIP] = DELTACIP( UTC, Name, Value ) calculates the adjustment to the
%   location of the CIP using additional options specified by one or more
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
%   DCIP: M-by-2 array of the adjustment [dDeltaX,dDeltaY] location of
%   the Celestial Intermediate Pole according to the International
%   Astronomical Union (IAU) 2000A resolutions, in radians.
%
%   DCIPError: M-by-2 array of the error for the adjustment of the
%   location of the Celestial Intermediate Pole according to the
%   International Astronomical Union (IAU) 2000A resolutions, in radians.
%
%   Examples:
%
%   Calculate the CIP adjustment for December 28, 2015.
%
%   mjd = mjuliandate( 2015, 12, 28 )
%   dCIP = deltaCIP( mjd )
%
%   Calculate the CIP adjustment and CIP adjustment error for December 28,
%   2015 and January 10, 2016 using the file aeroiersdata.mat.
%
%   mjd = mjuliandate( [2015 12 28;2016 1 10] )
%   [dCIP, dCIPErr] = deltaCIP( mjd, 'Source', 'aeroiersdata.mat' )
%
%
%   See also AEROREADIERSDATA, DCMECI2ECEF, LLA2ECI, ECI2LLA, ECI2AER,
%   MJULIANDATE, DELTAUT1, POLARMOTION.

% Copyright 2018 The MathWorks, Inc.

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

% Initialize output
out = zeros(length(MJD),2);
if errFlag
    outerr = zeros(length(MJD),2);
end

% Use final index for DCIP vector in stored data
enddxy = length(dataIERS.dxy);

for k=1:length(MJD)
    % Check if the date is within the range of the shipped data
    if MJD(k)>=dataIERS.mjd(enddxy)+1
        % This is the case when the selected modified Julian date is
        % later than 1992/1/1 but it is not listed in the data.
        switch action
            case 'error'
                error(message('aero:aerout1utc:tablePrediction',...
                    char(p.Results.source),num2str(MJD(k))));
            case 'warning'
                warning(message('aero:aerout1utc:tablePrediction',...
                    char(p.Results.source),num2str(MJD(k))));
        end
        
        out(k,:) = dataIERS.dxy(enddxy,:);
        if errFlag
            outerr(k,:) = dataIERS.dxyerr(enddxy,:);
        end
    elseif MJD(k)<dataIERS.mjd(1)
        % This is the case for dates previous to 1992/1/1.
        switch action
            case 'error'
                error(message('aero:aerout1utc:earlyPrediction',num2str(MJD(k))));
            case 'warning'
                warning(message('aero:aerout1utc:earlyPrediction',num2str(MJD(k))));
        end
        out(k,:) = dataIERS.dxy(1,:);
        if errFlag
            outerr(k,:) = dataIERS.dxyerr(1,:);
        end
    else
        % This is the case for the data available in the .MAT file
        if MJD(k)<dataIERS.mjd(enddxy)+1 && MJD(k)>=dataIERS.mjd(enddxy)
            % Special case for the last day of the file.
            out(k,:) = dataIERS.dxy(enddxy,:);
            ipFlag = dataIERS.dxyip{enddxy};
            if errFlag
                outerr(k,:) = dataIERS.dxyerr(enddxy,:);
            end
        else
            idx = find((MJD(k)<dataIERS.mjd));
            ipFlag = dataIERS.dxyip{idx(1)-1};
            out(k,:) = dataIERS.dxy(idx(1)-1,:);
            if errFlag
                outerr(k,:) = dataIERS.dxyerr(idx(1)-1,:);
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
% Convert from milliseconds of arc to radians
if errFlag
    varargout{1} = out*(pi/(648000*1000));
    varargout{2} = outerr*(pi/(648000*1000));
else
    varargout{1} = out*(pi/(648000*1000));
end
