function [T, a, P, rho] = atmosnonstd( h, atype, extreme, freq, varargin )
%  ATMOSNONSTD Use climatic data from MIL-STD-210 or MIL-HDBK-310.
%   [T, A, P, RHO] = ATMOSNONSTD(H, ATYPE, EXTREME, FREQ, EXTALT, ACTION,
%   SPEC) implements a portion of the climatic data of the MIL-STD-210C or
%   MIL-HDBK-310 worldwide air environment to 80 km (geometric or
%   approximately 262000 feet geometric) for absolute temperature,
%   pressure, density, and speed of sound for the input geopotential
%   altitude. 
%
%   Inputs for ATMOSNONSTD is:
%   H       :an array of M geopotential height in meters. 
%   ATYPE   :a string selecting the representation of 'Profile' or
%           'Envelope' for the atmospheric data. 
%           'Profile' is the realistic atmospheric profiles associated with
%           extremes at specified altitudes. 'Profile' is recommended for
%           simulation of vehicles vertically traversing the atmosphere or
%           when the total influence of the atmosphere is needed. 
%           'Envelope' uses extreme atmospheric values at each altitude.
%           'Envelope' is recommended for vehicles only horizontally
%           traversing the atmosphere without much change in altitude.  
%   EXTREME :a string selecting the atmospheric parameter which is the
%           extreme value. Atmospheric parameters that can be specified
%           are 'High temperature', 'Low temperature', 'High density', 'Low
%           density', 'High pressure' and 'Low pressure'. 'High pressure'
%           and 'Low pressure' are available only when ATYPE is 'Envelope'. 
%   FREQ    :a string selecting percent of time the extreme values would
%           occur. Valid values for FREQ include 'Extreme values', '1%',
%           '5%', '10%', and '20%'. 'Extreme values', '5%', and '20%' are
%           available only when ATYPE is 'Envelope'.  When using ATYPE of
%           'Envelope' and FREQ of '5%', '10%' and '20%' only the EXTREME
%           parameter selected (temperature, density, or pressure) has a
%           valid output.  All other parameter outputs will be zero. 
%   EXTALT  :a scalar value in kilometers selecting geometric altitude at
%           which the extreme values occur. EXTALT applies only when ATYPE
%           is 'Profile'.  Valid values for EXTALT include 5 (16404 ft), 10
%           (32808 ft), 20 (65617 ft), 30 (98425 ft), and 40 (131234 ft).
%   ACTION  :a string to determine action for out of range input. Specify if
%           out of range input invokes a 'Warning', 'Error', or no action
%           ('None'). The default is 'Warning'.
%   SPEC    :a string specifying the atmosphere model MIL-STD-210C or
%           MIL-HDBK-310: '210c' or '310'.  The default is '310'.  
%   
%   Outputs calculated for the non-standard atmosphere are: 
%   T       :an array of M temperature in kelvin.
%   a       :an array of M speed of sound in meters per second.
%   P       :an array of M air pressure in pascal.
%   rho     :an array of M air density in kilograms per meter cubed.
%
%   Limitations:
%
%   This function uses the metric data from MIL-STD-210 and MIL-HDBK-310
%   and has the limitations of MIL-STD-210 and MIL-HDBK-310. For more
%   information see the documentation.
%
%   Examples:
%
%   Calculate the non-standard atmosphere profile with high density
%   occurring 1% of the time at 5 kilometers from MIL-HDBK-310 at 1000
%   meters with warnings for out of range inputs: 
%      [T, a, P, rho] = atmosnonstd( 1000,'Profile','High density','1%',5 )
%
%   Calculate the non-standard atmosphere envelope with high pressure
%   occurring 20% of the time from MIL-STD-210C at 1000, 11000 and 20000
%   meters with errors for out of range inputs: 
%      [T, a, P, rho] = atmosnonstd([1000 11000 20000],'Envelope', ...
%                                    'High pressure','20%','Error','210c' )
%
%   See also ATMOSCIRA, ATMOSCOESA, ATMOSISA, ATMOSLAPSE.

%   Copyright 2000-2016 The MathWorks, Inc.

%   Limitation: All values are held below the geometric altitude of 0 m (0
%   feet) and above the geometric altitude of 80000 meters (approximately
%   262000 feet). The envelope atmospheric model has a few exceptions where
%   values are held below the geometric altitude of 1 kilometer
%   (approximately 3281 feet) and above the geometric altitude of 30000
%   meters (approximately 98425 feet). These exceptions are due to lack of
%   data in MIL-STD-210 or MIL-HDBK-310 for these conditions. 
%
%   In general, temperature values are interpolated linearly and density
%   values are interpolated logarithmically. Pressure and speed of sound
%   are calculated using a perfect gas relationship. The envelope
%   atmospheric model has a few exceptions where the extreme value is the 
%   only value provided as an output.  Pressure in these cases is 
%   interpolated logarithmically. These envelope atmospheric model 
%   exceptions apply to all cases of high and low pressure, high and low 
%   temperature, and high and low density, excluding the extreme values 
%   and 1% frequency of occurrence. These exceptions are due to lack of
%   data in MIL-STD-210 or MIL-HDBK-310 for these conditions.          
%
%   A limitation is that climatic data for the region south of 60 deg S
%   latitude is excluded from consideration in MIL-STD-210 or MIL-HDBK-310.
%
%   This function uses the metric version of data from the MIL-STD-210 or
%   MIL-HDBK-310 specifications. A limitation of this is some inconsistent
%   data between the metric and English data. Locations where these
%   inconsistencies occur are within the envelope data for low density, low
%   temperature, high temperature, low pressure, and high pressure.  The
%   most noticeable differences occur in the following values:
%  
%      For low density envelope data with 5% frequency, the density values 
%      in metric units are inconsistent at 4km and 18km and the density 
%      values in English units are inconsistent at 14km.
% 
%      For low density envelope data with 10% frequency, the density values 
%      in metric units are inconsistent at 18km and the density values in
%      English units are inconsistent at 14km.  
% 
%      For low density envelope data with 20% frequency, the density values
%      in English units are inconsistent at 14km.   
%
%      For low temperature envelope data with 20% frequency, the
%      temperature values at 20km are inconsistent.
%
%      For high pressure envelope data with 10% frequency, the pressure
%      values at 8km are inconsistent.
%
%   References:  
%   [1] Global Climatic Data for Developing Military Products
%   (MIL-STD-210C), 9 January 1987, Department of Defense, Washington, D.C.
%   [2] Global Climatic Data for Developing Military Products
%   (MIL-HDBK-310), 23 June 1997, Department of Defense, Washington, D.C.
%

narginchk(4, 7);

% Format inputs
atype = char(atype);
extreme = char(extreme);
freq = char(freq);

checkinputs();

action = 'warning';
spec = '310';
actionset = false;
specset = false;

switch nargin
    case 4
        % Envelope calculation no spec and no action
        % no processing required
    case 5
        switch lower(atype)
            case 'envelope'
                % Envelope calculation with spec or action
                if ischar( varargin{1} ) || isstring( varargin{1} )
                    specoraction(varargin{1});
                 else
                    error(message('aero:atmosnonstd:unknownOptionString'));
                end
            case 'profile'
                % Profile no spec and no action
                if isnumeric( varargin{1})
                    extalt = varargin{1};
                else
                    error(message('aero:atmosnonstd:unknownOptionNumeric'));
                end
            otherwise
                error(message('aero:atmosnonstd:unknownType', atype));
        end
    case 6
        switch lower(atype)
            case 'envelope'
                % Envelope calculation with spec and action
                if ((ischar(varargin{1}) || isstring(varargin{1})) && ...
                        (ischar( varargin{2} ) || isstring( varargin{2} )))
                    for i = [1 2]
                        specoraction(varargin{i});
                    end
                else
                    error(message('aero:atmosnonstd:unknownOptionStrings'));
                end
            case 'profile'
                % Profile without a spec or an action
                if (isnumeric( varargin{1}) && ...
                        (ischar( varargin{2} ) || isstring( varargin{2} )))
                    extalt = varargin{1};
                    specoraction(varargin{2});
                else
                    error(message('aero:atmosnonstd:unknownOptionNumberString'));
                end
            otherwise
                error(message('aero:atmosnonstd:unknownType', atype));
        end
    case 7
        % Profile calculation
        if (isnumeric( varargin{1})&& (ischar( varargin{2} ) || isstring( varargin{2} )) && ...
            (ischar( varargin{3} ) || isstring( varargin{3} )))
            if strcmpi(atype,'profile')
                extalt = varargin{1};
                for i = [2 3]
                    specoraction(varargin{i});
                end
            else
                error(message('aero:atmosnonstd:envelopeTooManyInputs'));
            end
        else
            error(message('aero:atmosnonstd:unknownOptionNumberStrings'));
        end
end


% volumetric mean radius of the earth in meters
R_EARTH = 6371010;

%  convert from geopotential altitude to geometric altitude
h_gm = (h*R_EARTH)./(R_EARTH - h);

% Enumerate the type
type = {'profile' 'envelope'};           
eatype = find(strcmpi(type,atype));

% Enumerate the frequency and the extreme
if eatype == 1
    type = {'1%' '10%'};
    etype = {'high temperature' 'low temperature' 'high density' 'low density'};           
else
    type = {'extreme values' '1%' '5%' '10%' '20%'};
    etype = {'high temperature' 'low temperature' 'high density' 'low density' ...
    'high pressure' 'low pressure'};           
    % define a dummy value for extalt
    extalt = 5;
end
efreq = find(strcmpi(type,freq));
eextreme = find(strcmpi(etype,extreme));
typep = [5 10 20 30 40];
eextalt = find(typep == extalt);

% Enumerate the action
type = {'none' 'warning' 'error'};
eaction = find(strcmpi(type,action));

% Enumerate the spec
type = {'' '310' '210c'};           
espec = find(strcmpi(type,spec));

% check for valid enumeration values for atype, extreme, freq and extalt
if ~(isempty(eatype) || isempty(eextreme) || isempty(efreq) || ...
                                                          isempty(eextalt))
    [T, a, P, rho] = nonstdmethods(h_gm, eatype, eextreme, efreq, eextalt,...
                                                          espec, eaction );
else
    index = cellfun('isempty',{eatype eextreme efreq eextalt});
    str = {atype extreme [freq '%'] num2str(extalt)};
    error(message('aero:atmosnonstd:badEnum', str{ index }));
end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkinputs()
        if ~isnumeric( h )
            % Altitude should be a numeric array.  Otherwise error.
            error(message('aero:atmosnonstd:notNumeric'));
        end
        if ~ischar( atype ) && ~isstring( atype )
            error(message('aero:atmosnonstd:atmosType'));
        end
        if ~ischar( extreme ) && ~isstring( extreme )
            error(message('aero:atmosnonstd:extreme'));
        end
        if ~ischar( freq ) && ~isstring( freq )
            error(message('aero:atmosnonstd:frequency'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function specoraction(str)
    str = char(str);
        switch lower(str)
            case {'210c' '310'}
                if specset
                    error(message('aero:atmosnonstd:specSet', spec));
                else
                    spec = lower(str);
                    specset = true;
                end
            case {'none' 'warning' 'error'}
                if actionset
                    error(message('aero:atmosnonstd:actionSet', action));
                else
                    action = lower(str);
                    actionset = true;
                end
            otherwise
                error(message('aero:atmosnonstd:unknownString', str));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

end
