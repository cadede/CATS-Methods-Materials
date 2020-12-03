function [T, rho] = atmosnrlmsise00(h,lat,lon,year,doy,sec,varargin)
%  ATMOSNRLMSISE00 Use NRLMSISE-00 atmosphere model.
%   [T, RHO] = ATMOSNRLMSISE00( H, LAT, LON, YEAR, DOY, SEC, LST, F107A,
%   F107, APH, FLAGS, ITYPE, OTYPE, ACTION ) implements the mathematical
%   representation of the 2001 United States Naval Research Laboratory Mass
%   Spectrometer and Incoherent Scatter Radar Exosphere (NRLMSISE-00),
%   of the MSIS(R) class model. NRLMSISE-00 calculates the neutral
%   atmosphere empirical model from the surface to lower exosphere (0 to
%   1,000,000 meters) with the option of including contributions from
%   anomalous oxygen which can affect satellite drag above 500,000 meters. 
%
%   Inputs for ATMOSNRLMSISE00 are:
%   H      :an array of M altitude in meters.
%   LAT    :an array of M geodetic latitude in degrees.
%   LON    :an array of M geodetic longitude in degrees.
%   YEAR   :an array of M year. Year is currently ignored in this model. 
%   DOY    :an array of M day of year. Day of year ranges from 1 to 365 (or
% 	       366).
%   SEC    :an array of M seconds in day in universal time (UT)
%   LST    :an array of M local apparent solar time (hours). To obtain a
%          physically realistic value, LST is set to (SEC/3600 + LON/15) by
%          default.  See Limitation section for more information.
%   F107A  :an array of M 81 day average of F10.7 flux (centered on doy).
%          If F107A is input, F107 and APH must also be input. The effects
%          of F107A are neither large nor well established below 80,000
%          meters, therefore the default value is set to 150. See
%          Limitation section for more information. 
%   F107   :an array of M daily F10.7 flux for previous day. If F107 is
%          input, F107A and APH must also be input. The effects of F107 are
%          neither large nor well established below 80,000 meters,
%          therefore the default value is set to 150. See Limitation
%          section for more information.
%   APH    :an array of M-by-7 of magnetic index information. If APH is
%          input, F107A and F107 must also be input. This information
%          consists of daily magnetic index (AP), 3 hour AP for current
%          time, 3 hour AP for 3 hours before current time, 3 hour AP for 6
%          hours before current time, 3 hour AP for 9 hours before current
%          time, average of eight 3 hour AP indices from 12 to 33 hours
%          prior to current time, and average of eight 3 hour AP indices
%          from 36 to 57 hours prior to current time. The effects of daily
%          magnetic index are neither large nor well established below
%          80,000 meters, therefore the default value is set to 4. See
%          Limitation section for more information.  
%   FLAGS  :a numerical array of 23 values for setting particular
%          variations in calculation the output.  Setting a value to 0.0
%          removes that value's effect on the output.  Setting a value to
%          1.0 applies the main and the cross term effects of that value
%          on the output.  Setting a value to 2.0 applies only the cross
%          term effect of that value on the output.  Additionally setting
%          FLAGS(9) = -1 uses the entire matrix APH rather than just
%          APH(:,1). The variations contained in FLAGS are ordered as
%          follows: 
%           FLAGS(1)  :F10.7 effect on mean  
%           FLAGS(2)  :Time independent
%           FLAGS(3)  :Symmetrical annual    
%           FLAGS(4)  :Symmetrical semi-annual
%           FLAGS(5)  :Asymmetrical annual   
%           FLAGS(6)  :Asymmetrical semi-annual
%           FLAGS(7)  :Diurnal               
%           FLAGS(8)  :Semi-diurnal
%           FLAGS(9)  :Daily AP             
%           FLAGS(10) :All UT, longitudinal effects
%           FLAGS(11) :Longitudinal        
%           FLAGS(12) :UT and mixed UT, longitudinal
%           FLAGS(13) :Mixed AP, UT, longitudinal     
%           FLAGS(14) :Ter-diurnal
%           FLAGS(15) :Departures from diffusive equilibrium
%           FLAGS(16) :All exospheric temperature variations         
%           FLAGS(17) :All variations from 120,000 meter temperature (TLB)
%           FLAGS(18) :All lower thermosphere (TN1) temperature variations           
%           FLAGS(19) :All 120,000 meter gradient (S) variations
%           FLAGS(20) :All upper stratosphere (TN2) temperature variations           
%           FLAGS(21) :All variations from 120,000 meter values (ZLB)
%           FLAGS(22) :All lower mesosphere temperature (TN3) variations           
%           FLAGS(23) :Turbopause scale height variations
%          The default values are 1.0 for all FLAGS.
%   OTYPE  :a string specifying if the total mass density output will
%          include anomalous oxygen ('Oxygen') or not ('NoOxygen'). The
%          default is 'NoOxygen'.
%   ACTION :a string to determine action for out of range input. Specify if
%          out of range input invokes a 'Warning', 'Error', or no action
%          ('None'). The default is 'Warning'.
%
%   Outputs calculated for the NRLMSISE-00 model are: 
%   T      :an array of M-by-2 values of temperatures.  These values are
%           exospheric temperature in Kelvin and temperature at altitude in
%           Kelvin.
%   RHO    :an array of M-by-9 values of densities.  These values are
%           HE number density in meters^-3, O number density in meters^-3,
%           N2 number density in meters^-3, O2 number density in meters^-3,
%           AR number density in meters^-3, total mass density in kilogram 
%           per meters cubed, H number density in meters^-3, N number
%           density in meters^-3, and Anomalous oxygen number density in
%           meters^-3. 
%
%   Limitation:
%
%   If array length, M, is 23 and all available inputs are not specified,
%   FLAGS will always be assumed to be set.
%
%   This function has the limitations of the NRLMSISE-00 model. For more
%   information see the documentation. 
%  
%   SEC, LST, and LON are used independently in the NRLMSISE-00 model and
%   are not of equal importance for every situation. For the most
%   physically realistic calculation these three variables are chosen to be
%   consistent by default (LST = SEC/3600 + LON/15). Departures from the
%   prior formula for LST can be included if available but are of minor
%   importance.
% 
%   F107 and F107A values used to generate the model correspond to the 10.7
%   cm radio flux at the actual distance of the Earth from the Sun rather
%   than the radio flux at 1 AU. The following site provides both classes
%   of values: ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
%
%   Examples:
%
%   Calculate the temperatures, densities not including anomalous oxygen
%   using NRLMSISE-00 model at 10000 meters, 45 degrees latitude, -50
%   degrees longitude, on January 4, 2007 at 0 UT using default values for
%   flux, magnetic index data, and local solar time with out of range
%   actions generating warnings:    
%      [T, rho] = atmosnrlmsise00( 10000, 45, -50, 2007, 4, 0)
%
%   Calculate the temperatures, densities not including anomalous oxygen
%   using NRLMSISE-00 model at 10000 meters, 45 degrees latitude, -50
%   degrees longitude, and at 25000 meters, 47 degrees latitude, -55
%   degrees longitude on January 4, 2007 at 0 UT using default values for
%   flux, magnetic index data, and local solar time with out of range
%   actions generating warnings:    
%      [T, rho] = atmosnrlmsise00( [10000; 25000], [45; 47], [-50; -55], [2007; 2007], [4; 4], [0; 0] )
%
%   Calculate the temperatures, densities including anomalous oxygen
%   using NRLMSISE-00 model at 10000 meters, 45 degrees latitude, -50
%   degrees longitude, on January 4, 2007 at 0 UT using default values for
%   flux, magnetic index data, and local solar time with out of range
%   actions generating errors:    
%      [T, rho] = atmosnrlmsise00( 10000, 45, -50, 2007, 4, 0, 'Oxygen', 'Error' )
%
%   Calculate the temperatures, densities including anomalous oxygen
%   using NRLMSISE-00 model at 100000 meters, 45 degrees latitude, -50
%   degrees longitude, on January 4, 2007 at 0 UT using defined values for
%   flux, and magnetic index data, and default local solar time with out of
%   range actions generating no message:    
%      aph = [17.375 15 20 15 27 (32+22+15+22+9+18+12+15)/8 (39+27+9+32+39+9+7+12)/8]
%      f107 = 87.7
%      nov_6days  = [78.6 78.2 82.4 85.5 85.0 84.1]
%      dec_31daymean = 84.5
%      jan_31daymean = 83.5
%      feb_13days = [89.9 90.3 87.3 83.7 83.0 81.9 82.0 78.4 76.7 75.9 74.7 73.6 72.7]
%      f107a = (sum(nov_6days) + sum(feb_13days) + (dec_31daymean + jan_31daymean)*31)/81
%      flags = ones(1,23)
%      flags(9) = -1
%      [T, rho] = atmosnrlmsise00( 100000, 45, -50, 2007, 4, 0, f107a, f107, aph, flags, 'Oxygen', 'None' )
%
%   See also ATMOSCIRA, ATMOSCOESA.

%   Copyright 2007-2018 The MathWorks, Inc.

%   Limitation:
%
%   Low-order spherical harmonics are used to describe the major variations
%   through out the atmosphere including latitude, annual, semiannual, and
%   simplified local time and longitude variations. 
%
%   Number densities for O, H, and N are set to zero below 72,500 meters.
%  
%   Exospheric temperature is set to global average for altitudes below
%   120,000 meters. The 120,000 meter gradient is left at global average
%   value for altitudes below 72,500 meters. 

%   References:  
%   [1] http://modelweb.gsfc.nasa.gov/atmos/nrlmsise00.html
%   [2] http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
%   [3] http://www.brodo.de/english/pub/nrlmsise/index.html

narginchk(6, 13);

otypestr = {'nooxygen' 'oxygen'}; 
actionstr = {'none' 'warning' 'error'};

checkinputs();

% wrap latitude and longitude if needed
[lat_wrapped, lat, lon] = wraplatitude( lat, lon, 'deg' );

% check and fix angle wrapping in longitude
[lon_wrapped, lon] = wraplongitude( lon, 'deg', '180' );

% calculate local solar time
lst = sec(1)/3600 + lon(1)/15;

% default values
f107a = 150;
f107 = 150;
aph = zeros(1,7);
aph(:,1) = 4;
flags = ones(1,24);
otype = 'nooxygen';
action = 'warning';
actionset = false;
otypeset = false;
f107set = false;

if nargin == 7
    checknargin7(varargin{:})
elseif nargin == 8
    checknargin8(varargin{:})
elseif nargin == 9
    checknargin9(varargin{:})
elseif nargin == 10
    checknargin10(varargin{:})
elseif nargin == 11
    checknargin11(varargin{:})
elseif nargin == 12
    checknargin12(varargin{:})
elseif nargin == 13
    checknargin13(varargin{:})
end

% enumerate string inputs
otypeidx = find(strcmpi(otypestr, otype));
actionidx = find(strcmpi(actionstr, action));

% check lengths of inputs against altitude
checkinputlength();

% check doy for valid range
checkdoy();

h      = convlength( h , 'm' , 'km' );
idoy   = int32(doy);
iyear  = int32(year);
iflags = int32(flags);

% check altitude for valid range
checkalt();

% check if f107a, f107, and aph should have been set
checkf107af107aph();

% check latitude and longitude for valid range
checklatlon();

[T,rho] = nrlmsisemethod(idoy(:), iyear(:), sec(:), lat(:), ...
                         lon(:), h(:), lst(:), f107a(:), f107(:), ...
                         aph(:,1), aph, iflags(:));

if (otypeidx == 2) 
    % include anomalous oxygen in total mass density
	 rho(:,6) = 1.66E-27 * (16.0 *  rho(:,9)) + rho(:,6);
end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkinputs( )
        if ~isnumeric( h )
            % altitude should be a numeric array.  Otherwise error.
            error(message('aero:atmosnrlmsise00:notNumericAltitude'));
        end
        if ~isnumeric( lat )
            % latitude should be a numeric array.  Otherwise error.
            error(message('aero:atmosnrlmsise00:notNumericLatitude'));
        end
        if ~isnumeric( lon )
            % longitude should be a numeric array.  Otherwise error.
            error(message('aero:atmosnrlmsise00:notNumericLongitude'));
        end
        if ~isnumeric( year )
            % year should be a numeric array.  Otherwise error.
            error(message('aero:atmosnrlmsise00:notNumericYear'));
        end
        if ~isnumeric( doy )
            % doy should be a numeric array.  Otherwise error.
            error(message('aero:atmosnrlmsise00:notNumericDOY'));
        end
        if ~isnumeric( sec )
            % sec should be a numeric array.  Otherwise error.
            error(message('aero:atmosnrlmsise00:notNumericSec'));
        end

    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkdoy( )
        if any(((doy > 365) & ~leapyear(year))|((doy > 366) & leapyear(year)) | ( doy < 1 ))
            
            % Create function handle array to handle messages based on action
            fh = {@() disp(''), ...
                @() warning(message('aero:atmosnrlmsise00:invalidDOY')), ...
                @() error(message('aero:atmosnrlmsise00:invalidDOY'))};
            % Call appropriate function
            fh{actionidx}()            
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkalt( )
        if any(( h > 1000 ) | ( h < 0 ))
            
            % Create function handle array to handle messages based on action
            fh = {@() disp(''), ...
                @() warning(message('aero:atmosnrlmsise00:invalidAltitude')), ...
                @() error(message('aero:atmosnrlmsise00:invalidAltitude'))};

            % Call appropriate function
            fh{actionidx}()            
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checklatlon()

        if lat_wrapped
            
            % Create function handle array to handle messages based on action
            fh = {@() disp(''), ...
                  @() warning(message('aero:atmosnrlmsise00:latitudeWrap')), ...
                  @() error(message('aero:atmosnrlmsise00:latitude90'))};

            % Call appropriate function
            fh{actionidx}()            
        end
        if lon_wrapped            
            
            % Create function handle array to handle messages based on action
            fh = {@() disp(''), ...
                  @() warning(message('aero:atmosnrlmsise00:longitudeWrap')), ...
                  @() error(message('aero:atmosnrlmsise00:longitude180'))};

            % Call appropriate function
            fh{actionidx}()                        
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function otypeoraction(str)
        str = char(str);
        switch lower(str)
            case otypestr
                if otypeset
                    error(message('aero:atmosnrlmsise00:otypeSet',otype));
                else
                    otype = lower(str);
                    otypeset = true;
                end
            case {'none' 'warning' 'error'}
                if actionset
                    error(message('aero:atmosnrlmsise00:actionSet',action));
                else
                    action = lower(str);
                    actionset = true;
                end
            otherwise
                error(message('aero:atmosnrlmsise00:unknownString',str));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function setf107af107aph(a, b, c)
        f107a = varargin{a};
        f107  = varargin{b};
        aph   = varargin{c};
        f107set = true;
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function lstorflag(input)
        if (numel(input) == 23)
            flags(2:24) = input;
        else
            lst = input;
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function setflagsandlst(input,input_lst)
        if (numel(input) == 23)
            flags(2:24) = input;
            lst = input_lst;
        else
            error(message('aero:atmosnrlmsise00:flagLength'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkf107af107aph()
        if (any(h > 80) & ~f107set) %#ok<AND2>

            % Create function handle array to handle messages based on action
            fh = {@() disp(''), ...
                @() warning(message('aero:atmosnrlmsise00:setf107af107aph')), ...
                @() error(message('aero:atmosnrlmsise00:setf107af107aph'))};

            % Call appropriate function
            fh{actionidx}()            
         end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkinputlength()
        lengths = [numel(h) numel(lat) numel(lon) numel(year) numel(doy) ...
            numel(sec) numel(lst) numel(f107a) numel(f107) size(aph,1)];
        len = max(lengths);
        lengthstr = {'H' 'LAT' 'LON' 'YEAR' 'DOY' 'SEC' 'LST' 'F107A' 'F107' 'APH'};

        badlength = find(lengths ~= len & lengths ~= 1);
        if ~isempty(badlength)
            error(message('aero:atmosnrlmsise00:wrongInputLength', ...
                lengthstr{badlength(1)}, lengthstr{find(lengths == max(lengths),1)}));
        end
        if (len ~= 1)
            % set all inputs of length 1 to max length.
            idx1 = find(lengths == 1);
            for j = idx1
                eval([lower(lengthstr{j}) '= ' ...
                      'ones(len,1)*' lower(lengthstr{j}) ';']);
            end
            % set value for lst if not entered by user
            if (lst(1) == sec(1)/3600 + lon(1)/15)
                lst = sec(:)/3600 + lon(:)/15;
            end
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checknargin7(varargin)
        if isnumeric(varargin{1})
            lstorflag( varargin{1} );
        else
            otypeoraction( varargin{1} );
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checknargin8(varargin)
        if isnumeric(varargin{1})
            if ischar(varargin{2}) || isstring(varargin{2})
                otypeoraction( varargin{2} );
                lstorflag( varargin{1} );
            else
                setflagsandlst(varargin{2},varargin{1})
            end
        else
            if ischar(varargin{2}) || isstring(varargin{2})
                for j = [1 2]
                    otypeoraction( varargin{j} );
                end
            else
                error(message('aero:atmosnrlmsise00:expect2Strings'));
            end
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checknargin9(varargin)
        if isnumeric(varargin{1})
            if ((ischar(varargin{2}) || isstring(varargin{2})) && ...
                    (ischar(varargin{3})|| isstring(varargin{3})))
                for j = [2 3]
                    otypeoraction( varargin{j} );
                end
                lstorflag( varargin{1} );
            elseif (isnumeric(varargin{2}) && isnumeric(varargin{3}))
                setf107af107aph(1, 2, 3);
            elseif ((numel(varargin{2}) == 23) && (ischar(varargin{3}) || isstring(varargin{3})))
                lst = varargin{1};
                flags(2:24) = varargin{2};
                otypeoraction( varargin{3} );
            else
                error(message('aero:atmosnrlmsise00:incorrectInput89'));
            end
        else
            error(message('aero:atmosnrlmsise00:incorrectInput7', nargin+6 ));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checknargin10(varargin)
        if (isnumeric(varargin{1}) && isnumeric(varargin{2}))
            if ((ischar(varargin{3}) || isstring(varargin{3})) && ...
                    (ischar(varargin{4})|| isstring(varargin{4})))
                for j = [3 4]
                    otypeoraction( varargin{j} );
                end
                setflagsandlst(varargin{2},varargin{1})
            elseif (isnumeric(varargin{3}) && isnumeric(varargin{4}))
                if ((size(varargin{3},2) == 7) && (numel(varargin{4}) == 23))
                    setf107af107aph(1, 2, 3);
                    flags(2:24) = varargin{4};
                elseif (size(varargin{4},2) == 7)
                    lst   = varargin{1};
                    setf107af107aph(2, 3, 4);
                else
                    error(message('aero:atmosnrlmsise00:incorrectInput9or10'));
                end
            elseif ((size(varargin{3},2) == 7) && (isnumeric(varargin{3}) && ...
                    (ischar(varargin{4})|| isstring(varargin{4}))))
                setf107af107aph(1, 2, 3);
                otypeoraction( varargin{4} );
            else
                error(message('aero:atmosnrlmsise00:incorrectInput910'));
            end
        else
            error(message('aero:atmosnrlmsise00:incorrectInput7or8', nargin+6 ));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checknargin11(varargin)
        if (isnumeric(varargin{1}) && ...
                isnumeric(varargin{2}) && isnumeric(varargin{3}))
            if ((size(varargin{3},2) == 7) && (ischar(varargin{4})|| isstring(varargin{4})) && ...
                    (ischar(varargin{5})|| isstring(varargin{5})))
                setf107af107aph(1, 2, 3);
                for j = [4 5]
                    otypeoraction( varargin{j} );
                end
            elseif (isnumeric(varargin{4}) && isnumeric(varargin{5}) && ...
                    (size(varargin{4},2) == 7) && (numel(varargin{5}) == 23))
                setflagsandlst(varargin{5},varargin{1})
                setf107af107aph(2, 3, 4);
            elseif (isnumeric(varargin{4}) && (ischar(varargin{5})|| isstring(varargin{5})))
                otypeoraction( varargin{5} );
                if (size(varargin{3},2) == 7) && (numel(varargin{4}) == 23)
                    setf107af107aph(1, 2, 3);
                    flags(2:24) = varargin{4};
                elseif (size(varargin{4},2) == 7)
                    lst    = varargin{1};
                    setf107af107aph(2, 3, 4);
                else
                    error(message('aero:atmosnrlmsise00:incorrectInput9or10'));
                end
            else
                error(message('aero:atmosnrlmsise00:incorrectInput1011'));
            end
        else
            error(message('aero:atmosnrlmsise00:incorrectInput78or9', nargin+6 ));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checknargin12(varargin)
        if (isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
                isnumeric(varargin{3}) && isnumeric(varargin{4}) && ...
                (ischar(varargin{6}) || isstring(varargin{6})))
            if ((size(varargin{4},2) == 7) && isnumeric(varargin{5}) && (numel(varargin{5}) == 23))
                lst    = varargin{1};
                setf107af107aph(2, 3, 4);
                flags(2:24) = varargin{5};
                otypeoraction( varargin{6} );
            else
                if (ischar(varargin{5}) || isstring(varargin{5}))
                    if (size(varargin{3},2) == 7) && (numel(varargin{4}) == 23)
                        setf107af107aph(1, 2, 3);
                        flags(2:24) = varargin{4};
                    elseif (size(varargin{4},2) == 7) 
                        lst    = varargin{1};
                        setf107af107aph(2, 3, 4);
                    else
                        error(message('aero:atmosnrlmsise00:incorrectInput9or10'));
                    end
                    for j = [5 6]
                        otypeoraction( varargin{j} );
                    end
                end
            end
        else
            error(message('aero:atmosnrlmsise00:incorrectInput78910or12', nargin+6 ));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   function checknargin13(varargin)
        if (isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
                isnumeric(varargin{3}) && ...
                (isnumeric(varargin{4}) && (size(varargin{4},2) == 7)) && ...
                (isnumeric(varargin{5}) && (numel(varargin{5}) == 23)) && ...
                (ischar(varargin{6}) || isstring(varargin{6})) && ...
                (ischar(varargin{7}) || isstring(varargin{7})))
            lst         = varargin{1};
            setf107af107aph(2, 3, 4);
            flags(2:24) = varargin{5};
            for j = [6 7]
                otypeoraction( varargin{j} );
            end
        else
            error(message('aero:atmosnrlmsise00:incorrectInputOrder'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end
