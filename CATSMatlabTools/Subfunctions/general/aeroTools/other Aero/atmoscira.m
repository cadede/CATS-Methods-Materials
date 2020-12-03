function [T, alt, ZW] = atmoscira( lat, ctype, coord, varargin )
%   ATMOSCIRA Use COSPAR International Reference Atmosphere 1986.
%   [T, ALT, ZW] = ATMOSCIRA( LAT, CTYPE, COORD ) implements
%   the mathematical representation of the Committee on Space Research
%   (COSPAR) International Reference Atmosphere (CIRA) from 1986.  The CIRA
%   1986 model provides a mean climatology of temperature, zonal wind, and
%   geopotential height or pressure with nearly pole-to-pole coverage
%   (80 degrees south to 80 degrees north) for 0 to 120 kilometers, 
%   encompassing the troposphere, middle atmosphere, and lower thermosphere. 
%   It can be used as a function of pressure or geopotential height. 
%
%   Inputs for ATMOSCIRA are:
%   LAT    :an array of M geodetic latitudes in degrees where north
%          latitude is positive, and south latitude is negative.
%   CTYPE  :a string selecting the representation of 'Pressure' or
%           'GPHeight' for the coordinate type to be used. 
%           'Pressure' uses pressure in pascal. 
%           'GPHeight' uses geopotential height in meters.
%   COORD  :either an array of M pressures in pascal or an array of M
%           geopotential height in meters depending on CTYPE selected.  
%   MTYPE  :a string selecting mean value type of 'Monthly' or 'Annual'.
%           'Annual' only valid when CTYPE is 'Pressure'. The default
%           value is 'Monthly'.
%   MONTH  :a scalar value selecting the month at which the mean values are
%           taken. MONTH applies only when MTYPE is 'Monthly'.  Valid
%           values for Month include 1 through 12 (January through
%           December). The default value is 1 (January).
%   ACTION :a string to determine action for out of range input. Specify if
%          out of range input invokes a 'Warning', 'Error', or no action
%          ('None'). The default is 'Warning'.
%
%   Outputs calculated for the CIRA 1986 model are: 
%   T      :an array of M temperatures in kelvin for 'Monthly' MTYPE, or
%           an array of M-by-7 values for temperatures for 'Annual' MTYPE.
%           These values are annual mean temperature in kelvin, annual
%           temperature cycle amplitude in kelvin, annual temperature cycle
%           phase in month of maximum, semi-annual temperature cycle
%           amplitude in kelvin, semi-annual temperature cycle phase in
%           month of maximum, ter-annual temperature cycle amplitude in
%           kelvin, and ter-annual temperature cycle phase in month of
%           maximum. 
%   ALT    :an array of M geopotential heights in meters or an array of M
%           air pressures in pascal depending on CTYPE for 'Monthly' MTYPE.
%           'Pressure' will output geopotential height and a CYPE of
%           A CTYPE of 'GPHeight' will output pressure.
%           an array of M-by-7 values for geopotential heights for 'Annual'
%           MTYPE. These values are annual mean geopotential heights in
%           meters, annual geopotential heights cycle amplitude in meters,
%           annual geopotential heights cycle phase in month of maximum,
%           semi-annual geopotential heights cycle amplitude in meters,
%           semi-annual geopotential heights cycle phase in month of
%           maximum, ter-annual geopotential heights cycle amplitude in
%           meters, and ter-annual geopotential heights cycle phase in
%           month of maximum.  This M-by-7 array is only defined for the
%           northern hemisphere (LAT > 0).
%   ZW     :an array of zonal winds in meters per second for 'Monthly'
%           MTYPE.  
%           an array of M-by-7 values for zonal winds for 'Annual' MTYPE.
%           These values are annual mean zonal winds in meters per second,
%           annual zonal winds cycle amplitude in meters per second, annual
%           zonal winds cycle phase in month of maximum, semi-annual zonal
%           winds cycle amplitude in meters per second, semi-annual zonal
%           winds cycle phase in month of maximum, ter-annual zonal winds
%           cycle amplitude in meters per second, and ter-annual zonal
%           winds cycle phase in month of maximum. 
%
%   Limitation:
%
%   This function uses a corrected version of the CIRA data files as
%   provided by J. Barnett in July 1990 in ASCII format.  
%
%   This function has the limitations of the CIRA 1986 model. The values
%   for the CIRA 1986 model are limited to the regions of 80 degrees south 
%   to 80 degrees north on the Earth and geopotential heights of 0 to 120 
%   kilometers. In each monthly mean data set, values at 80 degrees south 
%   for 101300 pascal or 0 meters were omitted since these levels are 
%   within the Antarctic land mass.
%   For zonal mean pressure in constant altitude coordinates, pressure data 
%   is not available below 20 kilometers; therefore this is the bottom level
%   of the CIRA climatology. 
%
%   Examples:
%
%   Calculate the mean monthly values for temperature, geopotential height
%   and zonal wind using CIRA 1986 model at 45 degrees latitude and 101300
%   pascal for January with out of range actions generating warnings:
%      [T, alt, zwind] = atmoscira( 45, 'Pressure', 101300 )
%
%   Calculate the mean monthly values for temperature, pressure, and zonal
%   wind using CIRA 1986 model at 45 degrees latitude and 20000 meters for
%   October with out of range actions generating warnings:  
%      [T, pres, zwind] = atmoscira( 45, 'GPHeight', 20000, 'Monthly', 10 )
%
%   Calculate the mean monthly values for temperature, pressure, and zonal
%   wind using CIRA 1986 model at 45 and -30 degrees latitude and 20000
%   meters for October with out of range actions generating errors:  
%      [T, pres, zwind] = atmoscira( [45 -30], 'GPHeight', 20000, 10, 'error' )
%
%   Calculate the mean monthly values for temperature, geopotential height,
%   and zonal wind using CIRA 1986 model at 45 degrees latitude and 2000
%   pascal and at -30 degrees latitude and 101300 pascal for September with
%   out of range actions generating warnings: 
%      [T, alt, zwind] = atmoscira( [45 -30], 'Pressure', [2000 101300], 9 )
%
%   Calculate the annual values for temperature, geopotential height,
%   and zonal wind using CIRA 1986 model at 45 degrees latitude and 2000
%   pascal with out of range actions generating warnings: 
%      [T, alt, zwind] = atmoscira( 45, 'Pressure', 2000, 'Annual' )
%
%   See also ATMOSCOESA, ATMOSISA, ATMOSLAPSE, ATMOSNONSTD, ATMOSNRLMSISE00,
%    ATMOSPALT.

%   Copyright 2006-2018 The MathWorks, Inc.

%   The CIRA data files contain the following data:
%
%   Temperature data as a function of geopotential height and latitude has
%   points for height of 0 to 120 kilometers at 5 kilometer intervals and
%   for latitudes between +/- 80 degrees at 10 degree intervals.  
%
%   Pressure data as a function of geopotential height and latitude has
%   points for height of 20 to 120 kilometers at 5 kilometer intervals and
%   for latitudes between +/- 80 degrees at 10 degree intervals.  
%
%   Zonal wind data as a function of geopotential height and latitude has
%   points for height of 0 to 120 kilometers at 5 kilometer intervals and
%   for latitudes between +/- 80 degrees at 10 degree intervals.  
%
%   Temperature data as a function of log-pressure and latitude has
%   points for log-pressure of 0.0 to 17.5 sh at .5 sh intervals where sh =
%   -ln(P/P0), P0 = 101300 pascal, 0.0025 pascal <= P <= 101300 pascal and
%   for latitudes between +/- 80 degrees at 10 degree intervals.
%
%   Geopotential height data as a function of log-pressure and latitude has
%   points for log-pressure of 0.0 to 17.5 sh at .5 sh intervals where sh =
%   -ln(P/P0), P0 = 101300 pascal, 0.0025 pascal <= P <= 101300 pascal and
%   for latitudes between +/- 80 degrees at 10 degree intervals.
%
%   Zonal Wind data as a function of log-pressure and latitude has
%   points for log-pressure of 0.0 to 17.5 sh at .5 sh intervals where sh =
%   -ln(P/P0), P0 = 101300 pascal, 0.0025 pascal <= P <= 101300 pascal
%   and for latitudes between +/- 80 degrees at 10 degree intervals.
%
%   Additionally supplemental data to these cross sections is included. The
%   climatology of temperature and zonal wind in constant pressure
%   coordinates was Fourier analyzed at each latitude and pressure level
%   and the results were presented in the form of latitude height sections
%   of the annual mean and the amplitudes and phases of the annual,
%   semiannual, and terannual harmonics. 

%   References:  
%   [1] Fleming, E. L., Chandra, S., Shoeberl, M. R., Barnett, J. J.,
%   'Monthly Mean Global Climatology of Temperature, Wind, Geopotential 
%   Height and Pressure for 0-120 km', NASA TM100697, February 1988.
%   [2] http://modelweb.gsfc.nasa.gov/atmos/cospar1.html
%

narginchk(3, 6);

ctypestr = {'pressure' 'gpheight'};
mtypestr = {'monthly' 'annual'};
actionstr = {'none' 'warning' 'error'};

checkinputs();

mtype = 'monthly';
month = 1;
action = 'warning';
latalt = lat;
coordp = coord;

if nargin == 4
    if ~ischar( varargin{1} ) && ~isstring( varargin{1} )
        month = varargin{1};
        checkmonth();
    else
        if any(strcmpi(mtypestr, varargin{1}))
            checkmtype( varargin{1} );
        else
            evalaction( varargin{1} );
        end
    end
elseif nargin == 5
    if ~ischar( varargin{1} ) && ~isstring( varargin{1} )
        month = varargin{1};
        checkmonth();
        evalaction( varargin{2} );
    else
        checkmtype( varargin{1} );
        if (isnumeric( varargin{2} ) && strcmpi(mtype,'monthly'))
            month = varargin{2};
            checkmonth();
        else
            if ~ischar( varargin{2} ) && ~isstring( varargin{2} )
                error(message('aero:atmoscira:nonActionString'));
            end
            evalaction( varargin{2} );
        end
    end
elseif nargin == 6
    checkmtype(varargin{1});
    if (isnumeric( varargin{2} ) && strcmpi(mtype,'monthly'))
        month = varargin{2};
        checkmonth();
    else
        if strcmpi(mtype,'monthly')
            error(message('aero:atmoscira:unknownMonth'));
        else
            error(message('aero:atmoscira:annualTooManyInputs'));
        end
    end
    evalaction( varargin{3} );
end

% determine appropriate data set
ctypeidx = find(strcmpi(ctypestr, ctype));
mtypeidx = find(strcmpi(mtypestr, mtype));
actionidx = find(strcmpi(actionstr, action));

% check latitude and coordinate to make sure they are within the correct bounds
checklatitude();
checkcoordinate();

% data sets
file = {{'aerocirapres.mat' 'aerocirapresann.mat'} ...
        {'aerociraalt.mat' 'aerociraalt.mat'}};

interpolateCIRA();

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkinputs( )
        ctype = char(ctype);
        if ~isnumeric( lat )
            % latitude should be a numeric array.  Otherwise error.
            error(message('aero:atmoscira:notNumericLatitude'));
        end

        if ~ischar( ctype )
            % coordinate type should be a string.  Otherwise error.
            error(message('aero:atmoscira:notStringCoordType'));
        end

        if ~isnumeric( coord )
            % coordinate should be a numeric array.  Otherwise error.
            error(message('aero:atmoscira:notNumericCoord'));
        end
        if ~any(strcmpi(ctypestr,ctype))
            error(message('aero:atmoscira:unknownStringCtype', ctype));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkmonth( )
        if any( month < 1 | month > 12 )
            % clip to region
            month = min(max(month,1),12);
            warning(message('aero:atmoscira:invalidMonth'));
        end
        if any(length(num2str(month)) > 2)
            % create integer
            month = floor(month);
            month(month == 0) = 1;
            warning(message('aero:atmoscira:nonIntegerMonth'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function evalaction( input )
        input = char(input);
        if any(strcmpi(actionstr, input))
            action = lower( input );
        else
            error(message('aero:atmoscira:unknownStringAction', input));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkmtype( input )
        input = char(input);
        if any(strcmpi(mtypestr, input))
            mtype = lower( input );
            if (strcmpi(ctype,'gpheight') && strcmpi(mtype,'annual'))
                mtype = 'monthly';
                warning(message('aero:atmoscira:annualWHeight'));
            end
        else
            error(message('aero:atmoscira:unknownStringMtype', input));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checklatitude( )
        if any(abs(lat) > 80 )
            
            % clip to region
            lat = min(max(lat,-80),80);
            latalt = lat;
                        
            % Create function handle array to handle messages based on action
            fh = {@() disp(''), ...
                @() warning(message('aero:atmoscira:warnInvalidLatitude')), ...
                @() error(message('aero:atmoscira:invalidLatitude'))};
            % Call appropriate function
            fh{actionidx}()            
        end
        if (( mtypeidx == 2 ) && any( latalt < 0 ))
            % Annual altitude only defined for northern hemisphere
            latalt = abs(latalt);
                        
            % Create function handle array to handle messages based on action
            fh = {@() disp(''), ...
                @() warning(message('aero:atmoscira:warnInvalidLatForAlt')), ...
                @() error(message('aero:atmoscira:invalidLatForAlt'))};

            % Call appropriate function
            fh{actionidx}()            
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkcoordinate( )
        
        coordbptab = {{[17.50 0] [17.50 0]} {[120 0] [120 20]}};
        
        if (ctypeidx == 1)
            % get log-pressure coordinate
            coord = -log(coord./101300);
            coordp = coord;
        else
            % convert from meter to kilometers
            coord = convlength(coord,'m','km');
            coordp = convlength(coordp,'m','km');
        end

        % select correct coordinate range
        coordbp = coordbptab{ctypeidx}{1};
        coordbppres = coordbptab{ctypeidx}{2};

        if any( coord > coordbp(1) | coord < coordbp(end) )

            % clip to region
            coord = min(max(coord,coordbp(end)),coordbp(1));
            coordp = coord;
            
            invalidcoord(coordbp,' all ');
        end

        if (any( coordp < coordbppres(end) ) && (ctypeidx == 2))
            % Geopotential height only to 20 km for pressure tables
            coordp( coordp < coordbppres(end) ) = coordbppres(end);

            invalidcoord(coordbppres,' pressure ');
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function invalidcoord(coordbp, tabstr)

        coordunittab = {'log-pressure' 'kilometers'};
        coordtypetab = {'pressure' 'geopotential height'};
        coordunit = coordunittab{ctypeidx};
        coordtype = coordtypetab{ctypeidx};

        % Create function handle array to handle messages based on action
        fh = {@() disp(''), ...
            @() warning(message('aero:atmoscira:warnInvalidCoord', coordtype, sprintf( '%g', coordbp( end ) ), sprintf( '%g', coordbp( 1 ) ), coordunit, tabstr)), ...
            @() error(message('aero:atmoscira:invalidCoord', coordtype, sprintf( '%g', coordbp( end ) ), sprintf( '%g', coordbp( 1 ) ), coordunit))};

        % Call appropriate function
        fh{actionidx}()
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function interpolateCIRA()

        % Predefine variables for static workspace
        T_table=0;ZW_table=0;ZH_table=0;
        coordbp=0;coordshort_bp=0;latbp=0;nlatbp=0;
        
        % load data
        load(file{ctypeidx}{mtypeidx});

        if ( mtypeidx == 1 )
            coordbp_cell = {{coordbp coordbp} {coordbp coordshort_bp}};

            T_tab    = T_table(:,:,month);
            ZW_tab   = ZW_table(:,:,month);
            alt_tab  = ZH_table(:,:,month);
            coordbp_tab  = coordbp_cell{ctypeidx}{1};
            coordbpp_tab = coordbp_cell{ctypeidx}{2};

            T   = interp2(latbp,coordbp_tab,T_tab,lat,coord);
            ZW  = interp2(latbp,coordbp_tab,ZW_tab,lat,coord);
            alt = interp2(latbp,coordbpp_tab,alt_tab,latalt,coordp);
            if (ctypeidx == 2)
                % convert from log-pressure to pascal
                alt = 101300*exp(-alt);
            end
        elseif ( mtypeidx == 2 )
            for k = 7:-1:1
                T(:,k)   = interp2(latbp,coordbp,T_table(:,:,k),lat,coord); 
                ZW(:,k)  = interp2(latbp,coordbp,ZW_table(:,:,k),lat,coord);
                alt(:,k) = interp2(nlatbp,coordbp,ZH_table(:,:,k),latalt,coord);
            end
        else
            error(message('aero:atmoscira:undefinedMtype', mtype))
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end

