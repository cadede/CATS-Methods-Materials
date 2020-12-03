function varargout = dcmeci2ecef(reduction,UTC,varargin)
%   DCMECI2ECEF Convert Earth-centered inertial (ECI) to Earth-centered
%   Earth-fixed (ECEF) coordinates.
%   DCM = DCMECI2ECEF( REDUCTION,UTC,DELTAAT,DELTAUT1,POLARMOTION,'ADDPARAMNAME',ADDPARAMVALUE )
%   calculates the position direction cosine matrix (ECI to ECEF), for
%   given set of time and geophysical data, UTC, DELTAAT, DELTAUT1, 
%   POLARMOTION, and DNUTATION or DCIP.
%
%   Inputs arguments for DCMECI2ECEF are:
%   REDUCTION:      String indicating the reduction process through which
%                   the function calculates the direction cosine matrix.
%                   It can either be IAU-76/FK5 (which uses the IAU 1976
%                   Precession Model and the IAU 1980 Theory of Nutation,
%                   which is no longer current but some programs still use
%                   this reduction) or IAU-2000/2006 (which uses the P03
%                   precession model). The reduction method you select
%                   determines the ADDPARAMNAME parameter pair
%                   characteristics. The IAU-76/FK5 method returns a
%                   transformation matrix that does not have the
%                   characteristics of a DCM due to the polar motion
%                   approximation. 
%   UTC:            Array of Universal Coordinated Time (UTC) in year,
%                   month, day, hour, minutes and seconds for which the
%                   function calculates the transformation matrix. Define
%                   the array as one of the following: Array with 1 row and
%                   6 columns or M by 6 array for M transformation
%                   matrices, one for each UTC date. Values for year,
%                   month, day, hour and minutes should be whole numbers.
%   DELTAAT:        Difference in seconds between the International Atomic
%                   Time (TAI) and UTC. It can be defined as either a
%                   scalar or a one dimensional array with M elements (if M
%                   UTC dates defined) for equal number of transformation
%                   matrices. The default value is an M by 1 null array. 
%   DELTAUT1:       Difference, in seconds, between UTC and Universal Time
%                   (UT1). It can be defined as either a scalar or a one
%                   dimensional array with M elements (if M UTC dates 
%                   defined) for equal number of transformation matrices.
%                   The default value is an M by 1 null array. 
%   POLARMOTION:    Polar displacement due to the motion of the Earth's
%                   crust, in radians, along the x and y axis. It can be
%                   defined as either a 1 by 2 array or an M by 2 (if M UTC
%                   dates defined) for equal number of transformation
%                   matrices. The default value is an M by 2 null array.
%
%   Input parameter Name/Value pair for 'ADDPARAMNAME' and ADDPARAMVALUE
%   are: 
%   'DNUTATION':    (IAU-76/FK5 reduction only) M by 2 array for the
%                   adjustment in radians to the longitude (dDeltaPsi) and
%                   obliquity (dDeltaEpsilon). The default value is an M by
%                   2 null array.
%   'DCIP':         (IAU-2000/2006 reduction only) M by 2 array for the
%                   adjustment in radians to the location of the Celestial
%                   Intermediate Pole (CIP) along the x (dDeltaX) and y
%                   (dDeltaY) axis. The default value is an M by 2 null
%                   array.
%
%   For historical values for DNUTATION and DCIP, see the International
%   Earth Rotation and Reference Systems Service (IERS) website
%   (http://www.iers.org) under the 'Earth Orientation Data' product. 
%
%   Output calculated by DCMECI2ECEF is:
%   DCM:    3 by 3 by M array where M is the number of dates for which the
%           function calculates the transformation matrix. 
%
%   Examples:
%   Calculate the position direction cosine matrix (ECI to ECEF), using the
%   IAU-76/FK5 reduction, for one UTC date with all parameters defined.
%
%   DCM = dcmeci2ecef('IAU-76/FK5',[2000 1 12 4 52 12.4],32,0.234,[-0.0682e-5 ...
%                    0.1616e-5],'dNutation',[-0.2530e-6 -0.0188e-6])
%
%   Calculate the position direction cosine matrix (ECI to ECEF), using the
%   IAU-2000/2006 reduction, for two UTC dates. All other parameters
%   default to null arrays.
%
%   DCM = dcmeci2ecef('IAU-2000/2006',[2000 1 12 4 52 12.4;2000 1 12 4 52 13])
%
%
%   See also LLA2ECEF, ECEF2LLA, GEOC2GEOD, GEOD2GEOC, LLA2ECI, ECI2LLA,
%   DELTAUT1, DELTACIP, POLARMOTION.

%   Copyright 2012-2018 The MathWorks, Inc.

%   References:
%   [1]Petit, Gerard and Brian Luzum. IERS Conventions (2010) - IERS
%   Technical note 36. Verlag des Bundesamts fur Kartographie und Geodasie.
%   Frankfurt am Main, 2010.
%   [2]Vallado, D. A., Fundamentals of Astrodynamics and Applications,
%      McGraw-Hill, New York, 1997.

%% Validate i/o
% Validate outputs
nargoutchk(0,1)
% Validate amount
narginchk(2,7)
% Validate date
validateattributes(UTC,{'numeric'},{'ncols',6,'real','finite','nonnan'})
% Assign date vectors
year = UTC(:,1);
month = UTC(:,2);
day = UTC(:,3);
hour = UTC(:,4);
min = UTC(:,5);
sec = UTC(:,6);
% Validate vectors
if any(year<1)
    error(message('aero:dcmeci2ecef:invalidYear'));
end
if any(month<1) || any(month>12)
    error(message('aero:dcmeci2ecef:invalidMonth'));
end
if any(day<1) || any(day>31)
    error(message('aero:dcmeci2ecef:invalidDay'));
end
if any(hour<0) || any(hour>24)
    error(message('aero:dcmeci2ecef:invalidHour'));
end
if any(min<0) || any(min>60)
    error(message('aero:dcmeci2ecef:invalidMin'));
end
if any(sec<0) || any(sec>60)
    error(message('aero:dcmeci2ecef:invalidSec'));
end
% Check for year fractions
if any(mod(year,1)~=0)
    year = fix(year);
end
if any(mod(month,1)~=0)
    month = fix(month);
end
if any(mod(day,1)~=0)
    day = fix(day);
end
if any(mod(hour,1)~=0)
    hour = fix(hour);
end
if any(mod(min,1)~=0)
    min = fix(min);
end
len = length(year);
% Parse and validate the inputs
p = inputParser;
addRequired(p,'reduction',@(x) validateattributes(x,{'char','string'},{'nonempty'}));
addRequired(p,'UTC',@(x) validateattributes(x,{'numeric'},{'ncols',6,'real',...
    'finite','nonnan'}));
addOptional(p,'deltaAT',zeros(len,1),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,1]}));
addOptional(p,'deltaUT1',zeros(len,1),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,1]}));
addOptional(p,'polarMotion',zeros(len,2),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,2]}));
addParameter(p,'dNutation',zeros(len,2),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,2]}));
addParameter(p,'dCIP',zeros(len,2),@(x) validateattributes(x,{'numeric'},...
    {'real','finite','nonnan','size',[len,2]}));
parse(p,reduction,UTC,varargin{:});
%Validate reduction
reduction = p.Results.reduction;
validReduction = {'IAU-2000/2006','IAU-76/FK5'};
reduction = lower(validatestring(reduction,validReduction));
%Validate that the additional parameter matches the reduction method
if any(strcmp(p.UsingDefaults,'dNutation')) && ~any(strcmp(p.UsingDefaults,'dCIP')) &&...
        ~strcmp(reduction,'iau-2000/2006')
    error(message('aero:dcmeci2ecef:invalidDNutation'));
elseif ~any(strcmp(p.UsingDefaults,'dNutation')) && any(strcmp(p.UsingDefaults,'dCIP')) && ...
        ~strcmp(reduction,'iau-76/fk5')
    error(message('aero:dcmeci2ecef:invalidDCIP'));
end
%Validate that both reduction methods are not defined (just one should be
%defined)
if ~any(strcmp(p.UsingDefaults,'dNutation')) && ~any(strcmp(p.UsingDefaults,'dCIP'))
    error(message('aero:dcmeci2ecef:invalidDNutationDCIP'))
end
%% Common time calculations
% Calculations determining the number of Julian centuries for terrestrial
% time (tTT) and UT1.

% Seconds for UT1
ssUT = sec + p.Results.deltaUT1; 
% Seconds for UTC
ssTT = sec + p.Results.deltaAT + 32.184; 
% Julian date for terrestrial time
jdTT = mjuliandate(year,month,day,hour,min,ssTT);
% Number of Julian centuries since J2000 for terrestrial time.
tTT = (jdTT - 51544.5)/36525;
tTT2 = tTT.*tTT;
tTT3 = tTT2.*tTT;

switch reduction
    case 'iau-2000/2006'
        %% CIO Based Transformation
        % This transformation is based on the IERS technical note 36. 
        % Additional time calculations
        tTT = tTT';
        tTT2 = tTT2';
        tTT3 = tTT3';
        tTT4 = tTT3.*tTT;
        tTT5 = tTT4.*tTT;
        % Julian date for UT1
        jdUT1 = mjuliandate(year,month,day,hour,min,ssUT);
        % Elapsed Julian days since J2000
        jdElapsed = jdUT1 - 51544.5;
        % Julian day fraction
        jdFraction = mod(jdElapsed,1);
        
        %% Polar motion
        % TIO locator
        sp = convang(-0.000047*tTT/3600,'deg','rad');
        % Transformation matrix for polar motion
        W = angle2dcm(sp,-p.Results.polarMotion(:,1),-p.Results.polarMotion(:,2),'ZYX');
        
        %% Earth rotation
        % Earth rotation angle
        thetaERA = mod(2*pi*(jdFraction + 0.7790572732640 + 0.00273781191135448*jdElapsed),2*pi);
        % Transformation matrix for earth rotation
        R = angle2dcm(thetaERA,zeros(len,1),zeros(len,1));
        
        %% Celestial Motion of the CIP
        % Arguments for lunisolar nutation
        mMoon = 485868.249036 + 1717915923.2178*tTT + 31.8792*tTT2 + 0.051635*tTT3 - 0.00024470*tTT4;
        mSun = 1287104.793048 + 129596581.0481*tTT - 0.5532*tTT2 + 0.000136*tTT3 - 0.00001149*tTT4;
        umMoon = 335779.526232 + 1739527262.8478*tTT - 12.7512*tTT2 - 0.001037*tTT3 + 0.00000417*tTT4;
        dSun = 1072260.703692 + 1602961601.2090*tTT - 6.3706*tTT2 + 0.006593*tTT3 - 0.00003169*tTT4;
        omegaMoon = 450160.398036 - 6962890.5431*tTT + 7.4722*tTT2 + 0.007702*tTT3 - 0.00005939*tTT4;
        % Arguments for planetary nutation
        lMercury = 4.402608842 + 2608.7903141574*tTT;
        lVenus = 3.176146697 + 1021.3285546211*tTT;
        lEarth = 1.753470314 + 628.3075849991*tTT;
        lMars = 6.203480913 + 334.06124267*tTT;
        lJupiter = 0.599546497 + 52.9690962641*tTT;
        lSaturn = 0.874016757 + 21.329910496*tTT;
        lUranus = 5.481293872 + 7.4781598567*tTT;
        lNeptune = 5.311886287 + 3.8133035638*tTT;
        pa = 0.02438175*tTT + 0.00000538691*tTT2;
        % Vector arrangement for series evaluation
        nutationV = mod([convang([mMoon; mSun; umMoon; dSun; omegaMoon]/3600,'deg','rad'); ...
            lMercury; lVenus; lEarth; lMars; lJupiter; lSaturn; lUranus; lNeptune; pa],2*pi);
                
        % Coordinates of the CIP in the ICRS
        % Series data (according to Table 5.2a, 5.2b, and Table 5.2d)
        load aeroCIP2006.mat
        % Polynomial part of X and Y
        X0 = -16617 + 2004191898*tTT - 429782.9*tTT2 - 198618.34*tTT3 + 7.578*tTT4 + 5.9285*tTT5;
        Y0 = -6951 - 25896*tTT - 22407274.7*tTT2 + 1900.59*tTT3 + 1112.526*tTT4 + 0.1358*tTT5;
        % Polynomial part of S
        S0 = 94 + 3808.65*tTT - 122.68*tTT2 - 72574.11*tTT3 + 27.98*tTT4 + 15.62*tTT5;
        % Series evaluation
        % For X:
        FX = zeros(length(aeroCIP2006.X),len);
        FX(1:1306,:) = ones(1306,len);
        FX(1307:1559,:) = repmat(tTT,[253 1]);
        FX(1560:1595,:) = repmat(tTT2,[36 1]);
        FX(1596:1599,:) = repmat(tTT3,[4 1]);
        FX(1600,:) = tTT4;
        argX = aeroCIP2006.X(:,4:17)*nutationV;
        X = sum(cell2mat(arrayfun(@(k) (aeroCIP2006.X(:,2).*sin(argX(:,k)) ...
            + aeroCIP2006.X(:,3).*cos(argX(:,k))).*FX(:,k),1:len,'UniformOutput',false)));
        % For Y:
        FY = zeros(length(aeroCIP2006.Y),len);
        FY(1:962,:) = ones(962,len);
        FY(963:1239,:) = repmat(tTT,[277 1]);
        FY(1240:1269,:) = repmat(tTT2,[30 1]);
        FY(1270:1274,:) = repmat(tTT3,[5 1]);
        FY(1275,:) = tTT4;
        argY = aeroCIP2006.Y(:,4:17)*nutationV;
        Y = sum(cell2mat(arrayfun(@(k) (aeroCIP2006.Y(:,2).*sin(argY(:,k)) ...
            + aeroCIP2006.Y(:,3).*cos(argY(:,k))).*FY(:,k),1:len,'UniformOutput',false)));
        % For S:
        FS = zeros(length(aeroCIP2006.S),len);
        FS(1:33,:) = ones(33,len);
        FS(34:36,:) = repmat(tTT,[3 1]);
        FS(37:61,:) = repmat(tTT2,[25 1]);
        FS(62:65,:) = repmat(tTT3,[4 1]);
        FS(66,:) = tTT4;
        argS = aeroCIP2006.S(:,4:11)*[nutationV(1:5,:);nutationV(7:8,:);nutationV(14,:)];
        S = sum(cell2mat(arrayfun(@(k) (aeroCIP2006.S(:,2).*sin(argS(:,k)) ...
            + aeroCIP2006.S(:,3).*cos(argS(:,k))).*FS(:,k),1:len,'UniformOutput',false)));
        % Adding the elements:
        X = X + X0;
        Y = Y + Y0;
        S = S + S0;
        % Convert from microarcseconds to radians
        X = convang(X*1e-6/3600,'deg','rad') + p.Results.dCIP(:,1)';
        Y = convang(Y*1e-6/3600,'deg','rad') + p.Results.dCIP(:,2)';
        S = convang(S*1e-6/3600,'deg','rad')-X.*Y/2;
        % Coordinates of the CIP
        E = atan2(Y,X);
        d = atan(sqrt((X.^2+Y.^2)./(1-X.^2-Y.^2)));
        % Transformation matrix for celestial motion of the CIP
        Q = angle2dcm(E,d,-E-S,'ZYZ');
        
        %% Matrix calculations
        tmp = arrayfun(@(k) ((W(:,:,k)*R(:,:,k)*Q(:,:,k))) ,1:len,'UniformOutput',false);
        DCM = reshape(cell2mat(tmp),3,3,len);
            
    case 'iau-76/fk5'
        %% IAU-76/FK-5 based transformation
        % This transformation is based on the process described by Vallado
        % (originally McCarthy). 
        
        % Additional time calculations:
        jdUT1 = mjuliandate(year,month,day);
        tUT1 = (jdUT1 - 51544.5)/36525;
        tUT12 = tUT1.*tUT1;
        tUT13 = tUT12.*tUT1;
       
        %% Sidereal time
        % Greenwich mean sidereal time at midnight
        thGMST0h = 100.4606184 + 36000.77005361*tUT1 + 0.00038793*tUT12 - 2.6e-8*tUT13;
        % Ratio of universal to sidereal time
        omegaPrec = 1.002737909350795 + 5.9006e-11*tUT1 - 5.9e-15*tUT12;
        % Elapsed universal time since midnight to the time of the
        % observation
        UT1 = hour*60*60 + min*60 + ssUT;
        % Greenwich mean sidereal time at time of the observation
        thGMST = mod(thGMST0h + (360/(24*60*60))*omegaPrec.*UT1,360);
        
        %% Nutation
        % Mean obliquity of the ecliptic
        epsilonBar = convang(23.439291 - 0.0130042*tTT - 1.64e-7*tTT2 + 5.04e-7*tTT3,'deg','rad');
        % Nutation angles obtained using JPL data
        nutationAngles = earthNutation(2400000.5+jdTT);
        dpsi = nutationAngles(:,1); %Nutation in Longitude
        depsilon = nutationAngles(:,2); %Nutation in obliquity
        % Adjustments to nutation angles (provided from real measurements)
        dpsi = dpsi + p.Results.dNutation(:,1);
        depsilon = depsilon + p.Results.dNutation(:,2);
        %The last two terms for equation of the equinoxes are only included
        %if the date is later than January 1, 1997 (MJD=50449)
        omegaMoon = convang(125.04455501 - (5*360 + 134.1361851)*tTT ...
                                    + 0.0020756*tTT2 + 2.139e-6*tTT3,'deg','rad');
        omegaMoon(jdUT1<50449) = 0;
        % Equation of the equinoxes
        equinoxEq = dpsi.*cos(epsilonBar) + convang(0.00264/3600,'deg','rad')*...
            sin(omegaMoon) + convang(0.000063/3600,'deg','rad')*sin(2*omegaMoon); 
        % Greenwhich apparent sidereal time
        thGAST = thGMST*pi/180 + equinoxEq;
        % Transformation matrix for earth rotation
        R = angle2dcm(thGAST,zeros(len,1),zeros(len,1));
        % True obliquity of ecliptic
        epsilon = epsilonBar+depsilon;
        % True equator to mean equinox date transformation matrix
        N = angle2dcm(epsilonBar,-dpsi,-epsilon,'XZX');
        %% Precession
        % Zeta, theta and z represent the combined effects of general
        % precession.
        zeta = convang((2306.2181*tTT + 0.30188*tTT2 + 0.017998*tTT3)/3600,'deg','rad'); 
        theta = convang((2004.3109*tTT - 0.42665*tTT2 - 0.041833*tTT3)/3600,'deg','rad');
        z = convang((2306.2181*tTT + 1.09468*tTT2 + 0.018203*tTT3)/3600,'deg','rad'); 
        % Mean equinox to celestial reference frame
        P = angle2dcm(-zeta,theta,-z,'ZYZ');
        %% Polar motion
        W = repmat(eye(3),1,len);
        W = reshape(W,3,3,len);
        W(1,3,:) = p.Results.polarMotion(:,1);
        W(3,1,:) = -p.Results.polarMotion(:,1);
        W(2,3,:) = -p.Results.polarMotion(:,2);
        W(3,2,:) = p.Results.polarMotion(:,2);
        
        %% Matrix calculations
        tmp = arrayfun(@(k) ((W(:,:,k))*R(:,:,k)*N(:,:,k)*P(:,:,k)),1:len,'UniformOutput',false);
        DCM = reshape(cell2mat(tmp),3,3,len);
end
%% Outputs
varargout{1} = DCM;
