function [ XYZ, H, DEC, DIP, F, DXDYDZ, DH, DDEC, DDIP, DF] = igrfmagm(height,lat,lon,dyear,varargin)
% IGRFMAGM Uses International Geomagnetic Reference Field
%  [XYZ, H, DEC, DIP, F, DXDYDZ, DH, DDEC, DDIP, DF] = 
%                                     IGRFMAGM( HEIGHT, LAT, LON, DYEAR, GEN )
%  Calculates the Earth magnetic field and secular variation at a specific 
%  location and time using different generations of International Geomagnetic 
%  Reference Field. 
%
%  Inputs required by IGRF are:
%   HEIGHT :Scalar value in meters. 
%   LAT    :Scalar geodetic latitude in degrees where north latitude is 
%          positive and south latitude is negative.
%   LON    :Scalar geodetic longitude in degrees where east longitude 
%          is positive and  west longitude is negative.
%   DYEAR  :Scalar decimal year.  Decimal year is the desired year in 
%          a decimal format that includes any fraction of the year that has 
%          already passed.
%   GEN    :Scalar with the generation of the International Geomagnetic
%          Reference Field. Available options are 11 and 12. Default option
%          is 12.
%
%  Output calculated for the Earth magnetic field and secular variation include:
%   XYZ    :Magnetic field vector in nanotesla (nT). Z is the vertical component (+ve down)
%   H      :Horizontal intensity in nanotesla (nT).
%   DEC    :Declination in degrees. (+ve east)
%   DIP    :Inclination in degrees. (+ve down)
%   F      :Total intensity in nanotesla (nT).
%   DXDYDZ :Secular variation in magnetic field vector in nT/year. Z is the 
%          vertical component (+ve down) 
%   DH     :Secular variation in horizontal intensity in nT/year.
%   DDEC   :Secular variation in declination in minutes/year. (+ve east)
%   DDIP   :Secular variation in inclination in minutes/year. (+ve down)
%   DF     :Secular variation in total intensity in nT/year.
%
%   Limitations:
%
%   This function is valid between the heights of -1000 meters to 600000
%   meters. 
%
%   The 11th generation is valid between the years of 1900 and 2015, and 
%   12th generation is valid between the years of 1900 and 2020.
%
%   This function has the limitations of the International Geomagnetic
%   Reference Field (IGRF). For more information, see the IGRF web site,
%   http://www.ngdc.noaa.gov/IAGA/vmod/igrfhw.html.   
%
%   Example:
%
%   Calculate the magnetic model 1000 meters over Natick, Massachusetts on 
%   July 4, 2005 using IGRF-12:
%      [XYZ, H, DEC, DIP, F] = igrfmagm(1000, 42.283, -71.35, decyear(2005,7,4),12)
%
%   See also DECYEAR, WRLDMAGM.

%   Copyright 2015-2018 The MathWorks, Inc.

%   Reference:
%   [1] The IGRF-12 can be found on the web at 
%       http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
%   [2] Blakely, R. J., "Potential Theory in Gravity & Magnetic Applications", 
%       Cambridge University Press, Cambridge UK, 1996. 

narginchk(4, 5);

% Check that input values are expected data types
locCheckInputTypes;

if isempty(varargin)
    gen = 12;
else
    validateattributes(varargin{1},{'numeric'},{'>=',11,'<=',12,'nonnan'});
    gen = varargin{1};
end

switch gen
    case 11
        validModels = 1900:5:2010;
        warning(message('aero:igrfmagm:invalidGen',...
            '<a href="matlab:doc(''igrfmagm'')">documentation</a>'));
    case 12
        validModels = 1900:5:2015;
end
% Check that input values are within acceptable ranges
locCheckInputValues;

% Find correct model to use.  Assumes epochs are all 5 years long
tempModelToUse = validModels(dyear < validModels)-5;
if isempty(tempModelToUse)
    modelToUse = validModels(end);
    modelNext = modelToUse;
else
    modelToUse =  tempModelToUse(dyear >= tempModelToUse);
    modelNext = modelToUse + 5;
end

% Load Coefficients 
igrfdataCurrent = locLoadCoefficients(modelToUse,gen);
if (modelToUse ~= modelNext)
    igrfdataNext = locLoadCoefficients(modelNext,gen);
else
    igrfdataNext = igrfdataCurrent;
end

% Test geodetic heights
if ( isnan(height) || height < igrfdataCurrent.minAlt  || height > igrfdataCurrent.maxAlt )
    error(message('aero:igrfmagm:invalidHeight', igrfdataCurrent.minAlt, igrfdataCurrent.maxAlt));
end

% Interpolate/Extrapolate coefficients
[gha, ghb, nMax] = locInterpolateOrExtrapolate;

[X,Y,Z,XTEMP,YTEMP,ZTEMP] = locCalcField(lat, lon, height, nMax, gha, ghb);

% Calculate rest of Magnetic Field and Secular Variation
[DX,DY,DZ] = locCalcMagFieldAndSV;

% Make corrections for geographic and magnetic poles
[DX,DY,DZ] = locCorrectAtPoles(DX,DY,DZ);

% Vectorize Northward, Eastward and Downward components
XYZ = [X Y Z];
DXDYDZ = [DX DY DZ];

%==========================================================================
    function locCheckInputTypes()
        if (ischar(height) || ischar(lat) || ischar(lon) || ischar(dyear) ||...
            isstring(height) || isstring(lat) || isstring(lon) || isstring(dyear))
            error(message('aero:igrfmagm:noChar'));
        end
        
        if (~isscalar(height) || ~isscalar(lat) || ~isscalar(lon) || ~isscalar(dyear))
            error(message('aero:igrfmagm:onlyScalar'));
        end
        
        if ~isreal([height,lat,lon,dyear])
            error(message('aero:igrfmagm:noComplex'));
        end
    end
%==========================================================================
    function locCheckInputValues()
        
        dateMin = validModels(1);
        dateMax = validModels(end) + 5;
        
        
        if ( isnan(dyear) ||(dyear < dateMin) || (dyear > dateMax))
            error(message('aero:igrfmagm:invalidDYear', dateMin, dateMax));
        end
        
        if ( isnan(lat) || lat < -90.0  || lat > 90.0  )
            error(message('aero:igrfmagm:invalidLatitude'));
        end
        
        if ( isnan(lon) || lon < -180.0  || lon > 180.0  )
            error(message('aero:igrfmagm:invalidLongitude'));
        end
end
%==========================================================================
    function [gha, ghb, nMax] = locInterpolateOrExtrapolate()
        % Interpolate/Extrapolate coefficients if necessary
        if ((igrfdataNext.minYear - igrfdataCurrent.minYear) == 0)
            % Extrapolation
            factor = (dyear - igrfdataCurrent.minYear);
            factorNext = (dyear+1 - igrfdataCurrent.minYear);
            % Two orders are equal
            kk =  igrfdataCurrent.MAXORD * (igrfdataCurrent.MAXORD + 2);
            nMax = igrfdataCurrent.MAXORD;
            gha(1:kk) = igrfdataCurrent.gh(1:kk) + factor * igrfdataCurrent.sv(1:kk);
            ghb(1:kk) = igrfdataNext.gh(1:kk) + factorNext * igrfdataNext.sv(1:kk);
        else
            % Interpolation
            factor = (dyear - igrfdataCurrent.minYear) / (igrfdataNext.minYear - igrfdataCurrent.minYear);
            factorNext = (dyear+1 - igrfdataCurrent.minYear) / (igrfdataNext.minYear - igrfdataCurrent.minYear);
            if (igrfdataCurrent.MAXORD == igrfdataNext.MAXORD)
                % Two orders are equal
                kk =  igrfdataCurrent.MAXORD * (igrfdataCurrent.MAXORD + 2);
                nMax = igrfdataCurrent.MAXORD;
            else
                % current order is less than the next
                kk = igrfdataCurrent.MAXORD * (igrfdataCurrent.MAXORD + 2);
                ll = igrfdataNext.MAXORD * (igrfdataNext.MAXORD + 2);
                
                gha((kk + 1):ll) = factor * igrfdataNext.gh((kk + 1):ll);
                ghb((kk + 1):ll) = factorNext * igrfdataNext.gh((kk + 1):ll);
                
                nMax = igrfdataNext.MAXORD;
            end
            gha(1:kk) = igrfdataCurrent.gh(1:kk) + factor * (igrfdataNext.gh(1:kk) - igrfdataCurrent.gh(1:kk));
            ghb(1:kk) = igrfdataCurrent.gh(1:kk) + factorNext * (igrfdataNext.gh(1:kk) - igrfdataCurrent.gh(1:kk));
        end
    end
%==========================================================================
    function [DX,DY,DZ] = locCalcMagFieldAndSV()
        % Calculate Horizontal Intensity
        H = sqrt(X*X + Y*Y);
        
        % Calculate Total Intensity
        F = sqrt(X*X + Y*Y + Z*Z);
        
        if (F < 0.0001)
            % Declination and Inclination cannot be determined
            DEC = NaN;
            DIP = NaN;
        else
            % Calculate Inclination
            DIP = atan2(Z,H);
            if (H < 0.0001)
                % Declination cannot be determined
                DEC = NaN;
            else
                if ((H+X) < 0.0001)
                    % Set Declination
                    DEC = pi;
                else
                    % Calculate Declination
                    DEC = 2*atan2(Y,H+X);
                end
            end
        end
        
        h2 = XTEMP*XTEMP + YTEMP*YTEMP;
        argument = h2;
        hTemp = sqrt(argument);
        argument = h2 + ZTEMP*ZTEMP;
        fTemp = sqrt(argument);
        
        if (fTemp < 0.0001)
            % If DEC and DIP cannot be determined, set to NaN
            dTemp = NaN;
            iTemp = NaN;
        else
            argument = ZTEMP;
            argument2 = hTemp;
            iTemp = atan2(argument,argument2);
            if (hTemp < 0.0001)
                dTemp = NaN;
            else
                hpx = hTemp + XTEMP;
                if (hpx < 0.0001)
                    dTemp = pi;
                else
                    argument = YTEMP;
                    argument2 = hpx;
                    dTemp = 2.0 * atan2(argument,argument2);
                end
            end
        end
        
        DDEC = convang((dTemp - DEC),'rad','deg');
        if (DDEC > 180.0)
            DDEC = DDEC - 360.0;
        end
        if (DDEC <= -180.0)
            DDEC = DDEC + 360.0;
        end
        DDEC = DDEC*60.0;
        
        DDIP = convang((iTemp - DIP),'rad','deg')*60;
        DEC = convang(DEC,'rad','deg');
        DIP = convang(DIP,'rad','deg');
        
        DH = hTemp - H;
        DX = XTEMP - X;
        DY = YTEMP - Y;
        DZ = ZTEMP - Z;
        DF = fTemp - F;
    end
%==========================================================================
    function [DX,DY,DZ] = locCorrectAtPoles(DX,DY,DZ)
        if (H < 100.0)
            % at magnetic poles
            DEC = NaN;
            DDEC = NaN;
        end
        
        if (90.0-abs(lat) <= 0.001)
            % at geographic poles
            X = NaN;
            Y = NaN;
            DEC = NaN;
            DX  = NaN;
            DY  = NaN;
            DDEC = NaN;
        end
    end

%==========================================================================
    function data = locLoadCoefficients(epoch,gen)
        if gen==11 && epoch==2010
            % This will load IGRF-11 data for the 2010 period
            data = load(['aeroigrf' num2str(gen) num2str(epoch) 'data.mat']);
        else
            % This will load the data for the current version of the IGRF
            data = load(['aeroigrf' num2str(epoch) 'data.mat']);
        end
    end
%==========================================================================
end

function [x,y,z,xtemp,ytemp,ztemp] = locCalcField(latitude, longitude, elev, nMax, gha, ghb)
% Calculates field components from spherical harmonic 

% WGS84 
Re = 6371.2;
a2 = 40680631.59;      
b2 = 40408299.98;           

% preallocate vectors for calculations
sl = zeros(14,1);
cl = zeros(14,1);
p = zeros(119,1);
q = zeros(119,1);

% initialize outputs
x = 0;
y = 0;
z = 0;
xtemp = 0;
ytemp = 0;
ztemp = 0;

% initialize counters
l = 1;
n = 0;
m = 1;

% convert to kilometers from meters
elevKM = elev*0.001;

slat_gd = sin( convang(latitude,'deg','rad') );

if ((90.0 - latitude) < 0.001)    
    %  300 ft. from North pole
    latitudeCrt = 89.999;
else   
    if ((90.0 + latitude) < 0.001)        
        %  300 ft. from South pole
        latitudeCrt = -89.999;
    else
        latitudeCrt = latitude;
    end
end

clat_gd = cos( convang(latitudeCrt,'deg','rad') );

rlon = convang(longitude,'deg','rad');
sl(1) = sin( rlon );
cl(1) = cos( rlon );

npq = (nMax * (nMax + 3)) / 2;

% convert from geodetic to geocentric coordinates
aa = a2 * clat_gd * clat_gd;
bb = b2 * slat_gd * slat_gd;
cc = aa + bb;
dd = sqrt( cc );
r = sqrt( elevKM * (elevKM + 2.0 * dd) + (a2 * aa + b2 * bb) / cc );
cd = (elevKM + dd) / r;
sd = (a2 - b2) / dd * slat_gd * clat_gd / r;
slat = slat_gd * cd - clat_gd * sd;
clat = clat_gd * cd + slat_gd * sd;

ratio = Re / r;
sqrt3 = sqrt( 3.0 );

p(1) = 2.0 * slat;
p(2) = 2.0 * clat;
p(3) = 4.5 * slat * slat - 1.5;
p(4) = 3.0 * sqrt3 * clat * slat;
q(1) = -clat;
q(2) = slat;
q(3) = -3.0 * clat * slat;
q(4) = sqrt3 * (slat * slat - clat * clat);

for k = 1:npq
    if (n < m)        
        m = 0;
        n = n + 1;
        power =  n + 2;
        rr = ratio^power;
        fn = n;
    end
    fm = m;
    if (k >= 5)
        if (m == n)
            aa = sqrt( (1.0 - 0.5/fm) );
            j = k - n - 1;
            p(k) = (1.0 + 1.0/fm) * aa * clat * p(j);
            q(k) = aa * (clat * q(j) + slat/fm * p(j));
            sl(m) = sl(m-1) * cl(1) + cl(m-1) * sl(1);
            cl(m) = cl(m-1) * cl(1) - sl(m-1) * sl(1);
        else            
            aa = sqrt( fn*fn - fm*fm );
            bb = sqrt( ((fn - 1.0)*(fn-1.0)) - (fm * fm) )/aa;
            cc = (2.0 * fn - 1.0)/aa;
            ii = k - n;
            j = k - 2 * n + 1;
            p(k) = (fn + 1.0) * (cc * slat/fn * p(ii) - bb/(fn - 1.0) * p(j));
            q(k) = cc * (slat * q(ii) - clat/fn * p(ii)) - bb * q(j);
        end
    end
    aa = rr * gha(l);
    aaTemp = rr * ghb(l);
    if (m == 0)
        x = x + aa * q(k);
        z = z - aa * p(k);
        
        xtemp = xtemp + aaTemp * q(k);
        ztemp = ztemp - aaTemp * p(k);
        l = l + 1;     
    else
        
        bb = rr * gha(l+1);
        cc = aa * cl(m) + bb * sl(m);
        x = x + cc * q(k);
        z = z - cc * p(k);
        
        bbTemp = rr * ghb(l+1);
        ccTemp = aaTemp * cl(m) + bbTemp * sl(m);
        xtemp = xtemp + ccTemp * q(k);
        ztemp = ztemp - ccTemp * p(k);

        if (clat > 0)
            y = y + (aa * sl(m) - bb * cl(m)) * fm * p(k)/((fn + 1.0) * clat);
            ytemp = ytemp + (aaTemp * sl(m) - bbTemp * cl(m)) * fm * p(k)/((fn + 1.0) * clat);
        else
            y = y + (aa * sl(m) - bb * cl(m)) * q(k) * slat;
            ytemp = ytemp + (aaTemp * sl(m) - bbTemp * cl(m)) * q(k) * slat;
        end
        
        l = l + 2;
    end
    m = m + 1;
end
x_old = x;
x = x * cd + z * sd;
z = z * cd - x_old * sd;

xtemp_old = xtemp;
xtemp = xtemp * cd + ztemp * sd;
ztemp = ztemp * cd - xtemp_old * sd;
end
