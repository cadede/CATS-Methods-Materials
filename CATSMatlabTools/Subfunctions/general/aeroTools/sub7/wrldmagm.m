function [XYZ, H, D, I, F] = wrldmagm(height, lat, lon, dyear, varargin )
%WRLDMAGM calculates the Earth's magnetic field at a specific location and
%time using the World Magnetic Model (WMM).
%
% [XYZ, H, D, I, F] = WRLDMAGM(HEIGHT, LAT, LON, DYEAR) calculates
% the Earth's magnetic field at a specific location and time using the
% World Magnetic Model (WMM) WMM2020.
%
% [XYZ, H, D, I, F] = WRLDMAGM(HEIGHT, LAT, LON, DYEAR, MODEL)
% calculates the Earth's magnetic field using World Magnetic Model MODEL.
%
% [XYZ, H, D, I, F] = WRLDMAGM(HEIGHT, LAT, LON, DYEAR, 'Custom', FILE)
% calculates the Earth's magnetic field using the World Magnetic Model
% defined in the WMM.cof file, FILE. WMM.COF files must be in their original
% form as provided by NOAA (http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml).
%
%   Inputs to wrldmagm are:
%   HEIGHT :Scalar value in meters. 
%   LAT    :Scalar geodetic latitude in degrees where north latitude is 
%           positive and south latitude is negative.
%   LON    :Scalar geodetic longitude in degrees where east longitude 
%           is positive and west is negative.
%   DYEAR  :Scalar decimal year.  Decimal year is the desired year in 
%           a decimal format to include any fraction of the year that has 
%           already passed.
%   MODEL  :Scalar value or string specifying the WMM version to use in
%           calculations. MODEL input can be one of the following: 
%            MODEL:                     Coefficients used:
%            '2020' or 2020             WMM2020     (2020-2025) (default)
%            '2015v2', '2015', or 2015  WMM2015v2   (2015-2020)
%            '2015v1'                   WMM2015(v1) (2015-2020) (deprecated)
%            '2010' or 2010             WMM2010     (2010-2015)
%            '2005' or 2005             WMM2005     (2005-2010)
%            '2000' or 2000             WMM2000     (2000-2005)
%            'Custom'                   Custom WMM.COF file
%   FILE  :String specifying the WMM coefficient file to use. WMM.COF
%          files must be in their original form as provided by NOAA
%          (http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml). FILE is
%          available when MODEL is 'Custom'.
%
%   Output calculated by wrldmagm are:
%   XYZ    :Magnetic field vector in nanotesla (nT). 
%   H      :Horizontal intensity in nanotesla (nT). 
%   D      :Declination in degrees. 
%   I      :Inclination in degrees. 
%   F      :Total intensity in nanotesla (nT).
%
%   Limitations:
%
%   This function has the limitations of the World Magnetic Model(WMM). The
%   WMM2020 is valid between -1 and 850km, as outlined in the World
%   Magnetic Model 2020 Technical Report. For more information see the
%   references section in the documentation.
%
%   WMM2015v2 superseded WMM2015(v1). Consider replacing WMM2015(v1) with
%   WMM2015v2 when used for navigation and other systems between January 1,
%   2015 and December 31, 2019. WMM2015v2 was released by NOAA in February,
%   2019 to correct performance degradation issues in the Artic region.
%   Therefore, it is acceptable to use WMM2015(v1) in systems below
%   55-degrees latitude in the Northern hemisphere for January 1, 2015 to
%   December 31, 2019.
%
%   Example:
%
%   Calculate the magnetic model 1000 meters over Natick, Massachusetts on 
%   July 4, 2020 using WMM2020:
%      [XYZ, H, D, I, F] = wrldmagm(1000, 42.283, -71.35,...
%          decyear(2020,7,4))
%
%      [XYZ, H, D, I, F] = wrldmagm(1000, 42.283, -71.35,...
%          decyear(2020,7,4), '2020')
%
%   Calculate the magnetic model 1000 meters over Natick, Massachusetts on
%   July 4, 2020, using downloaded WMM.COF file:
%      [XYZ, H, D, I, F] = wrldmagm(1000, 42.283, -71.35,...
%          decyear(2020,7,4), 'Custom', 'WMM.COF')
%
%   See also DECYEAR, IGRFMAGM.

%   Copyright 2000-2020 The MathWorks, Inc.

%   Limitations:  The WMM specification produces data that is reliable 
%   five years after the epoch of the model, which begins January 1, of the 
%   model year selected.  The WMM specification describes only the 
%   long-wavelength spatial magnetic fluctuations due to the Earth's core. 
%   Intermediate and short-wavelength fluctuations, contributed from the 
%   crustal field (the mantle and crust), are not included. Also, the 
%   substantial fluctuations of the geomagnetic field, which occur constantly
%   during magnetic storms and almost constantly in the disturbance field 
%   (auroral zones), are not included.

%   Reference:
%   [1] The WMM2020 can be found on the web at
%   http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml

narginchk(4, 6);
try [varargin{:}] = convertStringsToChars(varargin{:}); catch; end;

% Use 2020 epoch or catch numeric inputs for epoch
if (nargin < 5)
    varargin{1} = 2020;
end
switch lower(varargin{1})
    case {2020,'2020'}
        epoch = 2020;
        load aerowmm2020; %#ok<*LOAD>
    case {2015,'2015','2015v2'}
        warning(message('aero:wrldmagm:obsoleteEpoch','2020'));
        epoch = 2015;
        load aerowmm2015v2;
    case '2015v1'
        warning(message('aero:wrldmagm:obsoleteEpoch2015v1'));
        epoch = 2015;
        load aerowmm2015v1;
    case {2010,'2010'}
        warning(message('aero:wrldmagm:obsoleteEpoch','2010'));
        epoch = 2010;
        load aerowmm2010;
    case {2005,'2005'}
        warning(message('aero:wrldmagm:obsoleteEpoch','2005'));
        epoch = 2005;
        load aerowmm2005;
    case {2000,'2000'}
        warning(message('aero:wrldmagm:obsoleteEpoch','2000'));
        epoch = 2000;
        load aerowmm2000;
    case 'custom'
        if nargin > 5
            [~,~,ext] = fileparts(varargin{2});
            if ~strcmpi(ext,'.cof')
                error(message('aero:wrldmagm:invalidExtension'));
            end
            [c,dc,epoch,fm,fn,k,maxdef,~] = readWMMCoeff(varargin{2});
        else
            error(message('aero:wrldmagm:missingCustomFile'));
        end
    otherwise
        error(message('aero:wrldmagm:invalidEpoch'));
end

% Calculate the time difference from the epoch
dt = dyear - epoch;

if ((nargin < 5) && ( dt < 0.0 && dt >= -5.0 ))
    error(message('aero:wrldmagm:updateEpoch','WMM2020'));
end

if ( isnan(dt) || dt < 0.0 || dt > 5.0)
    error(message('aero:wrldmagm:invalidDYear'));
end

if ( isnan(height) || ( epoch <= 2010 && ( height < 0.0  || height > 1000000.0 ) ) )
    warning(message('aero:wrldmagm:invalidHeight', epoch));
elseif ( isnan(height) || ( epoch > 2010 && ( height < -1000.0  || height > 850000.0 ) ) )
    warning(message('aero:wrldmagm:invalidHeightAfter2010', epoch));
end

% validateattributes(lat, 'numeric', {'nonnan', 'finite'});
% validateattributes(lon, 'numeric', {'nonnan', 'finite'});

% wrap latitude and longitude if needed
[~, lat, lon] = wraplatitude( lat, lon, 'deg' );

% check and fix angle wrapping in longitude
[~, lon] = wraplongitude( lon, 'deg', '180' );

maxord = maxdef + 1;

tc = zeros(maxord,maxord);
dp = zeros(maxord,maxord);
pp = ones(maxord,1);
snorm = zeros(maxord*maxord,1);
sp = zeros(1,maxord); 
cp = zeros(1,maxord);
dp(1,1) = 0.0;
cp(1) = 1.0; 
snorm(1) = 1.0; 
pp(1) = 1.0;

% convert to kilometers
height = height*0.001;

a = 6378.137;
b = 6356.7523142;
re = 6371.2;
a2 = a*a;
b2 = b*b;
c2 = a2-b2;
a4 = a2*a2;
b4 = b2*b2;
c4 = a4 - b4;

rlon = deg2rad(lon);
rlat = deg2rad(lat);
srlon = sin(rlon);
srlat = sin(rlat);
crlon = cos(rlon);
crlat = cos(rlat);
srlat2 = srlat*srlat;
crlat2 = crlat*crlat;
sp(2) = srlon;
cp(2) = crlon;

% convert from geodetic coordinates to spherical coordinates
q = sqrt(a2-c2*srlat2);
q1 = height*q;
q2 = ((q1+a2)/(q1+b2))*((q1+a2)/(q1+b2));
ct = srlat/sqrt(q2*crlat2+srlat2);
st = sqrt(1.0-(ct*ct));
r2 = (height*height)+2.0*q1+(a4-c4*srlat2)/(q*q);
r = sqrt(r2);
d = sqrt(a2*crlat2+b2*srlat2);
ca = (height+d)/r;
sa = c2*crlat*srlat/(r*d);

for m = 3:maxord
    sp(m) = sp(2)*cp(m-1)+cp(2)*sp(m-1);
    cp(m) = cp(2)*cp(m-1)-sp(2)*sp(m-1);
end
aor = re/r;
ar = aor*aor;
br = 0.0; bt = 0.0; bp = 0.0; bpp = 0.0;
for n = 1:maxord-1
    ar = ar*aor;
    for m = 0:n
        %     Compute unnormalized associated legendre polynomials
        %     and derivatives via recursion relations
        if (n == m)
            snorm(n+1 + m*maxord) = st*snorm(n + (m-1)*maxord);
            dp(m+1, n+1) = st*dp(m, n)+ct*snorm(n +(m-1)*maxord);
        end
        if (n == 1 && m == 0)
            snorm(n+1) = ct*snorm(n);
            dp(m+1, n+1) = ct*dp(m+1, n)-st*snorm(n);
        end
        if (n > 1 && n ~= m)
            if (m > n-2)
                snorm(n-2+m*maxord) = 0.0;
                dp(m+1, n-1) = 0.0;
            end
            snorm(n+1 + m*maxord) = ct*snorm(n + m*maxord)-k(m+1, n+1)*snorm(n-1 + m*maxord);
            dp(m+1, n+1) = ct*dp(m+1, n) - st*snorm(n + m*maxord)-k(m+1, n+1)*dp(m+1, n-1);
        end

        %   Time adjust the Gauss coefficients
        tc(m+1, n+1) = c(m+1, n+1)+dt*dc(m+1, n+1);
        if (m ~= 0)
            tc(n+1, m) = c(n+1, m)+dt*dc(n+1, m);
        end

        %    Accumulate terms of the spherical harmonic expansions
        par = ar*snorm(n+1+m*maxord);
        if (m == 0)
            temp1 = tc(m+1, n+1)*cp(m+1);
            temp2 = tc(m+1, n+1)*sp(m+1);
        else
            temp1 = tc(m+1, n+1)*cp(m+1)+tc(n+1, m)*sp(m+1);
            temp2 = tc(m+1, n+1)*sp(m+1)-tc(n+1, m)*cp(m+1);
        end
        bt = bt-ar*temp1*dp(m+1, n+1);
        bp = bp+(fm(m+1)*temp2*par);
        br = br+(fn(n+1)*temp1*par);

        %   Special Case:  North/South geographic poles

        if (st == 0.0 && m == 1)
            if (n == 1)
                pp(n+1) = pp(n);
            else
                pp(n+1) = ct*pp(n)-k(m+1,n+1)*pp(n-1);
            end
            parp = ar*pp(n+1);
            bpp = bpp + (fm(m+1)*temp2*parp);
        end

    end
end
if (st == 0.0)
    bp = bpp;
else
    bp = bp/st;
end

%    Rotate magnetic vector components from spherical to geodetic
%    coordinates

bx = -bt*ca-br*sa;
by = bp;
bz = bt*sa-br*ca;

%   Compute declination (D), inclination (I) & total intensity (F)

bh = sqrt((bx*bx)+(by*by));
F = sqrt((bh*bh)+(bz*bz));
D = atan2(by,bx);
I = atan2(bz,bh);
%   Compute XYZ & H components of the magnetic field

X = F*cos(D)*cos(I);
Y = F*cos(I)*sin(D);
Z = F*sin(I);
XYZ = [X; Y; Z];

H = F*cos(I);
D = rad2deg(D);
I = rad2deg(I);
end

