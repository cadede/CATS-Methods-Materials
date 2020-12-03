function [gx, gy, gz]= gravitysphericalharmonic( p, varargin )
%  GRAVITYSPHERICALHARMONIC Implement a spherical harmonic representation
%   of planetary gravity. 
%   [GX, GY, GZ] = GRAVITYSPHERICALHARMONIC( P ) implements the mathematical
%   representation of spherical harmonic planetary gravity based on
%   planetary gravitational potential. Using P, a M-by-3 array of
%   Planet-Centered Planet-Fixed coordinates, GX, GY and GZ, arrays of M 
%   gravity values in the x-axis, y-axis and z-axis of the Planet-Centered
%   Planet-Fixed coordinates are calculated for planet using 120th degree 
%   and order spherical coefficients for EGM2008 by default. 
%
%   Alternate formats for calling spherical harmonic gravity are:
%   [GX, GY, GZ] = GRAVITYSPHERICALHARMONIC( P, DEGREE )   
%   [GX, GY, GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL )   
%   [GX, GY, GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL, DEGREE )   
%   [GX, GY, GZ] = GRAVITYSPHERICALHARMONIC( P, MODEL, DEGREE, ACTION )   
%   [GX, GY, GZ] = GRAVITYSPHERICALHARMONIC( P, 'Custom', DEGREE, {DATAFILE DFREADER}, ACTION )   
%
%   Inputs for spherical harmonic gravity are:
%   P        :a M-by-3 array of Planet-Centered Planet-Fixed coordinates in
%            meters where the z-axis is positive towards the North Pole. For
%            Earth this would be ECEF coordinates.
%   MODEL    :a string specifying the planetary model:
%            'EGM2008' (Earth), 'EGM96' (Earth), 'LP100K' (Moon), 'LP165P'
%            (Moon), 'GMM2B' (Mars), 'Custom', or 'EIGENGL04C' (EARTH). The 
%            default is 'EGM2008'. 
%   DEGREE   :a scalar value specifying the degree and order of the      
%            harmonic gravity model. For 'EGM2008', the maximum degree and
%            order is 2159 and the default degree and order is 120.   For
%            'EGM96', the maximum degree and order is 360 and the default
%            degree and order is 70.  For 'LP100K', the maximum degree and
%            order is 100 and the default degree and order is 60.  For
%            'LP165P', the maximum degree and order is 165 and the default
%            degree and order is 60.  For 'GMM2B', the maximum degree and
%            order is 80 and the default degree and order is 60.  For
%            'Custom', the default degree and order is the maximum degree. 
%            For 'EIGENGL04C', the maximum degree and order is 360 and the 
%            default degree and order is 70.
%   DATAFILE :a file containing the planetary gravitational parameter,
%            planet equatorial radius, maximum degree, and normalized 
%            spherical harmonic coefficient matrices.
%   DFREADER :a function handle to an MATLAB(R) function which reads
%            DATAFILE.  The MATLAB function must output planetary
%            gravitational parameter in meters cubed per second squared,
%            planet equatorial radius in meters, maximum degree, and the
%            normalized spherical harmonic coefficient matrices, C and S.
%   ACTION   :a string to determine action for out of range input. Specify
%            if out of range input invokes a 'Warning', 'Error', or no
%            action ('None'). The default is 'Warning'.
%
%   Output calculated for the spherical harmonic gravity includes:
%   GX     :an array of M gravity values in the x-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared.
%   GY     :an array of M gravity values in the y-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared. 
%   GZ     :an array of M gravity values in the z-axis of the
%          Planet-Centered Planet-Fixed coordinates in meters per second
%          squared. 
%
%   Limitations:                                                           
%
%   This function has the limitations of excluding the centrifugal effects
%   of planetary rotation, and the effects of a precessing reference frame.
%
%   Spherical harmonic gravity model is valid for radial positions greater
%   than the planet equatorial radius.  Using it near or at the planetary
%   surface can probably be done with negligible error.  The spherical
%   harmonic gravity model is not valid for radial positions less than
%   planetary surface. 
%
%   Examples:                                                              
%
%   Calculate the gravity in the x-axis at the equator on the surface of
%   Earth, using the 120 degree model of EGM2008 with warning actions:
%       gx = gravitysphericalharmonic( [-6378.1363e3 0 0] ) 
%
%   Calculate the gravity at 25000 meters over the south pole of Earth using
%   the 70 degree model of EGM96 with error actions: 
%       [gx, gy, gz] = gravitysphericalharmonic( [0 0 -6381.751e3], 'EGM96', 'Error' )   
%
%   Calculate the gravity at 15000 meters over the equator and 11000 meters
%   over the north pole using a 30th order GMM2B Mars model with warning
%   actions:
%       p  = [2412.648e3 -2412.648e3 0; 0 0 3376.2e3]
%       [gx, gy, gz] = gravitysphericalharmonic( p, 'GMM2B', 30, 'Warning' )   
%
%   Calculate the gravity at 15000 meters over the equator and 11000 meters
%   over the north pole using a 60th degree custom planetary model with no
%   actions:  
%       p       = [2412.648e3 -2412.648e3 0; 0 0 3376e3]
%       [gx, gy, gz] = gravitysphericalharmonic( p, 'custom', 60, ...
%                       {'GMM2BC80_SHA.txt' @astReadSHAFile}, 'None' )
%
%   See also GRAVITYWGS84, GRAVITYCENTRIFUGAL, GRAVITYZONAL, GEOIDEGM96

%   Copyright 2009-2018 The MathWorks, Inc.

%   References:  
%   [1] Vallado, D. A., "Fundamentals of Astrodynamics and Applications",
%       McGraw-Hill, New York, 1997.  
%   [2] NIMA TR8350.2: "Department of Defense World Geodetic System 1984,
%       Its Definition and Relationship with Local Geodetic Systems." 
%   [3] Konopliv, A. S., S. W. Asmar, E. Carranza, W. L. Sjogen, D. N.
%       Yuan., "Recent Gravity Models as a Result of the Lunar Prospector
%       Mission", Icarus, Vol. 150, no. 1, pp 1-18, 2001.                    
%   [4] Lemoine, F. G., D. E. Smith, D.D. Rowlands, M.T. Zuber, G. A.
%       Neumann, and D. S. Chinn, "An improved solution of the gravity
%       field of Mars (GMM-2B) from Mars Global Surveyor", J. Geophys. Res.,
%       Vol. 106, No. E10, pp 23359-23376, October 25, 2001.   
%   [5] Kenyon S., J. Factor, N. Pavlis, and S. Holmes, "Towards the Next
%       Earth Gravitational Model", Society of Exploration Geophysicists
%       77th Annual Meeting, San Antonio, Texas, September 23-28, 2007.
%   [6] Pavlis, N.K., S.A. Holmes, S.C. Kenyon, and J.K. Factor, "An Earth
%       Gravitational Model to Degree 2160: EGM2008", presented at the 2008
%       General Assembly of the European Geosciences Union, Vienna,
%       Austria, April 13-18, 2008. 
%   [7] Grueber, T., and A. Kohl, "Validation of the EGM2008 Gravity Field
%       with GPS-Leveling and Oceanographic Analyses", presented at the IAG
%       International Symposium on Gravity, Geoid & Earth Observation 2008,
%       Chania, Greece, June 23-27, 2008.
%   [8] Forste, C., Flechtner, F., Schmidt, R., Konig, R., Meyer, U.,
%       Stubenvoll, R., Rothacher, M., Barthelmes, F., Neumayer, H.,
%       Biancale, R., Bruinsma, S., Lemoine, J.M., Loyer, S., "A Mean
%       Global Gravity Field Model From the Combination of Satellite
%       Mission and Altimetry/Gravmetry Surface Data - EIGEN-GL04C",
%       Geophysical Research Abstracts, Vol. 8, 03462, 2006.
%       http://icgem.gfz-potsdam.de/ICGEM/
%   [9] Hill, K. A. (2007). Autonomous Navigation in Libration Point Orbits. 
%       Doctoral dissertation, University of Colorado, Boulder.
%       http://ccar.colorado.edu/geryon/papers/Misc/Hill_thesis.pdf     
%  [10] Gottlieb, R. G., "Fast Gravity, Gravity Partials, Normalized Gravity, 
%       Gravity Gradient Torque and Magnetic Field: Derivation, Code and Data," 
%       Technical Report NASA Contractor Report 188243, NASA Lyndon B. Johnson 
%       Space Center, Houston, TX, February 1993.
%  [11] Colombo, Oscar L., "Numerical Methods for Harmonic Analysis on the
%       Sphere", Reports of the department of Geodetic Science, Report No.
%       310, The Ohio State University, Columbus, OH, March 1981.
%  [12] Colombo, Oscar L., "The Global Mapping of Gravity with Two
%       Satellites", Netherlands Geodetic Commission, Vol. 7 No. 3, Delft,
%       The Netherlands, 1984.
%  [13] Jones, Brandon A. (2010). Efficient Models for the Evaluation and
%       Estimation of the Gravity Field. Doctoral dissertation, University
%       of Colorado, Boulder.
%       http://ccar.colorado.edu/geryon/papers/Misc/bajones_phd.pdf

narginchk(1, 5);

checkinputs();

% set default values
model  = 'EGM2008';
action = 'warning';

switch nargin                                                              
    case 2
        % set degree or model
        if ~(isnumeric( varargin{1} ) && isreal( varargin{1} ))
            if ~ischar( varargin{1} ) && ~isstring( varargin{1} )
                error(message('aero:gravitysphericalharmonic:inputTypeVar'));
            else
                if strcmpi( varargin{1}, 'custom')
                    narginchk(4, 5);
                else
                    % assign model or action
                    modeloraction( varargin{1} );
                end
            end
        else
            maxdeg = varargin{1};
        end
    case 3
        if ~ischar( varargin{2} ) && ~isstring( varargin{2} )
            % not setting action
            if (~ischar( varargin{1} ) && ~isstring( varargin{1} ))|| ...
              ~(isnumeric( varargin{2} ) && isreal( varargin{2} ))
                error(message('aero:gravitysphericalharmonic:inputTypeNoAction3'));
            end
            if strcmpi( varargin{1}, 'custom')
                narginchk(4, 5);
            end
            % set model and degree
            checkmodel( varargin{1} );
            maxdeg = varargin{2};    
        else
            if ~(isnumeric( varargin{1} ) && isreal( varargin{1} ))
                if ~ischar( varargin{1} ) && ~isstring( varargin{1} )
                    error(message('aero:gravitysphericalharmonic:inputTypeVar3'));
                end
                if strcmpi( varargin{1}, 'custom')
                    narginchk(4, 5);
                end
                % set model and action
                checkmodel( varargin{1} );
                checkaction( varargin{2} );              
            else
            % set degree and action
            maxdeg = varargin{1};
            checkaction( varargin{2} );
            end
        end
    case 4
        if (~ischar( varargin{1} ) && ~isstring( varargin{1} )) || ...
                ~(isnumeric( varargin{2} ) && isreal( varargin{2} ))             
            error(message('aero:gravitysphericalharmonic:inputTypeVar4'));
        else
            if ischar( varargin{3} ) || isstring( varargin{3} )
                % set model, degree and action
                checkmodel( varargin{1} );
                maxdeg = varargin{2};
                checkaction( varargin{3} );
            elseif iscell( varargin{3} ) && ...
                    (ischar(varargin{3}{1}) || isstring( varargin{3}{1} )) && ...
                    isa(varargin{3}{2}, 'function_handle') && ...
                    strcmpi( varargin{1}, 'custom')
                % set model, degree, data file and data file reader
                checkmodel( varargin{1} );
                maxdeg = varargin{2};
                datafile = varargin{3}{1};                                 
                dfreader = varargin{3}{2};
            else
                if strcmpi( varargin{1}, 'custom')
                    error(message('aero:gravitysphericalharmonic:inputTypeVar4Data'));
                else
                    error(message('aero:gravitysphericalharmonic:inputTypeVar4Action'));
                end
            end
        end
    case 5
        if (ischar( varargin{1} ) || isstring( varargin{1} )) &&...
                ~strcmpi( varargin{1}, 'custom')
            % This option is only for 'custom'
            error(message('aero:gravitysphericalharmonic:wrongModel'));
        end
        if (~ischar( varargin{1} ) && ~isstring( varargin{1} )) || ...
           ~(isnumeric( varargin{2} ) && isreal( varargin{2} )) || ...
           ~(iscell( varargin{3} ) && ...
           (ischar(varargin{3}{1}) || isstring( varargin{3}{1} )) && ...
             isa(varargin{3}{2}, 'function_handle')) || ...
           (~ischar( varargin{4} ) && ~isstring( varargin{4} ))
            error(message('aero:gravitysphericalharmonic:inputTypeVar5'));
        end
        % set model, degree, data file, data file reader, and action
        checkmodel( varargin{1} );
        maxdeg   = varargin{2};
        datafile = varargin{3}{1};                                         
        dfreader = varargin{3}{2}; 
        checkaction( varargin{4} );     
end
           
switch lower( model )                                                      
    case 'custom'
        if ~exist( datafile, 'file' )
            error(message('aero:gravitysphericalharmonic:noDataFile', datafile))
        end
        try
            [GM, Re, degree, C, S] = dfreader( datafile );
        catch MECustomFile
            try
                % try loading a mat-file
                dfreader( datafile );  % [GM, Re, degree, C, S] 
            catch MECustomFileLoad
                throwAsCaller(MECustomFile)
            end                                                            
            % check for the existence of correct variables in mat-file
            if ~exist( 'GM', 'var' ) 
                error(message('aero:gravitysphericalharmonic:noGM', datafile))
            end
            if ~exist( 'Re', 'var' ) 
                 error(message('aero:gravitysphericalharmonic:noRe', datafile))
            end
            if ~exist( 'degree', 'var' ) 
                 error(message('aero:gravitysphericalharmonic:noDegree', datafile))
            end
            if ~exist( 'C', 'var' ) 
                 error(message('aero:gravitysphericalharmonic:noC', datafile))
            end
            if ~exist( 'S', 'var' ) 
                 error(message('aero:gravitysphericalharmonic:noS', datafile))
            end
        end
        
        % check data type and sizes of custom data
        if ~isscalar( GM ) || ~isscalar( Re ) || ~isscalar( degree ) || ~isnumeric( GM ) || ~isnumeric( Re ) || ~isnumeric( degree ) 
            error(message('aero:gravitysphericalharmonic:notScalar'))
        end
        if  (~isnumeric(C) || ~isnumeric(S)) || (~all( size(C) == [ degree+1 degree+1 ] ) || ~all( size(S) == [ degree+1 degree+1 ] ))
            error(message('aero:gravitysphericalharmonic:wrongMatrix'))
        end
        
        % The recommended default degree for custom is unknown, so use the
        % maximum degree of custom unless maxdeg is input
        checkmaxdeg( degree, degree )                                                   
    case 'egm2008'                                                         
        % Earth
        % Load normalized coefficients and planetary constants
        load('aeroegm2008.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 120 ) %#ok<NODEF>
    case 'egm96'
        % Earth
        % Load normalized coefficients and planetary constants
        load('aeroegm96.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 70 ) %#ok<NODEF>
    case 'lp100k' 
        % Moon
        % Load normalized coefficients and planetary constants
        load('aerolp100k.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 60 ) %#ok<NODEF>
    case 'lp165p'
        % Moon
        % Load normalized coefficients and planetary constants
        load('aerolp165p.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 60 ) %#ok<NODEF>
    case 'gmm2b'
        % Mars
        % Load normalized coefficients and planetary constants
        load('aerogmm2b.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 60 ) %#ok<NODEF>
    case 'eigengl04c'
        % Earth
        % Load normalized coefficients and planetary constants
        load('aeroeigengl04c.mat') % [GM, Re, degree, C, S] 
        checkmaxdeg( degree, 70 ) %#ok<NODEF>     
end

% Compute geocentric radius
r = sqrt( sum( p.^2, 2 ));
                                                                           
% Check if geocentric radius is less than equatorial (reference) radius
if r < Re
    switch action
        case 'none'
            % no message
        case 'warning'
            warning(message('aero:gravitysphericalharmonic:lessThanEquatorialRadius', sprintf( '%g', Re )));
        case 'error'
            error(message('aero:gravitysphericalharmonic:lessThanEquatorialRadius', sprintf( '%g', Re )));
        otherwise
            error(message('aero:gravitysphericalharmonic:unknownActionRadius'));
    end
end

% Compute geocentric latitude
phic = asin( p(:,3)./ r );

% Compute lambda                                                           
lambda = atan2( p(:,2), p(:,1) );

smlambda = zeros( size(p,1), maxdeg + 1 );
cmlambda = zeros( size(p,1), maxdeg + 1 );

slambda = sin(lambda);
clambda = cos(lambda);
smlambda(:,1) = 0;
cmlambda(:,1) = 1;
smlambda(:,2) = slambda;
cmlambda(:,2) = clambda;

for m=3:maxdeg+1
    smlambda(:,m) = 2.0.*clambda.*smlambda(:, m-1) - smlambda(:, m-2);
    cmlambda(:,m) = 2.0.*clambda.*cmlambda(:, m-1) - cmlambda(:, m-2);
end

% Compute normalized associated legendre polynomials
[P,scaleFactor] = loc_gravLegendre( phic, maxdeg );

% Compute gravity in ECEF coordinates
[gx,gy,gz] = loc_gravityPCPF( p, maxdeg, P, C( 1:maxdeg+1, 1:maxdeg+1 ), ...
                                  S( 1:maxdeg+1, 1:maxdeg+1 ), smlambda, ...
                                  cmlambda, GM, Re, r,scaleFactor );

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkinputs( )
        if ~isnumeric( p )
            error(message('aero:gravitysphericalharmonic:notNumeric'));
        end
        
        if (size( p, 2) ~= 3)
            error(message('aero:gravitysphericalharmonic:wrongDimension'));
        end
   end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkmaxdeg( upperlimit, default )
        if exist('maxdeg','var') && maxdeg > upperlimit
            switch action
                case 'none'
                    % no passing of maxdeg messages
                case 'warning'
                    warning(message('aero:gravitysphericalharmonic:exceedMaxDegWarn'));
                case 'error'
                    error(message('aero:gravitysphericalharmonic:exceedMaxDegError'));
                otherwise
                    error(message('aero:gravitysphericalharmonic:unknownActionMaxDeg'));
            end
            maxdeg = upperlimit;
        else
            if ~exist('maxdeg','var')
            % maxdeg was not set
            maxdeg = default;            
            end
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkmodel( str )
        switch lower( str )
            case { 'egm2008', 'egm96', 'lp100k', 'lp165p', 'gmm2b', 'custom', 'eigengl04c' }
                model = lower( str );
            otherwise
                error(message('aero:gravitysphericalharmonic:unknownModel'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function checkaction( str )
        switch lower( str )
            case { 'error', 'warning', 'none' }
                action = lower( str );
            otherwise
                error(message('aero:gravitysphericalharmonic:unknownAction'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    function modeloraction( str )
        switch lower( str )
            case { 'egm2008', 'egm96', 'lp100k', 'lp165p', 'gmm2b', 'eigengl04c' }
                model = lower( str );
            case { 'error', 'warning', 'none' }
                action = lower( str );
            otherwise
                error(message('aero:gravitysphericalharmonic:unknownString'));
        end
    end
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [P,scaleFactor] = loc_gravLegendre( phi, maxdeg )
% loc_GRAVLEGENDRE internal function computing normalized associated 
% legendre polynomials, P, via recursion relations for spherical harmonic
% gravity 

P = zeros(maxdeg+3, maxdeg+3, length(phi));
scaleFactor = zeros(maxdeg+3, maxdeg+3, length(phi));
cphi = cos(pi/2-phi);
sphi = sin(pi/2-phi);

% force numerically zero values to be exactly zero
cphi(abs(cphi)<=eps) = 0;
sphi(abs(sphi)<=eps) = 0;
 
% Seeds for recursion formula
P(1,1,:) = 1;            % n = 0, m = 0;
P(2,1,:) = sqrt(3)*cphi; % n = 1, m = 0;
scaleFactor(1,1,:) = 0;
scaleFactor(2,1,:) = 1;
P(2,2,:) = sqrt(3)*sphi; % n = 1, m = 1;
scaleFactor(2,2,:) = 0;

for n = 2:maxdeg+2
    k = n + 1;
    for m = 0:n
        p = m + 1;
        % Compute normalized associated legendre polynomials, P, via recursion relations 
        % Scale Factor needed for normalization of dUdphi partial derivative
                
        if (n == m)           
            P(k,k,:) = sqrt(2*n+1)/sqrt(2*n)*sphi.*reshape(P(k-1,k-1,:),size(phi));
            scaleFactor(k,k,:) = 0;
        elseif (m == 0)
            P(k,p,:) = (sqrt(2*n+1)/n)*(sqrt(2*n-1)*cphi.*reshape(P(k-1,p,:),size(phi)) - (n-1)/sqrt(2*n-3)*reshape(P(k-2,p,:),size(phi)));
            scaleFactor(k,p,:) = sqrt( (n+1)*(n)/2);
        else
            P(k,p,:) = sqrt(2*n+1)/(sqrt(n+m)*sqrt(n-m))*(sqrt(2*n-1)*cphi.*reshape(P(k-1,p,:),size(phi)) - sqrt(n+m-1)*sqrt(n-m-1)/sqrt(2*n-3)*reshape(P(k-2,p,:),size(phi)));
            scaleFactor(k,p,:) = sqrt( (n+m+1)*(n-m));
        end
    end
end
end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
function [gx,gy,gz] = loc_gravityPCPF(p,maxdeg,P,C,S,smlambda,cmlambda,GM,Re,r,scaleFactor)
% loc_GRAVITYPCPF internal function computing gravity in planet-centered
% planet-fixed (PCEF) coordinates using PCPF position, desired
% degree/order, normalized associated legendre polynomials, normalized
% spherical harmonic coefficients, trigonometric functions of geocentric
% latitude and longitude, planetary constants, and radius to center of
% planet. Units are MKS.

rRatio   = Re./r;
rRatio_n = rRatio;

% initialize summation of gravity in radial coordinates
dUdrSumN      = 1;
dUdphiSumN    = 0;
dUdlambdaSumN = 0;

% summation of gravity in radial coordinates
for n = 2:maxdeg
    k = n+1;
    rRatio_n      = rRatio_n.*rRatio;
    dUdrSumM      = 0;
    dUdphiSumM    = 0;
    dUdlambdaSumM = 0;
    for m = 0:n
        j = m+1;
        dUdrSumM      = dUdrSumM + reshape(P(k,j,:),size(r)).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j)); 
        dUdphiSumM    = dUdphiSumM + ( (reshape(P(k,j+1,:),size(r)).*reshape(scaleFactor(k,j,:),size(r))) - p(:,3)./(sqrt(p(:,1).^2 + p(:,2).^2)).*m.*reshape(P(k,j,:),size(r))).*(C(k,j).*cmlambda(:,j) + S(k,j).*smlambda(:,j)); 
        dUdlambdaSumM = dUdlambdaSumM + m*reshape(P(k,j,:), size(r)).*(S(k,j).*cmlambda(:,j) - C(k,j).*smlambda(:,j));
    end
    dUdrSumN      = dUdrSumN      + dUdrSumM.*rRatio_n.*k;
    dUdphiSumN    = dUdphiSumN    + dUdphiSumM.*rRatio_n;
    dUdlambdaSumN = dUdlambdaSumN + dUdlambdaSumM.*rRatio_n;
end

% gravity in spherical coordinates
dUdr      = -GM./(r.*r).*dUdrSumN;
dUdphi    =  GM./r.*dUdphiSumN;
dUdlambda =  GM./r.*dUdlambdaSumN;

% gravity in ECEF coordinates
gx = ((1./r).*dUdr - (p(:,3)./(r.*r.*sqrt(p(:,1).^2 + p(:,2).^2))).*dUdphi).*p(:,1) ...
      - (dUdlambda./(p(:,1).^2 + p(:,2).^2)).*p(:,2); 
gy = ((1./r).*dUdr - (p(:,3)./(r.*r.*sqrt(p(:,1).^2 + p(:,2).^2))).*dUdphi).*p(:,2) ...
      + (dUdlambda./(p(:,1).^2 + p(:,2).^2)).*p(:,1); 
gz = (1./r).*dUdr.*p(:,3) + ((sqrt(p(:,1).^2 + p(:,2).^2))./(r.*r)).*dUdphi;

% special case for poles
atPole = abs(atan2(p(:,3),sqrt(p(:,1).^2 + p(:,2).^2)))==pi/2;
if any(atPole)
    gx(atPole) = 0;
    gy(atPole) = 0;
    gz(atPole) = (1./r(atPole)).*dUdr(atPole).*p((atPole),3);
end

end
