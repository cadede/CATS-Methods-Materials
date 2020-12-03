% Aerospace Toolbox
% Version 3.1 (R2019a) 23-Nov-2018
%
% Axes Transformations.
%   angle2dcm          - Create direction cosine matrix from rotation angles.
%   angle2quat         - Convert rotation angles to quaternion.
%   angle2rod          - Convert rotation angles to Rodrigues vector.
%   dcm2alphabeta      - Convert direction cosine matrix to angle of attack 
%                        and sideslip angle.
%   dcm2angle          - Create rotation angles from direction cosine matrix.
%   dcm2latlon         - Convert direction cosine matrix to geodetic 
%                        latitude and longitude.
%   dcm2quat           - Convert direction cosine matrix to quaternion.
%   dcm2rod            - Convert direction cosine matrix to Rodrigues vector.
%   dcmbody2wind       - Convert angle of attack and sideslip angle to 
%                        direction cosine matrix.
%   dcmecef2ned        - Convert geodetic latitude and longitude to 
%                        direction cosine matrix.
%   dcmeci2ecef        - Convert from Earth-centered inertial (ECI) to Earth-
%                        centered Earth-fixed (ECEF) direction cosine matrix.
%   ecef2lla           - Convert Earth-centered Earth-fixed (ECEF) 
%                        coordinates to geodetic coordinates.
%   eci2aer            - Convert Earth-centered inertial (ECI) coordinates to 
%                        azimuth elevation range coordinates.
%   eci2lla            - Convert Earth-centered inertial (ECI) coordinates to 
%                        geodetic coordinates.
%   flat2lla           - Estimate geodetic latitude, longitude, and altitude
%                        from flat Earth position. 
%   geod2geoc          - Convert geodetic latitude to geocentric latitude.
%   geoc2geod          - Convert geocentric latitude to geodetic latitude.
%   lla2ecef           - Convert geodetic coordinates to Earth-centered 
%                        Earth-fixed (ECEF) coordinates.
%   lla2eci            - Convert geodetic coordinates to Earth-centered
%                        Inertial (ECI) coordinates.
%   lla2flat           - Estimate flat Earth position from geodetic latitude, 
%                        longitude, and altitude. 
%   quat2angle         - Convert quaternion to rotation angles.
%   quat2dcm           - Convert quaternion to direction cosine matrix.
%   quat2rod           - Convert quaternion to Rodrigues vector.
%   rod2angle          - Convert Rodrigues vector to rotation angles.
%   rod2dcm            - Convert Rodrigues vector to direction cosine matrix.
%   rod2quat           - Convert Rodrigues vector to quaternion.
%      
% Environment.
%   atmoscira                - Use COSPAR International Reference Atmosphere 1986.
%   atmoscoesa               - Use 1976 COESA atmosphere.
%   atmoshwm                 - Use Horizontal Wind Model.
%   atmosisa                 - Use International Standard Atmosphere Model.
%   atmoslapse               - Use Lapse Rate Atmosphere Model.
%   atmosnonstd              - Use climatic data from MIL-STD-210 or MIL-HDBK-310.
%   atmosnrlmsise00          - Use NRLMSISE-00 atmosphere model.
%   atmospalt                - Calculate pressure altitude based on ambient 
%                              pressure.
%   geoidheight              - Implement a geopotential model to calculate geoid height
%   gravitycentrifugal       - Implement a centrifugal effect of planetary gravity.
%   gravitysphericalharmonic - Implement a spherical harmonic representation
%                              of planetary gravity. 
%   gravitywgs84             - Implement the 1984 World Geodetic System (WGS84)
%                              representation of Earth's gravity.
%   gravityzonal             - Implement a zonal harmonic representation of 
%                              planetary gravity.
%   igrfmagm                 - Use International Geomagnetic Reference Field.
%   wrldmagm                 - Use World Magnetic Model.
%
% File Reading.
%   datcomimport       - Bring USAF Digital DATCOM file into MATLAB.
%
% Flight Instruments.
%   uiaeroairspeed      - Create airspeed indicator instrument in MATLAB figure.
%   uiaeroaltimeter     - Create altimeter instrument in MATLAB figure.
%   uiaeroclimb         - Create climb rate indicator instrument in MATLAB figure.
%   uiaeroegt           - Create exhaust gas temperature (EGT) indicator instrument in MATLAB figure.
%   uiaeroheading       - Create heading indicator instrument in MATLAB figure.
%   uiaerohorizon       - Create artificial horizon indicator instrument in MATLAB figure.
%   uiaerorpm           - Create revolutions per minute (RPM) indicator instrument in MATLAB figure.
%   uiaeroturn          - Create turn coordinator indicator instrument in MATLAB figure.
%
% Animation.
%   Aero.Animation                 - Construct animation object.
%   Aero.Body                      - Construct body object for use with 
%                                    animation object.
%   Aero.Camera                    - Construct camera object for use with 
%                                    animation object.
%   Aero.FlightGearAnimation       - Construct FlightGear animation object.
%   fganimation                    - Construct FlightGear animation object.
%   Aero.Geometry                  - Construct 3-D geometry for use with 
%                                    animation object.
%   Aero.Node                      - Construct node object for use with
%                                    virtual reality animation object.
%   Aero.Viewpoint                 - Construct viewpoint object for use with
%                                    virtual reality animation object.
%   Aero.VirtualRealityAnimation   - Construct virtual reality animation
%                                    object.
%
% Flight Parameters.
%   airspeed          - Compute airspeed from velocity.
%   alphabeta         - Compute incidence and sideslip angles.
%   correctairspeed   - Calculate equivalent airspeed (EAS), calibrated 
%                       airspeed (CAS) or true airspeed (TAS) from one of 
%                       the other two airspeeds.
%   dpressure         - Compute dynamic pressure using velocity and density.
%   geocradius        - Estimate radius of ellipsoid planet at geocentric 
%                       latitude.
%   machnumber        - Compute Mach number using velocity and speed of sound.
%   rrdelta           - Compute relative pressure ratio.
%   rrsigma           - Compute relative density ratio.
%   rrtheta           - Compute relative temperature ratio.
%
% Quaternion Math.
%   quatconj          - Calculate the conjugate of a quaternion.
%   quatdivide        - Divide a quaternion by another quaternion.
%   quatexp           - Calculate the exponential of a quaternion.
%   quatinv           - Calculate the inverse of a quaternion.
%   quatinterp        - Interpolate between two quaternions.
%   quatlog           - Calculate the natural logarithm of a quaternion.
%   quatmod           - Calculate the modulus of a quaternion.
%   quatmultiply      - Calculate the product of two quaternions.
%   quatnorm          - Calculate the norm of a quaternion.
%   quatnormalize     - Normalize a quaternion.
%   quatpower         - Calculate the power of a quaternion.
%   quatrotate        - Rotate a vector by a quaternion.
%
% Time.
%   decyear           - Calculate decimal year.
%   deltaUT1          - Compute the difference between Universal Time (UT1) and 
%                       Coordinated Universal time (UTC) UT1-UTC.
%   juliandate        - Calculate Julian date.
%   leapyear          - Determine leap year.
%   mjuliandate       - Calculate modified Julian date.
%   tdbjuliandate     - Calculate Julian date for Barycentric Dynamical time.
%
% Unit Conversions.
%   convacc           - Convert from acceleration units to desired 
%                       acceleration units.
%   convang           - Convert from angle units to desired angle units.
%   convangacc        - Convert from angular acceleration units to desired 
%                       angular acceleration units.
%   convangvel        - Convert from angular velocity units to desired 
%                       angular velocity units.
%   convdensity       - Convert from density units to desired density units.
%   convforce         - Convert from force units to desired force units.
%   convlength        - Convert from length units to desired length units.
%   convmass          - Convert from mass units to desired mass units.
%   convpres          - Convert from pressure units to desired pressure units.
%   convtemp          - Convert from temperature units to desired 
%                       temperature units.
%   convvel           - Convert from velocity units to desired velocity units.
%
% Gas Dynamics.
%   flowisentropic    - Isentropic flow ratios.
%   flownormalshock   - Normal shock relations.
%   flowprandtlmeyer  - Prandtl-Meyer functions for expansion waves.
%   flowfanno         - Fanno line flow relations for 1-D fluid flow with friction.
%   flowrayleigh      - Rayleigh line relations for 1-D flow with heat transfer.
%
% Celestial Phenomena.
%   deltaCIP          - Calculate Celestial Intermediate Pole (CIP) adjustment 
%                       location.
%   earthNutation     - Implement JPL development ephemeris for Earth 
%                       nutation longitude and obliquity angles.
%   moonLibration     - Implement JPL development ephemeris for Moon libration 
%                       Euler angles.
%   planetEphemeris   - Implement JPL development ephemeris for specified 
%                       bodies of the Solar System.
%   polarMotion       - Calculate polar motion.
%
% Utilities.
%   aeroDataPackage   - Utility command to launch the application to obtain or 
%                       install data required by some Aerospace Toolbox
%                       functions.
%   aeroReadIERSData  - Utility command to import Earth Orientation Data from
%                       the U.S. Naval Observatory.

% Copyright 1990-2018 The MathWorks, Inc.
