Goldbogen Lab
Hopkins Marine Station
Stanford University

Tags-
2-4: Platypus style, dual side-by-side cams, 40 Hz accelerometry, some have paddlewheel data, usually ~3 hours video
5: minion style, single cam, 40 Hz accelerometry
6-12: froback style, dual front-back cams, 790 Hz or 400 Hz accelerometery, 50 Hz other data, GPS, up to 7 hours video
20-25: Platypus style, dual side-by-side cams, 400 Hz accelerometery, 50 Hz other data, GPS, up to 7 hours video
26: Stereo, forward looking cams.  400 Hz accelerometery, 50 Hz other data, GPS, up to 7 hours video.  No y-axis accelerometry
30-32 36, 37: 2k Cameras. 400 Hz accelerometery, 50 Hz other data, GPS, up to 10 hours video
40-44: Wireless cameras, with CATS specific cameras.  400 Hz accelerometery, 50 Hz other data, GPS, up to 8 hours video

prh files are generally downsampled to 10 Hz.  Can resample files using CATS toolbox.  

raw tag data has both local and UTC time.  Local time in raw data is the time on the computer for which the tags were downloaded. prh files are in whale local time.  Adjustments from tag time are in the INFO.timedif variable (in hours)

Notes: 
All axes are in a North-East-Down (x,y,z) orientation. (Dtag3- accel NEU, Mag NEU)
Facing up towards gravity gives positive acceleration.  (Dtag3- same)
Euler angle orientations are calculated in heading, pitch, roll order and are counter-clockwise rotations when facing the plane of rotation from the positive side of the perpendicular axis.
This means positive pitch is up, positive roll is a roll to the right, and heading is congruent to compass orientation.  (Dtag3- positive roll is to the left)

At- Acceleration data for each axis (x,y,z) in the tag's reference frame (g's)
Aw- Acceleration data for each axis rotated to the whale's reference frame
camon- logical index of where in the data the camera was on
DN- index of datenumbers corresponding to each data point (Local Time)
DV- index of datevectors corresponding to each data point (Local Time)
flownoise- DB value of flownoise in 66-94 Hz band
fs- sample rate of the data, in Hz
geoPtrack- [x,y,z] in meters of georeferenced pseudotrack
GPS- An index of Lat, Long for each position that is known.  This is mostly nans, except where a GPS position is known.  This file is trimmed to include only presumably accurate posisitions that are not obviously wrong (i.e. they fit an approximate pseudotrack).  GPS(1,:) is an approx. GPS of the area, used for calculating inc & dec & b.
GPSerr- for tagged derived GPS points (not for manual GPS points from focal follows), an error estimate of the std in m.
Gt- Gyroscope data (radians/sec) around each axis (x,y,z) in the tag's reference frame
Gw- Gyroscope data rotated to the whale's reference frame
GwUB- Gyroscope data in the whale's reference frame with the bias removed for each sample.  
head- heading of the whale (radians)
headgy- heading of the whale (same as above except when jerk was high, then used the gyros to calculate the orientation)
INFO- a structure with some metadata of the prh file creation.  Includes:
INFO.cal- first 8 values are bench calibrations, application: (At-aconst)*acal.  Also accounts for axis orientations in raw data that differ from NED.  CAL.Mag3d and p3d are calibration constants from spherical_cal and fix_pressure from Mark Johnson/Stacy DeRuiter tag_tools workshop
INFO.calperiod- datenumbers of the surfacings and dive used to calculate each W
INFO.calperiodI- indices corresponding to above
INFO.notes- information on the prh file creation
INFO.tagprh- the pitch, roll and heading (in radians) of the tag from the whale frame
INFO.tagslip- 4 indices of when the tag slip.  Wchange- indices of where a different W was calculated.  Speedperiods- indices of when a different speed calibration curve was used. Worigcalperiods were the indices of the first manually identified tag slips (probably not useful).
INFO.W- cell with an orientation matrix calculating the rotation of the tag onto the whale (i.e. Aw = At*W for each calibration period)
INFO.usePaddles- true if a paddle wheel with good speed calculations was present.
JigRMS- the accelerometer motion in DB re 1 g for all axes, used to calculate speed.JJ.
Light- light levels (lux in tags 2-4, uncalibrated raw in newer tags).  Maximum value recorded is 3500 in tags 2-4
LightIR- infrared light levels (uncalibrated units).
Mt- Magnetometer data for each axis (x,y,z) in the tag's reference frame (uT).
Mw- Magnetometer data for each axis in the whale's reference frame
p- pressure (depth) in meters
Ptrack- [x,y,z] in meters of the pseudotrack
Paddles- impulses/second of the paddlewheel (tags 2-4 only)
pitch- pitch of the whale (radians)
pitchgy- pitch of the whale (same as above except when jerk was high, then used the gyros to calculate the orientation)
roll- roll of the whale (radians)
rollgy- roll of the whale (same as above except when jerk was high, then used the gyros to calculate the orientation)
speed- a table of speed values from different calculations
speedstats- a structure containing info on the models used to create the speed table (see speedstats.info and Cade et al 2017).
speed.FN- speed from flow noise (only when the camera/hydrophone is on and enabled)
speed.SP- speed from change in depth over sine of pitch
speed.JJ- speed from the jiggling of the tag (Cade et al. 2017 JEB).  Probably best but be sure to exlude surfacings (spikes in speed at those points)
T- temperature
tagon- logical index indicating when the tag was on the whale
UTC- UTC offset (usually calculated when GPS was added).
viddeploy- the videos recording time when the tag was on the whale
vidDN- datenumbers of the start of the videos
vidDurs- length of the videos (seconds)
vidNam- names of the videos

