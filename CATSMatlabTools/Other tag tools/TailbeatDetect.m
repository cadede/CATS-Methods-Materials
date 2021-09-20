% Flukebeat Detector and Analyser
% William Gough 
% Version 1.0 (2021)
% Stanford University

% This function is designed to analyze a segment of inertial sensing data
% and uses a series of thresholds to determine when tailbeat fluking is occuring.
% Note: Please remove or interpolate all NaN values from your data and
% ensure that all data is time-synced (especially if using pressure/speed).

% The necessary inputs for this function are:
% Inertial Sensing Data Sample (accelerometer x-axis or gyroscope y-axis) (InertialData)
% Data Sampling Rate (fs)

% Optional inputs are (in order):
% Pressure (depth) data (Depth)
% Surface threshold that denotes the line (default is 2m) between a surface tailbeat and a deep tailbeat (SurfaceThresh).
% Speed data (Speed)
% Threshold defining the maximum allowable length (default is 100 indices i.e. 10 seconds) of a tailbeat (DurThresh).
% Threshold defining the minimum distance (default is 1.25 times the zero crossing value) a tailbeat's peaks must be away from the zero crossing value (HeightThresh).
% Threshold defining the minimum magnitude (default is 1) of a tailbeat (MagThresh).
% Threshold (default is 10) used to compare the area under the curve for the upstroke and downstroke of a tailbeat to ensure that it is not too top- or bottom-heavy (ClarityThresh).

function [FlukingOverlay Tailbeats] = TailbeatDetect(InertialData,fs,Depth,Pitch,Roll,Heading,SurfaceThresh,Speed,DurThresh,HeightThresh,MagThresh,ClarityThresh, Lowpass, Highpass) % need switch to tell if gyro or accl for InertialData

NaNs = find(isnan(InertialData)==1); % This finds all of the nan positions, making it easier to align the inertial sensing data with the depth and speed data, if using either.
InertialData(NaNs,:) = []; % If you have not removed NaNs from your inertial sensing data, this will do that for you.

% This segment tests to see which inputs have been provided and creates
% empty variables or default values for the rest.
if nargin < 2 || isempty(fs); fs = 10; disp ('Set Sample Rate To 10'); end % If no sample rate is give, will default to 10Hz.
if nargin < 3 || isempty(Depth); Depth = []; disp('No Depth Data Provided'); end % If no depth data is give, will provide empty variable.
if nargin < 4 || isempty(Pitch); Pitch = []; disp('No Pitch Data Provided'); end % If no pitch data is give, will provide empty variable.
if nargin < 5 || isempty(Roll); Roll = []; disp('No Roll Data Provided'); end % If no roll data is give, will provide empty variable.
if nargin < 6 || isempty(Heading); Heading = []; disp('No Heading Data Provided'); end % If no heading data is give, will provide empty variable.
if nargin < 7 || isempty(SurfaceThresh); SurfaceThresh = 2; end % Threshold that denotes surface vs. deep tailbeats.
if nargin < 8 || isempty(Speed); Speed = []; disp('No Speed Data Provided'); end % If no speed data is give, will provide empty variable.
if nargin < 9 || isempty(DurThresh); DurThresh = 100; disp('Using Default Threshold Values'); end % If no thresholds are provided, will use default values (10s for tailbeat maximum duration).
if nargin < 10 || isempty(HeightThresh); HeightThresh = 0; end % Height of tailbeat (minumum distance away from zero value) set to 0 as placeholder - will be changed later.
if nargin < 11 || isempty(MagThresh); MagThresh = 1; end % This threshold makes sure that the magnitude of the tailbeat is above a value - set to 1 as a default.
if nargin < 12 || isempty(ClarityThresh); ClarityThresh = 10; end % This threshold ensures that the upstroke and downstroke are within a given multiple of one another - set to 10 as a default.
if nargin < 13 || isempty(Lowpass); Lowpass = 0.5; end % This threshold sets the lowpass frequency and will filter out any high frequency noise from your dataset - set to 0.5 Hz as a default.
if nargin < 14 || isempty(Highpass); Highpass = 0.1; end % This threshold sets the highpass frequency and will filter out any low frequency body motions from your dataset - set to 0.1 Hz as a default.

if ~isempty(Depth); Depth(NaNs,:) = []; end
if ~isempty(Speed); Speed(NaNs,:) = []; Speed(isnan(Speed)) = 0; end

% This segment will create and run a low-pass filter and high-pass filter
% over the inertial sensing data (InertialData) and provide a filtered
% version of the same data called (InertialDataFilt).
N = 10; % order
Fpasslow = Lowpass; % lowpass filter passband frequency
Fstoplow = Lowpass+0.1; % lowpass filter stopband frequency
Wpasslow = 0.1; % passband weight
Wstoplow = 0.1; % stopband weight
denslow  = 20; % density factor
blow  = firpm(N, [0 Fpasslow Fstoplow fs/2]/(fs/2), [1 1 0 0], [Wpasslow Wstoplow], ...
       {denslow});
Hdlow = dfilt.dffir(blow);
InertialDataFilt = filtfilt(blow,1,InertialData); % lowpass filters the dataset

Fstophigh = Highpass-0.1;    % Highpass Filter Stopband Frequency
Fpasshigh = Highpass;   % Highpass Filter Passband Frequency
Dstophigh = 0.0001;          % Stopband Attenuation
Dpasshigh = 0.057501127785;  % Passband Ripple
denshigh  = 20;              % Density Factor
[N, Fo, Ao, W] = firpmord([Fstophigh, Fpasshigh]/(fs/2), [0 1], [Dstophigh, Dpasshigh]);
bhigh  = firpm(N, Fo, Ao, W, {denshigh});
Hdhigh = dfilt.dffir(bhigh);
InertialDataFilt = filtfilt(bhigh,1,InertialDataFilt);

% This segment will create and run a low-pass filter over the speed data,
% if it is present.
if ~isempty(Speed);
    N = 35; % order
    Fpasslow = 0.09; % lowpass filter passband frequency
    Fstoplow = 0.2; % lowpass filter stopband frequency
    Wpasslow = 0.4; % passband weight
    Wstoplow = 0.4; % stopband weight
    denslow  = 60; % density factor
    blow  = firpm(N, [0 Fpasslow Fstoplow fs/2]/(fs/2), [1 1 0 0], [Wpasslow Wstoplow], ...
           {denslow});
    Hdlow = dfilt.dffir(blow);
    SpeedFilt = filtfilt(blow,1,Speed); % lowpass filters the jiggle speed
else
    SpeedFilt = [];
end

%This segment will start set a "zero crossing" value and define tailbeats
%based on that value.
zeroavg = mean(InertialDataFilt); % This sets the zero crossing value for use throughout the inertial sensing data.
if HeightThresh == 0; HeightThresh = mean(zeroavg)*1.25; end % If the tailbeat height threshold was not given, reset it here to be 1.25 times the mean value of the filtered inertial data.

abovezeroavg = [];
flkStarts = [];
flkStarts(1) = 0;
lookahead = 3;

for a=1:length(InertialDataFilt); % This creates an array (abovezeroavg) for segments when the inertial sensing data is above the zero average (proxy for fluking - "zero crossing").
    if (InertialDataFilt(a)) >= zeroavg;
        abovezeroavg(a,1) = 1;
    else
        abovezeroavg(a,1) = 0;
    end
end
clear a;

for a=2:length(InertialDataFilt)-lookahead; % This creates an array with zero crossings (flkStarts) that will denote the beginnings positions of tailbeats
    if abovezeroavg(a-1) == 0 && abovezeroavg(a) == 1 && abovezeroavg(a+lookahead) == 1;
        flkStarts(a,1) = 1;
    else
        flkStarts(a,1) = 0;
    end
end
clear a;

flkStarts(end+1:end+3) = 0; % Adds zeros to the end of flkStarts to make sure it is the same length as the inertial sensing data array

flukebeatstart = [];
flukebeatend = [];
useableflukebeats = [];

for a=1:length(flkStarts); % This will create a list of starting and ending indices for each tailbeat found in flkStarts
    if flkStarts(a) == 0;
        continue
    else
        flukebeatstart(end+1,1) = a;
        for b=(a+1):length(flkStarts);
            if flkStarts(b) == 0;
                continue
            else
                flukebeatend(end+1,1) = b;
                break
            end
        end
    end
end
clear a b;

flukebeatend(end+1,1) = length(flkStarts);
useableflukebeats = [flukebeatstart flukebeatend];

% This segments tests each tailbeat against the given thresholds values to
% determine if it is a viable, symmetrical tailbeat that should be saved or
% if it should be removed from the dataset.
for a=1:length(useableflukebeats(:,1));
    if minus((useableflukebeats(a,2)),(useableflukebeats(a,1))) >= DurThresh; % This tests the length of the tailbeat against the duration threshold.
        useableflukebeats(a,3:4) = 0;
    else
        tempforint = InertialDataFilt((useableflukebeats(a,1)):(useableflukebeats(a,2)));
        useableflukebeats(a,3) = trapz(tempforint); % This runs a trapezoidal integration on the tailbeat.
        useableflukebeats(a,4) = trapz(abs(tempforint)); % This runs another trapezoidal integration on the tailbeat to determine the approximate magnitude of the flukebeat.
        [pks,locs] = findpeaks(InertialDataFilt(useableflukebeats(a,1):useableflukebeats(a,2)));
        [pksabs,locsabs] = findpeaks(abs(InertialDataFilt(useableflukebeats(a,1):useableflukebeats(a,2))));
        highpks = find(pksabs > (HeightThresh));
        if useableflukebeats(a,3) <= ClarityThresh && ... % This tests to make sure the tailbeat isn't too top-heavy.
           useableflukebeats(a,3) >= -ClarityThresh && ... % This tests to make sure the tailbeat isn't too bottom-heavy.
           useableflukebeats(a,4) >= MagThresh && ... % This tests to make sure the tailbeat magnitude is above a threshold.
           length(locs) == 1 && ... % This tests to make sure there is only one positive peak in the tailbeat.
           length(locsabs) == 2 && ... % This tests to make sure there are exactly two peaks (one positive, one negative) in a tailbeat.
           isempty(highpks) == 0; % tests to make sure the tailbeat is at least as high as a threshold in the positive direction.
           continue
        else
            useableflukebeats(a,3:4) = 0;
        end
        continue
    end
end
clear a;

indices = find(useableflukebeats(:,3) == 0); % This removes tailbeats that are not within desired boundaries.
useableflukebeats(indices,:) = [];
    
isfluking = zeros(length(InertialDataFilt),1);

for flukephase = useableflukebeats.' % This creates an array denoting when the whale is (1) and is not (0) actively fluking.
    isfluking(flukephase(1):flukephase(2)) = 1;
end

flukingnans = find(~isfluking);
flukingnums = find(isfluking);
isfluking(flukingnans) = NaN; % changes 0s into NaNs
isfluking(flukingnums) = zeroavg; % changes 1s into 0s
FlukingOverlay = [InertialDataFilt isfluking]; % This appends the array of active fluking (isfluking) onto the filtered inertial sensing data (InertialDataFilt).

% This segment adds calculate the oscillatory frequency of each flukebeat
% (useableflukebeats column 7) and the overall mean oscillatory frequency
% (dsf) in Hz.
useableflukebeats(:,5) = minus((useableflukebeats(:,2)),(useableflukebeats(:,1)));
useableflukebeats(:,6) = useableflukebeats(:,5) ./ 10;
useableflukebeats(:,7) = 1 ./ useableflukebeats(:,6);
dsf = mean(useableflukebeats(:,7));

% If depth and/or speed data are provided, this segment determines the mean
% depth, relation to surface threshold, and mean speed values for each
% tailbeat.
for a=1:length(useableflukebeats(:,1));
    if ~isempty(Depth);
        useableflukebeats(a,8) = mean(Depth(useableflukebeats(a,1):useableflukebeats(a,2)));
        if useableflukebeats(a,8) < SurfaceThresh;
            useableflukebeats(a,9) = 0;
        else
            useableflukebeats(a,9) = 1;
        end
    end
    if ~isempty(SpeedFilt);
        useableflukebeats(a,10) = mean(SpeedFilt(useableflukebeats(a,1):useableflukebeats(a,2)));
    end
end

% This segment displays the mean oscillatory frequency (dsf) and number of
% tailbeats for two sets of tailbeats: 1) all tailbeats and 2) non-surface
% tailbeats (if depth data is provided).
disp('Dominant Stroking Frequency: ');
disp(dsf);
disp('Number of Tailbeats Measured: ');
disp(length(useableflukebeats(:,3)));

useableflukebeatsDeep = [];
if ~isempty(Depth);
    useableflukebeatsDeep = useableflukebeats(find(useableflukebeats(:,9) == 1),:);
    dsfdeep = mean(useableflukebeatsDeep(:,7));
    disp('Dominant Stroking Frequency (Only Beats Below Surface Threshold): ');
    disp(dsfdeep);
    disp('Number of Tailbeats Measured (Only Beats Below Surface Threshold): ');
    disp(length(useableflukebeatsDeep(:,3)));
end

useableflukebeatsSteady = [];
if

% This segment creates the Tailbeats variable with the start position, end
% position, oscillatory frequency, mean depth (if depth data is
% provided), and mean speed (if speed data is provided) for each tailbeat
% or each non-surface tailbeat (if depth data is provided).
Tailbeats = [];
if ~isempty(Depth);
    Tailbeats = array2table(useableflukebeatsDeep(:,[1,2,7,8]));
    Tailbeats.Properties.VariableNames = {'BeatStartInd','BeatEndInd','BeatOsFreq','BeatAvgDepth'};
    disp('Depth Data Detected - Tailbeats variable will only include tailbeats below the Surface Threshold');
else
    Tailbeats = array2table(useableflukebeats(:,[1,2,7]));
    Tailbeats.Properties.VariableNames = {'BeatStartInd','BeatEndInd','BeatOsFreq'};
    disp('No Depth Data Detected - Tailbeats variable will include all tailbeats found');
end

if ~isempty(Speed);
    if ~isempty(Depth);
        Tailbeats.BeatAvgSpeed = useableflukebeatsDeep(:,10);
    else
        Tailbeats.BeatAvgSpeed = useableflukebeats(:,10);
    end
    disp('Speed Data Detected - Adding Mean Speeds to Tailbeats variable')
end

close all