function [RMS, timeoffset] = CATSrms(aud,fs,df)
%input an audio structure from mmread of CATS video.  df was recommended in
%d3rms.  Because our audio rate is much lower than theirs (22000 vs. 48-192
%kHz, and they use 100, 10 seems reasonable.  output is "decibels" re the
%number of uPa that would give a reading of 1 in the microphone.
% timeoffset is in seconds

if nargin<3; df = 10; end
if nargin<2; fs = 5; end
afs = aud.rate;
nov = round((1-1/fs)*afs/df); %allow for 0.8*fs samples of overlap (for 5 hz)...so the one second analysis window slides forward by 1/fs sec per time ->  fs hz data
intbins = true;
if abs(nov - (1-1/fs)*afs/df) > .01 % if your bin size isn't an integer, you have a problem
    try F = factor(afs/fs);
        [d,I] = min(abs(F-df));
        if d/df <.5 && afs/df > 1100 %if the adjusted decimation factor is between 50% and 150% of the old one and you still get a sampling rate above 1100 Hz, roll with it
            df = F(I);
            nov = round((1-1/fs)*afs/df);
            disp('adjusted decimation factor of audio file to result in integer analysis bins');
        else
            error('decimation factor chosen does not result in integer 1 second bins');
        end
    catch
        intbins = false;
    end
end
if sum(aud.data(:,1)==0) == length(aud.data) ||...
        sum(diff(aud.data(:,1))==0)>.95*length(aud.data(:,1)) ||...
        sum(isnan(aud.data(:,1)))>.95*length(aud.data(:,1))
    disp('audiodata column 1 is empty, using audio channel 2');
    x = aud.data(:,2);
elseif size(aud.data,2) == 1 || sum(diff(aud.data(:,2))==0)>.95*length(aud.data(:,2)) || sum(isnan(aud.data(:,2)))>.95*length(aud.data(:,2))
x = aud.data(:,1); % just use one channel of the data
else
    maxRMS1 = nanmean(aud.data(:,1)); maxRMS2 = nanmean(aud.data(:,2)); 
    figure(200); clf; plot(aud.data); legend('Channel 1', 'Channel 2');
    if maxRMS1<maxRMS2; x = aud.data(:,1); disp('First column is lower gain, using that audio channel');
    else x = aud.data(:,2); disp('Second column is lower gain, using that audio channel');
    end
end
timeoffset = aud.totalDuration-length(x)/afs; %the amount of time after the start of the video that the audio starts
xd = decdc(x,df); %x decimated
RMS = ones(floor((aud.totalDuration-1)*fs),1);
RMS2 = ones(floor((aud.totalDuration-1)*fs),1);
%below from d3rms by Alison and Stacy
xdf = fir_nodelay(xd,128,[66/(afs/df/2) , 94/(afs/df/2)]); % x decimated and bandpass filtered between 66 and 94 Hz
n = afs/df; %analyse 1 sec of acoustic data at a time -- there are afs/df samples per second in the decimated data
if intbins
    [X,~,~]= buffer(xdf,n,nov,'nodelay') ; % n is the number of samples per chunk; nov is the amount of overlap in samples, goes until the last complete second can be put in a bin (will end on a mulitple of .2 seconds)
else
%     nov = (1-1/fs)*afs/df;
    X = nan(n,length(RMS));
    k = 1;
    nk = 1;
    kplus = round(nk/fs*afs/df);
    ktot = 0;
    for i = 1:length(X(1,:))
        X(:,i) = xdf(k:k+n-1);
        k = k+kplus;
%         ktot = ktot+kplus;
%         if ktot >= n
%             ktot = 0;
%             nk = 0;
%         end
        nk = nk+1;
        kplus = round(nk/fs*afs/df-round((nk-1)/fs*afs/df));
    end
end

% FYI, the second output to buffer are any values leftover that weren't at
% least .2 seconds long.  The third output is the last .8 seconds that were
% included in the previous bin.

timeoffset = timeoffset + .5; % add .5 so that the RMS flow noise bins are centered (one second long each)
for j = 1:size(X,2) %process one column of the buffered data at a time
    y = X(:,j); %the data to use this iteration of the loop (1 second worth of acoustic data)
    RMS(j) = 20*log10(std(y));% + SensTag3; % calculate rms level in dB re 1 uPa
    RMS2(j) = 20*log10(sqrt(mean(y.^2))); %standard definition.  d3rms uses std which implies they subtract the mean of each second from each value.  Unsure why.
%     if isnan(RMS2(j))
%      disp(3)
%     end
    
    clear y %just in case, clean up!
end



