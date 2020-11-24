function [newaud, shortmovies] = fixmmreadaud(aud, auddur, dispwarnings)
% mmread can only read up to 16 bit data, so this adjusts for that.  Also
% drops the cats convention of having a second, empty, hydrophone.
% auddur is duration of audio data in seconds
% dispwarnings, set to true if you want to display the warnings for fixed data, false if you do not
global movieNum
shortmovies = [];

totalDuration = auddur;
% if aud.rate == 96000;
%     if dispwarnings
%     disp('Recalculating audio rate, 96000 seems more than it can handle');
%     end
% %     aud.rate = round(sum(cellfun(@length, aud.frames(40:end-40)))/4/sum(diff(aud.times(40:end-39)))); audrate = aud.rate
%     if median(cellfun(@length, aud.frames(40:end-40))) == mean(cellfun(@length, aud.frames(40:end-40)))
%         numsamps = median(cellfun(@length, aud.frames(40:end-40)));
%         timedifs = round(diff(aud.times(40:end-39))*100000)/100000;
%         [ns,bins] = hist(timedifs,unique(timedifs));
%         [~,idx] = sort(-ns);
%         timedifs = bins(idx); timedifs = timedifs(1:2);
%         aud.rate = round(numsamps/4/mean(timedifs)); audrate = aud.rate
%     else error('Can''t identify frame lengths in aud');
%     end
% end
if (totalDuration*1.98<size(aud.data,1)/aud.rate && totalDuration*2.02>size(aud.data,1)/aud.rate) ... % for new tags that are recording in 32 bits but the video read only has 16 bit capabilities
        || aud.rate>70000 % account for really high audio rates that can't be handled by the tag processor
    y = aud.data(1:2:end,1)/2^15+aud.data(1:2:end,2);
    aud.data = y; %zeros(size(y,1),1)
    aud.nrChannels = 1;
    aud.bits = 32;
    clear y;
elseif totalDuration*1.9<size(aud.data,1)/aud.rate && totalDuration*2.1>size(aud.data,1)/aud.rate % if it's close but maybe something is wrong with the timestamps of the videos
    if dispwarnings
    disp(['May be a problem with video # ' num2str(movieNum) '.  May be shorter than expected- suggest re-downloading if possible']);
    end
    shortmovies = [shortmovies movieNum];
    y = aud.data(1:2:end,1)/2^15+aud.data(1:2:end,2);
    aud.data = y;%zeros(size(y,1),1)
    aud.nrChannels = 1;
    aud.bits = 32;
    clear y;
end
if size(aud.data,2) == 2; aud.data = aud.data(:,1); end
if abs(totalDuration-size(aud.data,1)/aud.rate)>6 && aud.rate<70000; error(['whole audio was not recorded in video num ' num2str(movieNum) '. Total length should be ' num2str(totalDuration) ' sec but only ' num2str(size(aud.data,1)/aud.rate) ' sec recorded.  RECOMMEND: use VLC to record a new version of the video in mp4 format, maintaining the video codec, but setting the audio codec to MPEG audio, 128 bit rate, and a sample rate that is an integer multiple of your original sample rate.  Other possible causes are listed here, but try VLC first: Sometimes the fix is restarting matlab (I don''t know why). Another possible error is that the audio sampling rate was not 8, 16 or 32-bit.  The third possibility, if the times are close or if one is close to double of the other, is that the video is a bit corrupted and is shorter than it should be.  Something probably went wrong with the download.']);
elseif totalDuration-size(aud.data,1)/aud.rate>1/25; % if there's more than a 25th of a second missing in the audio (more than one video video frame)
    da = diff(aud.times);
    %             numpoints = cellfun(@length,aud.frames);
    if aud.bits == 8; error ('aud.bits issue'); end
    dadiff = [da median(da)]-cellfun(@length, aud.frames)/2/(aud.bits/16)/aud.rate;  % should account for 16 vs 32 bit problem
    normaldiff = median(cellfun(@length, aud.frames)/2/(aud.bits/16)/aud.rate);
    daI = find(dadiff>1.01*percentile(dadiff,.9)); % find the time stamps that have gaps or that did not record all audio
%     if aud.times(end)-aud.times(1)+3>aud.totalDuration; % if you have the whole audio file
%         clipDur = aud.totalDuration;
%     else clipDur = aud.times(end)-aud.times(1)+length(aud.frames(end))/aud.rate;
%     end
%     if aud.rate>70000; aud.rate = round(size(aud.data,1)/(totalDuration - sum(dadiff(daI)))); disp(['Final Audio rate: ' num2str(aud.rate)]); audrate = aud.rate; end
    
    oi = zeros(round(aud.rate*totalDuration),1); % dummy file
    if ~isempty(daI)
        if daI(end) == length(dadiff); daI(end) = []; end
        if dispwarnings
            disp(['Padded ' num2str(sum(dadiff(daI))) ' seconds interspersed among audio file ' num2str(movieNum) ' with silence to match video length.  If this number is bigger than 0.5, recommend redownloading this video as a few frames might have been missed.  (it also could mean the audio sampling rate is too high for the video processor to handle and frames are skipped, currently this happens for anything over 48 kHz)']);
        end
        audI = 1; oiI = 1;
        for iii = 1:length(daI)
%             fixI = round(aud.times(daI(iii))*aud.rate+length(aud.frames{daI(iii)})/2/(aud.bits/16)+1-aud.times(1)*aud.rate);
%             aud.data = [aud.data(1:min(fixI-1,size(aud.data,1)),:); zeros(floor(dadiff(daI(iii))*aud.rate),size(aud.data,2)); aud.data(fixI:end,:)];
            bI = min(round(aud.times(daI(iii))*aud.rate+length(aud.frames{daI(iii)})/2/(aud.bits/16)-aud.times(1)*aud.rate),length(aud.data)-audI+oiI);
            oi(oiI:bI,1) = aud.data(audI:(bI-oiI)+audI);
            audI = (bI-oiI)+audI + 1;
            oiI = floor(dadiff(daI(iii))*aud.rate) + bI + 1;
        end
    else
        if dispwarnings; disp(['Padded ' num2str(totalDuration-size(aud.data,1)/aud.rate) ' seconds at the end of audio file ' num2str(movieNum) ' with silence to match video length.']); end
        %             aud.data = [aud.data; zeros(floor((totalDuration-size(aud.data,1)/aud.rate)*aud.rate),size(aud.data,2));];
        oi(1:length(aud.data),1) = aud.data;
    end
    aud.data = oi; clear oi;
    %             padding = totalDuration-size(aud.data,1)/aud.rate;
    %             aud.data = [zeros(round(padding*aud.rate),1); aud.data];
    %             disp(['Padded first ' num2str(padding) ' sec of video ' num2str(movN(n)) ' with silence to match video length']);
elseif totalDuration+1/25<size(aud.data,1)/aud.rate && dispwarnings; disp([num2str(abs(totalDuration-size(aud.data,1)/aud.rate)) ' s more audio exists than video in this frame']);
end
newaud = aud;