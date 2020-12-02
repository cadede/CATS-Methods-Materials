global movieNum
shortmovies = [];
if ~isempty(wavfiles)
    aud = struct();
    [aud.data,aud.rate,aud.bits] = wavread(wavstr{1},[1 10]);
else
    [~,aud] = mmread([movieloc m2{find(~cellfun(@(x) strcmp(x(end-2:end),'raw'),m2),1)}],[],[3 4],true);
end
audrate = aud.rate;
if aud.rate == 96000; warning('audio rates higher than 48 kHz may have problems, suggesting downsampling the audio using a different software first'); end
%     disp('Recalculating audio rate, 96000 seems more than it can handle');
%     if median(cellfun(@length, aud.frames(40:end-40))) == mean(cellfun(@length, aud.frames(40:end-40)))
%         numsamps = median(cellfun(@length, aud.frames(40:end-40)));
%         timedifs = round(diff(aud.times(40:end-39))*100000)/100000;
%         [ns,bins] = hist(timedifs,unique(timedifs));
%         [~,idx] = sort(-ns);
%         timedifs = bins(idx); timedifs = timedifs(1:2);
%         aud.rate = round(numsamps/4/mean(timedifs)); audrate = aud.rate
%     else error('Can''t identify frame lengths in aud');
%     end
% %     aud.rate = round(sum(cellfun(@length, aud.frames(40:end-40)))/4/sum(diff(aud.times(40:end-39)))); audrate = aud.rate
% end

% LOOK AT videos from 4k to see if zeros are being added unncessarily to
% those files.

 audioloc = movieloc; 
 checkaudiofs = false; % may have to enable this flag if audio sample rates are funny
for n = 1:length(movies) %for some reason if this is in with the next for loop it messes up
    clear aud; clear totalDuration;
    wavfile = []; sm = [];
    if ~strcmp(movies{n}(end-2:end),'raw')
        movieNum = movN(n);
        if ~checkaudiofs; audioend = movies{n}(end-2:end); end
        if ~isempty(wavstr); wavfile = wavstr{cellfun(@(x) strcmp(movies{n}(1:end-3),x(1:end-3)),wavfiles)}; end
        if ~isempty(wavfile)
            aud = struct();
            [aud.data,aud.rate,aud.bits] = wavread(wavfile);
            vid= mmread([movieloc movies{n}], [1 3]); %just audio
            aud.nrChannels = size(aud.data,2);
            aud.totalDuration = size(aud.data,1)/aud.rate;
            if aud.nrChannels == 2 && sum(aud.data(:,2)==0) == length(aud.data(:,2)); aud.data(:,2) = []; aud.nrChannels = aud.nrChannels - 1; end
            totalDuration = aud.totalDuration;
            if size(aud.data(:,1),1) == 0 || vid.totalDuration - aud.totalDuration > .5;
                warning(['Video ' num2str(movieNum) ' duration is ' num2str(vid.totalDuration) ', while its wav file is ' num2str(aud.totalDuration); 's']);
                sm = movN(n);
            end
        else
            [~,aud]= mmread([audioloc movies{n}(1:end-3) audioend], [],[],true); %just audio
            totalDuration = aud.totalDuration;
            [aud, sm] = fixmmreadaud(aud,totalDuration,true);
            if checkaudiofs;
                if abs(aud.rate-audrate)>.1 % if the audio rate of the new movies is not the same as the old movies
                    if abs(round(aud.rate/audrate)-aud.rate/audrate) > .1 || audrate>aud.rate;
                        error(['Files from which audio should be read must be sampled at a rate that is an integer multiple of the original movie rate (' num2str(audrate) ' Hz) but were instead sampled at ' num2str(aud.rate) ' Hz.']);
                    else
                        aud.data = decdc(aud.data,aud.rate/audrate);
                        disp(['Warning: resampled audio from ' num2str(aud.rate) ' Hz to ' num2str(audrate) ' Hz to match original audio file']);
                        aud.rate = audrate;
                    end
                else audrate = aud.rate;
                end
            else
                audrate = aud.rate;
            end
        end
        shortmovies = [shortmovies sm];
    else
        fid = fopen([movieloc movies{n}]); 
        y = fread(fid,'int16'); bits = 16; fclose(fid);
        if sum(y(2:2:end)==0) ~= length(y)/2
            fid = fopen([movieloc movies{n}]); 
            y = fread(fid,'int32'); bits = 32; fclose(fid);
            if sum(y(2:2:end)==0) ~= length(y)/2
               error('Audio bit rate unknown (not 16 or 32 bits)'); 
            end
        end
        
        y = y/2^(bits-1);
        aud.data = y(1:2:end); % y(2:2:end)
        aud.rate = audrate;
        aud.bits = bits;
        aud.nrChannels = 1;
        aud.totalDuration = size(aud.data,1)/aud.rate;
        totalDuration = aud.totalDuration;
        clear y;
    end
    
    if isempty(wavfile)
        try wavwrite(aud.data,aud.rate,aud.bits,[DIR movies{n}(1:end-4) '.wav'])
        catch %assumes the error was that the file was too large
            for i = 1:aud.totalDuration/60/60+1 % save one hour increment files
                fname = [movies{n}(1:end-4) 'hour' num2str(i) '.wav'];
                newDN = datestr(datenum(fname(min(regexp(fname,'-'))+1:max(regexp(fname,'-'))-1),'yyyymmdd-HHMMSS')+(i-1)*1/24,'yyyymmdd-HHMMSS');
                fname(min(regexp(fname,'-'))+1:max(regexp(fname,'-'))-1) = newDN;
                wavwrite(aud.data((i-1)*aud.rate*60*60+1:i*aud.rate*60*60,:),aud.rate,aud.bits,[DIR fname])
            end
        end
    else try movefile(wavfile,[DIR movies{n}(1:end-4) '.wav']); catch; warning(['could not move file ' movies{n}(1:end-4) '.wav into audioData directory']); end
    end
%     if audrate>70000
%         aud.data = decdc(aud.data,2); % decimate to make them reasonable
%         aud.rate = aud.rate/2;
%     end
    lastwarn('');
    try
        save([DIR movies{n}(1:end-4) 'audio.mat'],'aud','totalDuration');
        if ~isempty(lastwarn)
            error(lastwarn);
        end
    catch %v7.3 allows for bigger files, but makes a freaking huge file if used when you don't need it
        save([DIR movies{n}(1:end-4) 'audio.mat'],'aud','totalDuration','-v7.3');
        disp(['Made a version 7.3 file (audio ' movies{n}(1:end-4) ' was big)']);
    end
end
disp('Check wav files for accuracy/length, then they can be deleted (keep the .mat files for audio analysis).');%  If the wav files are the wrong length, reimport the audio files.  The frametimes etc. should still be accurate and should not need to be redone.');
if ~isempty(shortmovies); disp(['Audio lengths are off in videos: ' num2str(shortmovies) '.  This may suggest a problem with the download, recommend redownloading if possible. This message will repeat.']); end