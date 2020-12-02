function [DB,AUD] = getflownoise(audiodir,vars)

names = fieldnames(vars);
for i = 1:length(names)
    eval([names{i} ' = vars.' names{i} ';']);
end


DB = nan(size(Depth,1),1);
DNaud = DN(1):1/fs/24/60/60:DN(end);
DBaud = nan(size(DNaud));
audiofiles = dir([audiodir '*audio.mat']); % first look for audio files that have already been converted to mat files
if isempty(audiofiles) % check for other wav files not from movies (acousonde/dtag etc.)
    wavfiles = dir([audiodir '*.wav']);
    if ~isempty(wavfiles);
        for n = 1:length(wavfiles) %for some reason if this is in with the next for loop it messes up
            clear aud;
%             [~,aud]= mmread([audiodir '' wavfiles(n).name], [],[],true); %just audio
aud = struct();
            [aud.data,aud.rate,aud.bits] = wavread([audiodir '' wavfiles(n).name]);
            aud.nrChannels = size(aud.data,2);
            aud.totalDuration = size(aud.data,1)/aud.rate;
            lastwarn('');
            try
                save([audiodir '' wavfiles(n).name(1:end-4) 'audio.mat'],'aud');
                if ~isempty(lastwarn)
                    error(lastwarn);
                end
            catch %v7.3 allows for bigger files, but makes a freaking huge file if used when you don't need it
                save([audiodir '' wavfiles(n).name(1:end-4) 'audio.mat'],'aud','-v7.3');
                disp('Made a version 7.3 file (audio was big)');
            end
        end
        noaud = false;
        audiofiles = dir([audiodir '*audio.mat']); % first loo
    else noaud = true;
    end
else noaud = false;
end

if ~nocam || ~noaud
    for i = 1:length(audiofiles)
        load([audiodir '' audiofiles(i).name]);
        AUD.rate = aud.rate;
        AUD.bits = aud.bits;
        AUD.nrChannels = aud.nrChannels;
        if aud.rate == 22050
            DBdf = 7; %7 is the decimation factor that allows for integer subsampling bins
        elseif aud.rate == 48000 || aud.rate == 12000 || aud.rate == 24000
            DBdf = 15;
        elseif aud.rate == 25811
            DBdf = 53;
        else error('new sampling rate, check your decimation factor to ensure an integer bin');
        end
        try
        [DBt, offset] = CATSrms(aud,fs,DBdf); %offset is in seconds
        catch
            disp(['audio ' num2str(audiofiles(i).name) ' could not be read']);
            continue;
        end
%         if tagnum == 50 || tagnum == 51
%              vidnum = str2num(audiofiles(i).name(regexp(audiofiles(i).name,'_')+1:max(regexp(audiofiles(i).name,'\d'))));
%              if isempty(vidnum); vidnum = str2num(audiofiles(i).name(regexp(audiofiles(i).name,'C')+1:max(regexp(audiofiles(i).name,'\d')))); end
%         elseif kitten
            vidnum = audiofiles(i).name(regexp(audiofiles(i).name,'\d'));
            vidnum = str2num(vidnum(end-3:end));
            disp(['Video number ' num2str(vidnum) ' being read, sample rate is ' num2str(aud.rate) ' Hz']);

%         elseif tagnum<20 && tagnum>12
%             vidnum = i;
%         else
%             vidnum = str2num(audiofiles(i).name(regexp(audiofiles(i).name,'\d')));
%         end
        try audStart = vidDN(vidnum)+offset/24/60/60;
        catch
            if i == 1; disp('assuming audiofiles start at beginning of file and are continuous');
                audStart = data.Date(1)+data.Time(1)+offset/24/60/60;
            elseif (vidnum>length(vidDN) || isempty(vidDN)) && ~exist('lastDur','var')
                warning(['audio ' num2str(vidnum) ' does not seem to be on whale, skipping']);
                continue
            else
                audStart = audStart+lastDur;
            end
            lastDur = aud.totalDuration/24/60/60;
        end
        [~,I] = min(abs(DNaud-audStart));
        if I~=1 % if the audio starts before the video starts
            DBaud(I:I+length(DBt)-1) = DBt;
        else
            [~,I] = min(abs(DNaud-(audStart+length(DBt)/fs/24/60/60))); % from the end
            DBaud(1:I) = DBt(length(DBt)-I+1:end);
        end
        for j = 1:length(DBt)
            tI = find(DN>=audStart+(j-1)*1/fs/24/60/60-(1/fs/2)/24/60/60 & DN<audStart+(j-1)*1/fs/24/60/60+1/fs/2/24/60/60); % find the times within 1/fs/2 seconds of the calculated DB value (basically acounts for the difference in sample size, and non-linearly related sample sizes)
            DB(tI) = DBt(j);
        end
        disp(['Audio file ' num2str(i) ' completed']);
        % should be faster but doesn't seem to be, there is also a difference
        % between above for some reason
        %     Xtime = DNaud(I:I+length(DBt)-1);
        %     %     JiggleRMS = nan(size(DN));
        %     k = 1; [~,j] = min(abs(DN-(Xtime(1)-1/fs/2/24/60/60))); DB(j) = DBt(k);
        %     for k = 1:length(Xtime);
        %         j2 = find(DN(j:min(j+fs,length(DN)))<=Xtime(k)+1/fs/2/24/60/60,1,'last')+j-1;% find the times that are within 1/fs/2 seconds of Xtime(k)
        %         if isempty(j2); [~,j] = min(abs(DN-(Xtime(k)-1/fs/2/24/60/60)));[~,j2] = min(abs(DN-(Xtime(k)+1/fs/2/24/60/60))); end
        %         DB(j:j2) = DBt(k);
        %         j = j2+1;
        %     end
        
    end

end
DB(isinf(DB)) = nan;