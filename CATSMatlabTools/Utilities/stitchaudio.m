function stitchaudio(audiodir,whaleName,starttime,vidDN,prhloc)

% input the location of the audiofiles and start times of each file
% stitches together all audio files into a continuous wav file that matches
% the start time of the prh file.  inputs 2:5 are optional, if not present
% will ask you to import a prh file


if nargin<2 || isempty(vidDN); d = pwd; cd(audiodir); [prhfile,prhloc] = uigetfile('*.mat','select prh file with vidDN, DN and INFO vars');
    cd (d);
    load([prhloc prhfile],'vidDN','INFO','DN');
    starttime = DN(1);
    whaleName = INFO.whaleName;
end
if ~strcmp(audiodir(end),'\'); audiodir = [audiodir '\']; end
if ~strcmp(prhloc(end),'\'); prhloc = [prhloc '\']; end
% names = fieldnames(vars);
% for i = 1:length(names)
%     eval([names{i} ' = vars.' names{i} ';']);
% end
audiofiles = dir([audiodir '*audio.mat']); % first look for audio files that have already been converted to mat files
if isempty(audiofiles) % check for other wav files not from movies (acousonde/dtag etc.)
    wavfiles = dir([audiodir '*.wav']);
    if ~isempty(wavfiles)
        for n = 1:length(wavfiles) %for some reason if this is in with the next for loop it messes up
            clear aud;
            %             [~,aud]= mmread([audiodir '' wavfiles(n).name], [],[],true); %just audio
            aud = struct();
            try [aud.data,aud.rate,aud.bits] = wavread([audiodir '' wavfiles(n).name]);
                    catch
            audioI = audioinfo([audiodir '' wavfiles(n).name]);
            [aud.data,aud.rate] = audioread([audiodir '' wavfiles(n).name]);
            aud.bits = audioI.BitsPerSample;
        end
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
    else; noaud = true;
    end
else; noaud = false;
end
AUD = struct();
audStart = 0;
if isempty(vidDN); disp('assuming audiofiles are continuous with no gaps'); lastDur = 0; end

audiodata = zeros(0,1);
firstaudio = true;
if ~noaud
    for i = 1:length(audiofiles)
        vidnum = audiofiles(i).name(regexp(audiofiles(i).name,'\d'));
        vidnum = str2num(vidnum(end-3:end));
        %         if audStart>DN(find(tagondec,1,'last')) || ((vidnum>length(vidDN) || isempty(vidDN)) && ~exist('lastDur','var'))
        %                 warning(['audio ' num2str(vidnum) ' does not seem to be on whale, skipping']);
        %                 continue
        %         end
        load([audiodir '' audiofiles(i).name],'aud');
        if i == 1
            AUD.rate = aud.rate;
            AUD.bits = aud.bits;
            AUD.nrChannels = size(aud.data,2);%aud.nrChannels;
        else if AUD.rate~=aud.rate; error('Sample Rates of Audio do not match'); end
            if AUD.nrChannels>size(aud.data,2); warning(['Reducing wav file to ' size(aud.data,2) ' channel(s)']);
                AUD.nrChannels = size(aud.data,2);
                audiodata = audiodata(:,1:size(aud.data,2));
            end
        end
        disp(['Audio number ' num2str(vidnum) ' being read, sample rate is ' num2str(aud.rate) ' Hz']);
        
        %         elseif tagnum<20 && tagnum>12
        %             vidnum = i;
        %         else
        %             vidnum = str2num(audiofiles(i).name(regexp(audiofiles(i).name,'\d')));
        %         end
        afile = audiofiles(i).name;
        try audStart = vidDN(vidnum);%+offset/24/60/60;
        catch
            try audStart = datenum(afile(min(regexp(afile,'-'))+1:max(regexp(afile,'-'))-1),'yyyymmdd-HHMMSS-fff');
            catch; audStart = datenum(afile(min(regexp(afile,'-'))+1:max(regexp(afile,'-'))-1),'yyyymmdd-HHMMSS'); disp('could only read timestamp to nearest second');
            end
        end
        
        if i == 1 || firstaudio; STARTTIME = audStart;
            if STARTTIME>starttime; audiodata = zeros(round(AUD.rate*(STARTTIME-starttime)*24*60*60),AUD.nrChannels); STARTTIME = starttime; end
            if STARTTIME<starttime; II = round((starttime-STARTTIME)*24*60*60*AUD.rate)+1; if II>length(aud.data); continue; end
                aud.data = aud.data(II:end,:); STARTTIME = starttime; audStart = STARTTIME;
            end
        end
        %             if i == 1; disp('assuming audiofiles start at beginning of file and are continuous');
        %                 audStart = DN(1) +offset/24/60/60;
        %             elseif (vidnum>length(vidDN) || isempty(vidDN)) && ~exist('lastDur','var')
        %                 warning(['audio ' num2str(vidnum) ' does not seem to be on whale, skipping']);
        %                 continue
        %             else
        %                 audStart = audStart+lastDur;
        %             end
        %             lastDur = aud.totalDuration/24/60/60;
        %     end
        
        I = round((audStart-STARTTIME)*24*60*60*AUD.rate)+1;
%         [~,I] = min(abs(DNaud-audStart));
%         if I~=1 % if the audio starts before the video starts
            audiodata(I:I+length(aud.data)-1,1:AUD.nrChannels) = aud.data(:,1:AUD.nrChannels);
%         else
%             [~,I] = min(abs(DNaud-(audStart+length(DBt)/fs/24/60/60))); % from the end
%             audiodata(1:I) = DBt(length(DBt)-I+1:end);
%         end
%         for j = 1:length(DBt)
%             tI = find(DN>=audStart+(j-1)*1/fs/24/60/60-(1/fs/2)/24/60/60 & DN<audStart+(j-1)*1/fs/24/60/60+1/fs/2/24/60/60); % find the times within 1/fs/2 seconds of the calculated DB value (basically acounts for the difference in sample size, and non-linearly related sample sizes)
%             DB(tI) = DBt(j);
%         end
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
        firstaudio = false;
    end
   if length(audiodata)>20*60*60*24000 % split in to parts < 4 GB each
       k = 1; 
       for i = 1:20*60*60*24000:length(audiodata)
           audiowrite([prhloc,whaleName,'_fullaudiopart',num2str(k),'-',datestr(starttime+(i-1)/AUD.rate/24/60/60,'yyyymmdd-HHMMSS-fff'),'.wav'],audiodata(i:min(i+20*60*60*24000-1,length(audiodata)),:),AUD.rate);
           k = k+1;
       end
   else
       audiowrite([prhloc,whaleName,'_fullaudio-',datestr(starttime,'yyyymmdd-HHMMSS-fff'),'.wav'],audiodata,AUD.rate);
   end

end
disp('full audio wav file created (start time is the same as prh file starttime)');