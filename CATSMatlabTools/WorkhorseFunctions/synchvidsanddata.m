function [camon,audon,vidDN,vidDurs,nocam,tagslip,vidadj,audstart] = synchvidsanddata(data,headers,tagon,viddata,Hzs,DN,ODN,fs,CAL,nocam,synchusingvidtimestamps,useFrames)

% this function looks to see if any adjustment is needed for the video and
% data times, based on inputs from the meta data xls file, and then also
% provides the option to synch videos and data using surfacing behavior
% data or some other marker of time from the videos that can be matched
% with data.

dbstop if error;
global fileloc filename
if nargin<12; useFrames = false; end %this is a legacy switch for if you enter framenumbers into the excel sheet instead of times
if sum(diff(data.Pressure)<.001) == length(data.Pressure); nopress = true; else nopress = false; end
 Atemp = ([data.Acc1 data.Acc2 data.Acc3]-repmat(CAL.aconst,size(data,1),1))*CAL.acal; %temp Acceleration file for guessing at tagslip location
% synch 
try if viddata.vid4k
        synchaudio = input('Synch audio files with data (for CATS tags, select yes if audio was downloaded as a single audio file and then split in an earlier step)? (1 = yes, 2 = no) ');
        
    else audstart = [];
    end
catch
    audstart = []; viddata.vid4k = false;
    disp('To synch an audio file, add a true ''vid4k'' boolean variable to your movieTimes file');
end
cf = pwd;
if viddata.vid4k && ~isempty(synchaudio) && synchaudio == 1
    timestamp = input('Use tag on animal time as synch point? (1 = yes, 2 = no) ');
    if timestamp == 2; timestamp = input('Enter time of synch point in date vector format ([yyyy mm dd HH MM SS]): ');
        timestamp = datenum(timestamp);
    else timestamp = find(tagon,1); timestamp = DN(timestamp);
    end
    try cd([fileloc 'AudioData\']); catch; cd(fileloc); end
    [~,I] = min(abs(DN-timestamp)); I = I-fs*10:I+fs*10;
    [audfile,audloc] = uigetfile('*.wav','Get first audio file from "AudioData" folder (assumes all audiofiles from here are consecutive with no gaps)');
    if exist([audloc audfile(1:end-4) 'audio.mat'],'file'); load([audloc audfile(1:end-4) 'audio.mat']); else
        aud = struct();
        try [aud.data,aud.rate,aud.bits] = wavread([audloc audfile]);
        catch
            audioI = audioinfo([audloc audfile]);
            [aud.data,aud.rate] = audioread([audloc audfile]);
            aud.bits = audioI.BitsPerSample;
        end
    end
    FF = figure(100); clf;
    audioI = round((timestamp-DN(1))*24*60*60*aud.rate); audioI = audioI-aud.rate*10:audioI+aud.rate*10;
    sp1 = subplot(2,1,1); hold on;
    ax = plotyy(DN(I),data.Pressure(I),DN(I),Atemp(I,3)); legend('Depth','Az')
    set(ax(1),'ydir','rev')
    set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
    set(ax,'xlim',[DN(I(1)) DN(I(end))]);
    t = title('Click on synch location in upper graph (will select peak within 1s)');
    set(ax,'nextplot','add')
    subplot(2,1,2);
    pI = (1:length(audioI))/aud.rate;
    plot(pI,aud.data(audioI));
    xlabel('Seconds')
    xlim([pI(1) pI(end)])
    [x,~] = ginput(1);
    [~,x] = min(abs(DN-x));
    [i1,mag] = peakfinder(Atemp(x-fs:x+fs,3));
    [~,i] = max(mag); i1 = i1(i);
    x = x-fs-1+i1;
    plot(ax(2),DN(x),mag(i),'rs','markerfacecolor','r')
    set(ax,'xlim',[DN(x-10*fs) DN(x+10*fs)])
    set(t,'string','');
    subplot(2,1,2); hold on;
    title('Click on synch location in lower graph (will select peak within 1s)');
    [x2,~] = ginput(1);
    x2 = x2*aud.rate;
    x2 = audioI(1)-1+round(x2);
    [i1,mag] = peakfinder(aud.data(x2-aud.rate:x2+aud.rate));
    [~,i] = max(mag); i1 = i1(i);
    x2 = x2-aud.rate-1+i1;
    plot(pI(x2-audioI(1)+1),mag(i),'rs','markerfacecolor','r')
    audioI2 = x2-audioI(1)+1-10*aud.rate:x2-audioI(1)+1+10*aud.rate;
    xlim([audioI2(1) audioI2(end)]/aud.rate)
    audstart = DN(x)-x2/aud.rate/24/60/60;
    disp(['Old Audio start time: ' datestr(DN(1), 'dd-mmm-yyyy HH:MM:SS.fff')]);
    disp(['New Audio Start time: ' datestr(audstart,'dd-mmm-yyyy HH:MM:SS.fff') '. Will be asked to rename audio files after prh process is complete.']);
%     close(FF);
else
    audstart = nan;
end

cd(cf); 

try if sum(data.Camera) == 0; nocam = true; else nocam = false; end; catch; disp('Could not automatically detect camera status, maintaining user input'); end
if strcmp(headers{5,1},'Time Style'); rowstart = 8; else rowstart = 7; end 
times = cell2mat(headers(rowstart:end,2:3));
%if no movie files
if nocam
    tagslip = [1 1]; %tagslipC = 1; % confidence of tagslip.  1 if you see it move on a video, 0 if you estimate based on max jerk
    camon = false(size(DN)); audon = false(size(DN));
    if exist('pconst','var'); p = (data.Pressure-pconst)*pcal;
        %     [~,peakmag] = peakfinder(p,3,min(p)+6,-1); p = p-max(min(peakmag),max(peakmag)-2);
    else p = data.Pressure;
    end
 vidDN = nan; vidDurs = 0;
else
     nocam = false;
    
    names =fieldnames(viddata);
    for ii = 1:length(names)
        eval([names{ii} ' = viddata.' names{ii} ';']);
    end
    if any(vidDN<1); vidDN = vidDN+floor(DN(1)); 
        warning('No date number on vidDN, adding the following:');
        for ii = 2:length(vidDN);  if vidDN(ii)<vidDN(ii-1); vidDN(ii:end) = vidDN(ii:end)+1; end; try disp(['vid number ' num2str(ii) ' date: ' datestr(vidDN(ii),'yyyy-mmm-dd')]); catch; end; end
    end
     
    if synchusingvidtimestamps
        UTC = Hzs.UTC;
        vidtimedif = cell2mat(headers(3,2));
        vidDNorig = vidDN;
        if vidtimedif ~=0; warning('Time Difference does not = 0. Do video starttimes need to be adjusted by timedif?');
            pp = input('1 = yes, 2 = no  ');
            if pp == 2; vidtimedif = 0;
            end
        end
        if data.Date(1)<datenum([2017 09 01]) % old versions of tag where time stamps were in UTC
            vidDN = vidDN + data.Date(1)+UTC/24+vidtimedif/24;
            disp('Old file detected, assuming time stamps are in UTC');
        elseif vidDN(find(~isnan(vidDN),1))<10 % confirm that new vidDN takes into account full date
            vidDN = vidDN + data.Date(1)+vidtimedif/24;
        else vidDN = vidDN + vidtimedif/24;
        end
        if sum(vidDN<DN(1)) + sum(vidDN>DN(end))>1; warning('multiple vidDN are outside range of data. This may be okay if more videos are included in movieTimes outside of those on whale (so just continue), but check that any timedif in xls header file should apply to both data and video.  If, for example, data were downloaded on a computer in a different timezone than for which the tag was programmed, video and data will be offset by different amounts and will have to be adjused for, potentially by unhighlighting the line above this one'); end
    end
    
      if isempty (frameTimes); error('Put frametimes from movies in the same folder as the movies'); end
        whaleName = char(headers(1,2));
        timedif = cell2mat(headers(3,2)); % The number of hours the tag time is behind (if it's in a time zone ahead, use -).  Also account for day differences here (as 24*# of days ahead)
        try
            videonum = cell2mat(headers(rowstart:end,1)); todel = find(~isnan(videonum),1,'last')+1:length(videonum);
            types = headers(rowstart:end,4);
            times = cell2mat(headers(rowstart:end,2:3));
            times(times==0) = nan;
            surfI = find(strcmp(types,'Surface'));
            videonum(todel) = []; times(todel,:) = []; types(todel) = [];
            if any(any(times>1)); useI = true; else useI = false; end% if you don't enter the indices, but use the mm:ss.ff format
            timesM = times(strcmp(types,'Manual'),2);
            % old line, times should be in days since video start
%             if kitten;  times = (times - repmat((vidDN([videonum videonum])- floor(vidDN([videonum videonum])))',size(times,1),1))*24*60*60; end
        catch
            tagslip = [1 1];% tagslipC = 1; % confi
            times = [nan nan];
            timesM = [];
        end
        
        % this is legacy code, useful if instead of putting timestamps in
        % the excel sheet, you'd rather put in frame numbers
        if any(any(~isnan(times))) && ~synchusingvidtimestamps && useFrames%&& ~kitten
            if ~useI
                addframes = round((times*24*60*60-floor(times*24*60*60))*100); %number of frames past the seconds value
                addframes3 = round((times*24*60*60-floor(times*24*60*60))*1000);
                if ~any(any(abs(addframes3-addframes*10)>=1)) % tries to determine if you entered 3 digit times or 2 digit frame numbers.
                    times = floor(times*24*60*60); % convert to seconds without the "frames" decimal
                    addF = true;
                else
                    times = times*24*60*60;
                    addF = false;
                end
            else addF = true;
            end
            if addF
                for i = 1:length(videonum)
                    for j = 1:2
                        if ~useI
                            [~,tI] = find(frameTimes{videonum(i)}>times(i,j),1,'first'); %frame 0 is the first frame that starts with those seconds.  In the video, there may be some (maybe 6?) missing frames in the first "0" second that are recorded in frameTimes.
                            tI = tI+addframes(i,j);
                        else tI = times(i,j);
                        end
                        if ~isnan(tI); times(i,j) = frameTimes{videonum(i)}(tI); end
                    end
                end
            end
        end
        if nopress; p = data.Pressure; else
            if exist('pconst','var'); p = (data.Pressure-pconst)*pcal;
                %             [~,peakmag] = peakfinder(p,3,min(p)+6,-1); p = p-max(min(peakmag),max(peakmag)-2);
            else p = data.Pressure;
            end

        end
        % took out DN because it is automatically imported from the
        % function
%         DN = data.Date+data.Time+timedif/24;

        %adjust the video date numbers by assuming video 1 started with the start
        %of the tag
        ODNa = ODN + timedif/24;
        if ~synchusingvidtimestamps
            if (vidDN(1)<DN(1)-1 || vidDN(1)>DN(end)+1) % if the vidDN don't relate at all to the data times
                if ~isnan(vidDN(1)) % legacy- if there is no info on the video start times, assume it starts just after data turns on as a first case
                    warning('Video 1 has a start time apparently out of the data, assuming it starts 5 seconds after data turns on');
                    vidDN = vidDN - (vidDN(1)-(ODNa-floor(ODNa))+5/24/60/60); %assumes the 1st video starts 5 seconds after you turn on the data
                end
                vidDN = vidDN +floor(DN(1)); % account for the date
            end
            vidDNorig = vidDN;
        end
        timesM = timesM + floor(DN(1));
        times(strcmp(types,'Manual'),2) = timesM;
        
        %         if isempty(strfind(movies{1},'CATS')) && isempty(strfind(movies{1},'AW')) %kitten is the new wireless cam without timestamps on the videos
        %             kitten = true; else kitten = false;
        %         end
        %         if kitten
        %             II = find(strcmp(types,'Start'));
        %             for i = 1:length(II)
        %                 vidDN(videonum(II(i))) = times(II(i),2)+floor(DN(1)) - times(II(i),1);
        %             end
        %         end
        
        tagslip = [1 1];% tagslipC = 1; % confidence of tagslip.  1 if you see it move on a video, 0 if you estimate based on max jerk
        njerkTemp = (9.81*fs)*sqrt(diff(Atemp).^2*ones(3,1)); %temp jerk for guessing at tag slip
        
        
%         timeDN = times/24/60/60; %times was in seconds
        timeDN = times; % now time stamps are always time (not seconds), and timeDN is time since video start
        timesO = times; 
        lastgood = 1; FIX = [];
        for i = 1:length(frameTimes)
            if ~synchusingvidtimestamps
                surfI = find(videonum == i & strcmp(types,'Surface'));
                manI = find(videonum == i & strcmp(types,'Manual')); % if the whale does not surface, can do a manual calibration using another factor (tagoff, no paddlewheel to paddlewheel audio etc.)
                noneI = find(videonum == i & strcmp(types,'None'));
                if isempty(frameTimes{i}); continue; end
                if strcmp(headers{5,1},'Time Style') && strcmp(headers{5,2},'Embedded'); timeDN(surfI,:) = timeDN(surfI,:) - (vidDN(i)-floor(vidDN(i))); end
                times = timeDN*24*60*60;
                if (~isempty(surfI) || ~isempty(manI))&&(~synchusingvidtimestamps)
                    if ~isempty(surfI)
%                         if i == min(videonum) || i == videonum(find(cellfun(@(x) strcmp(x,'Surface'),types),1));
%                             dives = finddives2(p,fs,min(p)+1,min(p)+2,true);
%                             
%                             I = round(find(p>1,1,'first')+5*fs); %time for the camera to turn on after the camera is deployed
%                         else [~,I] = min(abs(DN-(vidDN(i)+timeDN(surfI(1),1)))); end
                        [~,I] = min(abs(DN-(vidDN(i)+timeDN(surfI(1)))));
                        peakLoc = nan(size(surfI));
                        if any(isnan(p)); p = fixgaps(p); p(isnan(p)) = 0; end
                        smoothp = runmean(p,round(fs/4));
                        oi = peakfinder(smoothp(round(I-5*fs):end),.5,6,-1); %find the first "peak" (surfacing) shallower than 6 m that is >.5m higher than surrounding areas.  Give 30 seconds leeway to ensure you hit it
                        oi = oi(oi>1); peakLoc(1) = oi(1)+I-5*fs-1; % greater than 1 to get rid of first peak if the whale is descending
                        isokay = false;
                        while ~isokay %double check calculated values
                            peakLoc(2:end) = peakLoc(1)+round((times(surfI(2:end),1)-times(surfI(1),1))*fs);
                            peakLoc = round(peakLoc);
                            viewI = max(1,round(I-30*fs)); % -times(surfI(1),1)*fs
                            viewI = viewI:min(round(viewI+vidDurs(i)*fs+30*fs),length(p));
                            figure(i); clf; plot(DN(peakLoc),min(p(viewI))+.05,'rs','markersize',18,'markerfacecolor','r','markeredgecolor','k');hold on;  set(i,'windowstyle','docked');
                            plot(DN(viewI),p(viewI));
                            set(gca,'ydir','rev','xlim',[DN(viewI(1)) DN(viewI(end))]);
                            oi = datestr(get(gca,'xtick'),'HH:MM:SS'); set(gca,'xticklabel',oi);
                            title(['Video Number ' num2str(i)]);
                            text(DN(viewI(1)),.9*max(p(viewI)),'If boxes correspond to surfacings, press enter, else click on the 1st surfacing that corresponds to the recorded value','verticalalignment','bottom','horizontalalignment','left','fontsize',14);
                            text(DN(viewI(1)),.92*max(p(viewI)),'Right click to zoom in','verticalalignment','top','horizontalalignment','left','fontsize',12);
                            %                     text(DN(viewI(1)),max(p(viewI)),'Press C to try to calibrate camera times (stupid FroBack)','verticalalignment','bottom','horizontalalignment','left','fontsize',14);
                            [I,~,button] = ginput(1);
                            if button == 3;
                                xlim([I-3/60/24 I+3/60/24]); oi = datestr(get(gca,'xtick'),'HH:MM:SS'); set(gca,'xticklabel',oi);
                                [I,~] = ginput(1);
                            end
                            %                     if button == 99;
                            %                        text(DN(viewI(1)),max(p(viewI)),'Press C to try to calibrate camera times (stupid FroBack)','verticalalignment','bottom','horizontalalignment','left','fontsize',14);
                            %
                            %                     end
                            if isempty(I); isokay = true;
                            else [~,I] = min(abs(DN-I)); peakLoc(1) = I; end
                        end
                    end
                    adjust = nan(size(surfI));
                    if any(diff(times(surfI,:),[],2))>10 || mean(diff(times(surfI,:),[],2))>5; longsurfacing = true; else longsurfacing = false; end
                    
                    % new to account for sleeping whales
                    if longsurfacing %if the surfacings are more than 5 seconds due to sleeping etc.
                        for j = 1:length(surfI)
                            [~,oi] = min(smoothp(peakLoc(j)-3*fs:peakLoc(j)+3*fs));
                            peakLoc(j) = round(oi+peakLoc(j)-3*fs-1); %adjust the peak to be the max of the smoothed data
                            surfT = round(fs*diff(times(surfI(j),:)));
                            [~,oi] = min(abs(smoothp(peakLoc(j)-surfT:peakLoc(j))-smoothp(peakLoc(j):peakLoc(j)+surfT)));
                            oi1 = find(smoothp(1:peakLoc(j))>smoothp(peakLoc(j))+0.5,1,'last');
                            oi2 = find(smoothp(peakLoc(j):end)>smoothp(peakLoc(j))+0.5,1,'first')+peakLoc(j)-1;
                            %                         adjust(j) = timeDN(surfI(j),1)+vidDN(i)-DN(peakLoc(j)-surfT+oi-1); %the DN difference that the video says compared to the data
                            oi1 = round(oi1+(oi2-oi1)/2-surfT/2); %start of surfacing
                            adjust(j) = timeDN(surfI(j),1) + vidDN(i) - DN(oi1);
                        end
                        
                        %old
                    else
                        for j = 1:length(surfI)
                            [~,oi] = min(smoothp(peakLoc(j)-3*fs:peakLoc(j)+3*fs));
                            peakLoc(j) = round(oi+peakLoc(j)-3*fs-1); %adjust the peak to be the max of the smoothed data
                            surfT = round(fs*diff(times(surfI(j),:)));
                            [~,oi] = min(abs(smoothp(peakLoc(j)-surfT:peakLoc(j))-smoothp(peakLoc(j):peakLoc(j)+surfT)));
                            adjust(j) = timeDN(surfI(j),1)+vidDN(i)-DN(peakLoc(j)-surfT+oi-1); %the DN difference that the video says compared to the data
                        end
                    end
                    for j = 1:length(manI)
                        adjust(end+1) = vidDN(i)+timeDN(manI(j),1)-times(manI(j),2);
                    end
                    %             disp(datestr(abs(adjust),'MM:SS.FFF'));
                    if any(abs(diff(adjust))>1/24/60/60);
                        SIGNS = sign(adjust);
                        ADJ = datestr(abs(adjust),':MM:SS.FFF');
                        for ii = 1:length(adjust); if SIGNS(ii) <0; ADJ(ii,1) = '-'; else ADJ(ii,1) = ' '; end; end
                        disp(['Adjustment indicated by each surfacing in video ' num2str(i)]);
                        disp(ADJ);
                        error('Surface time adjustments do not agree within 1 second');
                    end
                    disp('Adjustments (s):')
                    disp(num2str(round(adjust*24*60*60*1000)/1000));
                    if any(abs(adjust)>10/24/60/60) && i~=min(videonum); disp(['Large adjustment in video num ' num2str(i)]); if adjust<0; sn = num2str(sign(adjust)); sn = sn(1); else sn = ''; end; disp([sn datestr(abs(adjust),'MM:SS.FFF')]); end
                    adjust = mean(adjust);
                    disp(['Mean' num2str(i)]);
                    disp(datestr(abs(adjust),'MM:SS.FFF'));
                    vidDN(i) = vidDN(i)-adjust;
                    lastgood = i;
                end
                if ~isempty(noneI) || (~isempty(frameTimes{i}) && isempty(surfI) && isempty(manI)) % if there were no ways to calibrate, just calibrate based on video stamp difference between this video and the last
                    if any(abs(vidDN-vidDNorig)>0) % if there have been any adjustments, use them, else will have to use the next one.
                        vidDN(i) = vidDNorig(i)-vidDNorig(lastgood)+vidDN(lastgood);
                    else FIX = [FIX; i];
                    end
                end
            end
            %add in tag slips now that the video time has been adjusted
            if ~isempty(videonum)
                slipI = find(videonum == i & strcmp(types,'Slip'));
                for j = 1:length(slipI)
                    if isnan(times(slipI(j),1))
                        sI = vidDN(videonum(find(videonum < i,1,'last')))+vidDurs(videonum(find(videonum < i,1,'last')))/24/60/60; % okay to do i -1 because you shouldn't note a slip if it's the first video
%                         tagslipC = [tagslipC 0];
                        eI = vidDN(i);
                    else
                        if strcmp(headers{5,1},'Time Style') && strcmp(headers{5,2},'Embedded')
                            if ~synchusingvidtimestamps
                                times(slipI(j),:) = timeDN(slipI(j),:) + floor(vidDN(i));
                            else
                                times(slipI(j),:) = times(slipI(j),:) + floor(vidDN(i));
                            end
                            sI = times(slipI(j),1); %   tagslipC = [tagslipC 1];
                            eI = times(slipI(j),2);
                        else
                            sI = times(slipI(j),1)+vidDN(i); %   tagslipC = [tagslipC 1];
                            eI = times(slipI(j),2)+vidDN(i);
                        end
                    end
                    [~,sI] = min(abs(DN-sI)); [~,eI] = min(abs(DN-eI));
%                     [~,oi] = max(njerkTemp(sI:eI));
%                     tagslip = [tagslip oi+sI-1];
                    tagslip = [tagslip; [sI eI]];
                end
            end
        end
        if  ~synchusingvidtimestamps
            for ii = 1:length(FIX)
                i = FIX(ii);
                ADJ = find(abs(vidDN(i+1:end)-vidDNorig(i+1:end))>0,1,'first')+i; % if there have been any adjustments, use them, else will have to use the next one.
                vidDN(i) = vidDNorig(i)-vidDNorig(ADJ)+vidDN(ADJ);
                adjust = -vidDNorig(ADJ)+vidDN(ADJ);
                disp(['Adj for ' num2str(i)]);
                disp(datestr(abs(adjust),'MM:SS.FFF'));
            end
        end
        
        camon = false(size(p)); audon = false(size(p));
        for i = 1:length(vidDN)
            if isempty(frameTimes{i}) 
                if isnan(vidDN(i)); continue; end
                [~,sI] = min(abs(DN-vidDN(i)));
                [~,eI] = min(abs(DN-(vidDN(i)+vidDurs(i)/60/60/24))); %round(sI + vidDurs(i)*fs);
                audon(sI:eI) = true; %min(eI,length(p))) = true;
                continue;
            end
            [~,sI] = min(abs(DN-vidDN(i)));
            [~,eI] = min(abs(DN-(vidDN(i)+vidDurs(i)/60/60/24))); %round(sI + vidDurs(i)*fs);
            camon(sI:eI) = true; %min(eI,length(p))) = true;
        end
        if ~isempty(audstart) && ~isnan(audstart)
            [~,I] = min(abs(DN-audstart));
            audon(I:end) = true;
        end
        
    if any(isnan(vidDN(~cellfun(@isempty,vidNam)))) || any(diff(vidDN(~cellfun(@isempty,vidNam)))<0) || any(isnan(vidDurs(~cellfun(@isempty,vidNam))))
        error('Check vidDN and vidDurs!  Nans found.  Likely there are video files in the raw folder that are before or after the deployment, so you may need to open "movieTimes" in a different matlab version and add nans or empty cells to where the movies that are not part of the deployment should be.  Another thing to check is "videonum".  If "videonum" is longer than your # of videos, you may need to adjust your excel header files to ensure that any blanks spaces are actually deleted.');
    end
end
vidadj = vidDN-vidDNorig;
