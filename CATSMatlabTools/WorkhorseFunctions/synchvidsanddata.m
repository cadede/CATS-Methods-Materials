function [camon,audon,vidDN,nocam,tagslip] = synchvidsanddata(data,headers,Hzs,DN,ODN,fs,CAL,synchusingvidtimestamps,useFrames)
global fileloc filename
if nargin<10; useFrames = false; end %this is a legacy switch for if you enter framenumbers into the excel sheet instead of times
if sum(diff(data.Pressure)<.001) == length(data.Pressure); nopress = true; else nopress = false; end

try if sum(data.Camera) == 0; nocam = true; else nocam = false; end; catch; disp('Could not automatically detect camera status'); end
 times = cell2mat(headers(7:end,2:3));
%if no movie files
if ((size(headers,1)<7 || (isnan(headers{7,1})||headers{7,1} == 0))) || nocam
%     GPS = cell2mat(headers(2,2:3)); %from above file
%     whaleName = char(headers(1,2));
%     timedif = cell2mat(headers(3,2)); % The number of hours the tag time is behind (if it's in a time zone ahead, use -).  Also account for day differences here (as 24*# of days ahead)
%     DN = data.Date+data.Time+timedif/24;
%     DV = datevec(DN);
    tagslip = [1 1]; %tagslipC = 1; % confidence of tagslip.  1 if you see it move on a video, 0 if you estimate based on max jerk
    camon = false(size(DN)); audon = false(size(DN));
    if exist('pconst','var'); p = (data.Pressure-pconst)*pcal;
        %     [~,peakmag] = peakfinder(p,3,min(p)+6,-1); p = p-max(min(peakmag),max(peakmag)-2);
    else p = data.Pressure;
    end
%     try
%         p = surfadj(p,fs,min(p)+4,min(p)+2, ceil(nanmean(diff(times,[],2))/2)); % uses Ann's surface adjustment to "0 out" the pressure sensor;
%         disp('Used Anne''s method');
%     catch %if times does not exist
%         try
%             times = [1 4];
%             p = surfadj(p,fs,min(p)+4,min(p)+2, ceil(nanmean(diff(times,[],2))/2)); % use
%             disp('Used Anne''s method');
%         catch
%             p = Depth_Correction(p,fs);
%             disp('Used Angie''s method');
%         end
%     end
    
%       nocam = true;
%     if ~exist('magcalon','var'); magcalon = magcaloff; magconston = magconstoff; end
else
   
    useold = false;
    nocam = false;
    load([fileloc filename(1:end-4) 'movieTimes.mat']); %load frameTimes and videoDur from the movies, as well as any previously determined info from previous prh makings with different decimation factors
%     if (tagnum >= 40 && tagnum < 50)||tagnum>51 
%         kitten = true; else kitten = false;
%     end
%     kitten = false;
%     if isempty(strfind(movies{1},'CATS')) && isempty(strfind(movies{1},'AW')) && isempty(strfind(movies{1},'MP4'))
%         kitten = true; else kitten = false;
%     end
    if synchusingvidtimestamps;
%         try Hzs = importdata([fileloc filename(1:end-3) 'txt']);
%         catch; try Hzs = importdata([fileloc 'raw\' filename(1:end-3) 'txt']);
%             catch; try Hzs = importdata([fileloc 'raw\' filename(1:end-5) '.txt']); catch; try Hzs = importdata([fileloc filename(1:end-5) '.txt']);
%                     catch; try Hzs = importdata([fileloc(1:end-2) 'raw\' filename(1:end-3) 'txt']); catch; Hzs = importdata([fileloc(1:end-2) 'raw\' filename(1:end-5) '.txt']);
%                         end;end;end;end
%         end; Hzs = Hzs.textdata;
%         UTC = Hzs{find(~cellfun(@isempty, strfind(Hzs,'utc')),1,'first')};
%         UTC = str2num(UTC(regexp(UTC,'=')+1:end));
        UTC = Hzs.UTC;
        timedif = cell2mat(headers(3,2));
        vidDNorig = vidDN;
        if data.Date(1)<datenum([2017 09 01]) % old versions of tag where time stamps were in UTC
            vidDN = vidDN + data.Date(1)+UTC/24+timedif/24;
        elseif vidDN(find(~isnan(vidDN),1))<10 % confirm that new vidDN takes into account full date
            vidDN = vidDN + data.Date(1)+timedif/24;
        else vidDN = vidDN + timedif/24;
        end
%         vidDN = vidDN-timedif/24;
        if sum(vidDN<DN(1)) + sum(vidDN>DN(end))>1; warning('multiple vidDN are outside range of data. This may be okay if more videos are included in movieTimes outside of those on whale (so just continue), but check that any timedif in xls header file should apply to both data and video.  If, for example, data were downloaded on a computer in a different timezone than for which the tag was programmed, video and data will be offset by different amounts and will have to be adjused for, potentially by unhighlighting the line above this one'); end
    end
    
    if useold; try load([fileloc filename(1:end-4) 'Info.mat']); catch; end; end
    if isempty (frameTimes); error('Put frametimes from movies in the same folder as the movies'); end
    % if ~(exist('GPS','var') && exist('vidDN','var')&&exist('camon','var') && exist('tagslip','var')) % if you've saved previous versions in the Info file, just use those
    if useold && exist('vidDN','var') && exist('camon','var')
        videonum = cell2mat(headers(7:end,1)); todel = find(~isnan(videonum),1,'last')+1:length(videonum);
        videonum(todel) = [];
    else
%         GPS = cell2mat(headers(2,2:3)); %from above file
        whaleName = char(headers(1,2));
        timedif = cell2mat(headers(3,2)); % The number of hours the tag time is behind (if it's in a time zone ahead, use -).  Also account for day differences here (as 24*# of days ahead)
        try
            videonum = cell2mat(headers(7:end,1)); todel = find(~isnan(videonum),1,'last')+1:length(videonum);
            types = headers(7:end,4);
            times = cell2mat(headers(7:end,2:3));
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
%             try
%                 p = surfadj(p,fs,min(p)+4,min(p)+2, ceil(nanmean(diff(times,[],2))/2)); % uses Ann's surface adjustment to "0 out" the pressure sensor;
%                 disp('Used Anne''s method');
%             catch %if times does not exist
%                 try
%                     times2 = [1 4];
%                     p = surfadj(p,fs,min(p)+4,min(p)+2, ceil(nanmean(diff(times2,[],2))/2)); % use
%                     disp('Used Anne''s method');
%                 catch
%                     p = Depth_Correction(p,fs);
%                     disp('Used Angie''s method');
%                 end
%             end
        end
        DN = data.Date+data.Time+timedif/24;
%         DV = datevec(DN);
        %adjust the video date numbers by assuming video 1 started with the start
        %of the tag
        ODNa = ODN + timedif/24;
        if ~synchusingvidtimestamps
            if ~isnan(vidDN(1))
                vidDN = vidDN - (vidDN(1)-(ODNa-floor(ODNa))+5/24/60/60); %assumes the 1st video starts 5 seconds after you turn on the data
            end
            vidDN = vidDN +floor(DN(1)); % account for the date
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
        Atemp = (filterCATS([data.Acc1 data.Acc2 data.Acc3],ceil(fs/8),round(fs),.05)-repmat(CAL.aconst,size(data,1),1))*CAL.acal; %temp Acceleration file for guessing at tagslip location
        njerkTemp = (9.81*fs)*sqrt(diff(Atemp).^2*ones(3,1)); %temp jerk for guessing at tag slip
        
        timeDN = times/24/60/60; %times was in seconds
        lastgood = 1; FIX = [];
        for i = 1:length(frameTimes)
            if ~synchusingvidtimestamps
                surfI = find(videonum == i & strcmp(types,'Surface'));
                manI = find(videonum == i & strcmp(types,'Manual')); % if the whale does not surface, can do a manual calibration using another factor (tagoff, no paddlewheel to paddlewheel audio etc.)
                noneI = find(videonum == i & strcmp(types,'None'));
                if isempty(frameTimes{i}); continue; end
                if (~isempty(surfI) || ~isempty(manI))&&(~synchusingvidtimestamps)
                    if ~isempty(surfI)
                        if i == min(videonum) || i == videonum(find(cellfun(@(x) strcmp(x,'Surface'),types),1)); I = round(find(p>1,1,'first')+5*fs); %time for the camera to turn on after the camera is deployed
                        else [~,I] = min(abs(DN-(vidDN(i)+timeDN(surfI(1),1)))); end
                        peakLoc = nan(size(surfI));
                        smoothp = runmean(p,round(fs/4));
                        oi = peakfinder(smoothp(round(I-5*fs):end),.5,6,-1); %find the first "peak" (surfacing) shallower than 6 m that is >.5m higher than surrounding areas.  Give 5 seconds leeway to ensure you hit it
                        oi = oi(oi>1); peakLoc(1) = oi(1)+I-5*fs-1; % greater than 1 to get rid of first peak if the whale is descending
                        isokay = false;
                        while ~isokay %double check calculated values
                            peakLoc(2:end) = peakLoc(1)+round((times(surfI(2:end),1)-times(surfI(1),1))*fs);
                            peakLoc = round(peakLoc);
                            viewI = max(1,round(I-30*fs-times(surfI(1),1)*fs));
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
                    disp('Adjustments:')
                    disp(datestr(adjust,'HH:MM:SS.FFF'));
                    if any(abs(adjust)>10/24/60/60) && i~=min(videonum); disp(['Large adjustment in video num ' num2str(i)]); if adjust<0; sn = num2str(sign(adjust)); sn = sn(1); else sn = ''; end; disp([sn datestr(abs(adjust),'MM:SS.FFF')]); end
                    adjust = mean(adjust);
                    disp(['Mean' num2str(i)]);
                    disp(datestr(abs(adjust),'MM:SS.FFF'));
                    vidDN(i:end) = vidDN(i:end)-adjust;
                    lastgood = i;
                end
                if ~isempty(noneI) % if there were no ways to calibrate, just calibrate based on video stamp difference between this video and the last
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
                    else sI = times(slipI(j),1)+vidDN(i); %   tagslipC = [tagslipC 1];
                        eI = times(slipI(j),2)+vidDN(i);
                    end
                    [~,sI] = min(abs(DN-sI)); [~,eI] = min(abs(DN-eI));
%                     [~,oi] = max(njerkTemp(sI:eI));
%                     tagslip = [tagslip oi+sI-1];
                    tagslip = [tagslip; [sI eI]];
                end
            end
        end
        if  ~synchusingvidtimestamps %~kitten ||
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
        
    end
    if any(isnan(vidDN(~cellfun(@isempty,vidNam)))) || any(diff(vidDN(~cellfun(@isempty,vidNam)))<0) || any(isnan(vidDurs(~cellfun(@isempty,vidNam))))
        error('Check vidDN and vidDurs!  Nans found.  Likely there are video files in the raw folder that are before or after the deployment, so you may need to open "movieTimes" in a different matlab version and add nans or empty cells to where the movies that are not part of the deployment should be.  Another thing to check is "videonum".  If "videonum" is longer than your # of videos, you may need to adjust your excel header files to ensure that any blanks spaces are actually deleted.');
    end
end
