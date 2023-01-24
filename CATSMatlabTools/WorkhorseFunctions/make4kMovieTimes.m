function make4kMovieTimes(dur,folder,ripAudio,timestamps,whaleID)
dbstop if error;

% this function reads the time off each video frame for high resolution
% data.  May also work with lower resolution videos but has not been tested

% timestamps- signals whether to read embedded timestamps on the video or
% just read the encoded timestamps (i.e. if timestamps do not exist )
% folder - optional.  Can set a default folder to select a tag guide or check for videos.
% ripAudio - optional.  Default is true.  Rips the audio from video files
% into a separate folder (AudioData).  Set to false if you've already done
% this and just want to get the video time stamps (e.g. if the process got
% interupted and needed to be restarted).
% vidNums - optional.  if not present, looks at all videos.  If it is
% present, look for an existing movieTimes file and fill in the videos listed in "vidnums"
% audioonly- optional, reads audio start times from wav files and converts
% to audio data if no videos exist

% if nargin<6 || isempty(timewarn); timewarn = 0.1; end
if nargin<3 || isempty(ripAudio); ripAudio = false; end
if nargin<2 || isempty(folder);  folder = 'E:\CATS\'; end
if nargin<1 || isempty(dur); dur = 15; end
if ~exist('audioonly','var') || isempty(audioonly); audioonly = false; end

cf = pwd; try cd(folder); catch; end
titl = 'select FIRST movie file from deployment';
  if ~ispc; menu(titl,'OK'); end
[movies,movieloc] = uigetfile('*.*',titl,'multiselect','off');
cd(movieloc);
titl = 'select .bin file (just for the name)';
  if ~ispc; menu(titl,'OK'); end
[datafile,dataloc] = uigetfile('*.bin',titl,'multiselect','off');
cd(cf);
if ischar(movies); movies = {movies}; end
m1 = 1;

if ~exist('whaleID','var') || isempty(whaleID)
    whaleID = [];
    slashes = union(regexp(movieloc,'\'),regexp(movieloc,'/'));
    for i = 1:length(slashes)
        try
            if strcmp(movieloc(slashes(i)+9),'-')
                datenum(movieloc(slashes(i)+3:slashes(i)+8),'yymmdd');
                whaleID = movieloc(slashes(i)+1:min(slashes(i)+11,slashes(i+1)-1));
            end
        catch
            continue
        end
    end
end
%
disp(['ID = ' whaleID]);
% grab the movie files from the directory
m2 = dir(movieloc); m2 = {m2.name}; todel = false(size(m2)); todel(1:2) = true; for i = 3:length(m2); if length(m2{i})<=3 || ~(strcmp(m2{i}(end-2:end),'raw')||(strcmp(m2{i}(end-2:end),'wav')&&audioonly)||strcmpi(m2{i}(end-2:end),'mov')||strcmpi(m2{i}(end-2:end),'mp4')); todel(i) = true; end; end; m2(todel) = [];
mfirst = find(cellfun(@(x) strcmp(x,movies{1}), m2));
% these try catch help choose which files to import
% try % find all movies in the water to get all relevant audio
%     m2times = nan(size(m2));
%     for i = 1:length(m2); m2times(i) = datenum(m2{i}(min(regexp(m2{i},'-'))+1:max(regexp(m2{i},'-'))-1),'yyyymmdd-HHMMSS'); end
%     rootDIR = folder; loc1 = strfind(movieloc,'CATS'); if ~isempty(loc1); rootDIR = movieloc(1:loc1+4); end
%     try [~,~,txt] = xlsread([rootDIR 'TAG GUIDE.xlsx']); catch;  titl = 'Get Tag Guide to find tag off and recovery times'; if ~ispc; menu(titl,'OK'); end; [filename, fileloc] = uigetfile('*.xls*',titl); [~,~,txt] = xlsread([fileloc filename]); end
%     rows = find(~cellfun(@isempty, cellfun(@(x) strfind(x,whaleID),txt(:,1),'uniformoutput',false)));
%     if isempty(rows); error('Could not find whaleID in TAG GUIDE'); end
%     if length(rows)>1; let = input('letter of interest? (1 = a, 2 = b, 3 = c, etc.) '); lets = 'abcdefghijklmnopqrstuvwxyz'; whaleID = [whaleID lets(let)]; rows = find(~cellfun(@isempty, cellfun(@(x) strfind(x,whaleID),txt(:,1),'uniformoutput',false))); end
%     col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Recover_Time'),txt(3,:),'uniformoutput',false)));
%     RecTime = max(datenum(txt(rows,col)));
%     m1 = find(cellfun(@(x) strcmp(x,movies{1}), m2));
%     mlast = find(m2times<RecTime,1,'last');
%     movies = m2(m1:mlast); 
%     disp(['Got last audio file from Tag Recovery Time, last file in water is ' movies{end}]);
% catch
%     cd(movieloc);
%     titl = 'No recovery time in Tag Guide, select last movie or audio before recovery (to get all wav files)';
%     if ~ispc; menu(titl,'OK'); end
%     [mlast,movieloc] = uigetfile('*.*',titl,'multiselect','off');
%     m1 = find(cellfun(@(x) strcmpi(x,movies{1}), m2));
%     mlast = find(cellfun(@(x) strcmpi(x,mlast), m2));
%     movies = m2(m1:mlast);
%     disp(['last file selected in water is ' movies{end}]);
% end

try % find all movies on whale to get frameTimes
    m2times = nan(size(m2));
    for i = 1:length(m2); m2times(i) = datenum(m2{i}(min(regexp(m2{i},'-'))+1:max(regexp(m2{i},'-'))-1),'yyyymmdd-HHMMSS'); end
    rootDIR = strfind(movieloc,'CATS'); rootDIR = movieloc(1:rootDIR+4);
    try [~,~,txt] = xlsread([rootDIR 'TAG GUIDE.xlsx']); catch; try [~,~,txt] = xlsread([fileloc filename]); catch; cd(fileloc); titl = 'Get Tag Guide to find tag off and recovery times'; if ~ispc; menu(titl,'OK'); end; [filename, fileloc] = uigetfile('*.xls*',titl); [~,~,txt] = xlsread([fileloc filename]); end; end
    rows = find(~cellfun(@isempty, cellfun(@(x) strfind(x,whaleID),txt(:,1),'uniformoutput',false)));
    col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Tag_Off'),txt(3,:),'uniformoutput',false)));
    RecTime = max(datenum(txt(rows,col)));
    mlast = find(m2times<RecTime,1,'last');
    movies = m2(mfirst:mlast); 
    LastOnMovie = m2{mlast};    
    disp(['Got last movie file from Tag Off time, last file on whale is ' LastOnMovie]);
catch
    cd(movieloc);
    titl = 'No tag off time in Tag Guide, select last movie or audio before tag fell off whale (to get movie times)';
    if ~ispc; menu(titl,'OK'); end; 
    [mlast,movieloc] = uigetfile('*.*',titl,'multiselect','off');
    mlast = find(cellfun(@(x) strcmp(x,mlast), m2));
    LastOnMovie = m2{mlast};
    movies = m2(mfirst:mlast); 
    disp(['last file selected on whale is ' LastOnMovie]);
end

% this is CATS specific, read the video numbers off the video names
movN = nan(size(movies));
for n = 1:length(movies)
    lastnum = regexpi(movies{n},'.mp4');
    if isempty(lastnum); lastnum = regexpi(movies{n},'.mov'); end
    if isempty(lastnum); lastnum = regexpi(movies{n},'.raw'); end
     if isempty(lastnum)&&audioonly; lastnum = regexpi(movies{n},'.wav'); end
    movN(n) = str2num(movies{n}(lastnum-3:lastnum-1));
%     else audN(n) = str2num(movies{n}(lastnum-3:lastnum-1));
%     end
end
mfirst = movN(1);


% One version of CATS video files did not have video numbers (only times).  This adds them in. (Should be unnecessary now).
% if any(diff(movN)<1) || any (diff(movN)>40)
%     rename = input('Rename files to include video numbers? (Y or N) ','s');
%     if regexpi(rename,'y')
%         D = dir([movieloc '*.mp4']); % get the total lengths of the other files as well as the stored file time
%         if isempty(D); D = dir([movieloc '*.mov']); end
%         vidNam = cellstr(vertcat(D.name));
%         for n = 1:length(vidNam)
%             VN = [vidNam{n}(1:end-4) '-' sprintf('%08.f',n) vidNam{n}(end-3:end)];
%            copyfile([movieloc vidNam{n}],[movieloc VN]);
%            delete([movieloc vidNam{n}]);
%            vidNam{n} = VN;
%         end
%         for n = 1:length(movies)
%             movies{n} = vidNam{~cellfun(@isempty,cellfun(@(x) strfind(x,movies{n}(1:end-4),vidNam),'uniformoutput',false))};
%             lastnum = regexpi(movies{n},'.mp4');
%             if isempty(lastnum); lastnum = regexpi(movies{n},'.mov'); end
%             movN(n) = str2num(movies{n}(lastnum-3:lastnum-1));
%         end
%     end
% end

frameTimes = cell(max(movN),1);
vidDurs = nan(size(frameTimes));
if strcmp(dataloc(end-3:end-1),'raw'); dataloc = dataloc(1:end-4); end
DIR = [dataloc 'AudioData\'];
if ~exist (DIR,'dir')
    mkdir (DIR);
end

% if timestamps; oframeTimes = frameTimes; end

warning('off','all');
%
if  ripAudio 
    % see if there are any wavfiles
    [D,F] = subdir(dataloc);
    wavfiles = cell(0,1);
    wavstr = wavfiles;
    D1 = dir(dataloc); id = [D1.isdir]; D1 = {D1.name}'; D1 = D1(~id);
    for i = 1:length(D1); D1{i} = [movieloc D1{i}]; end
    disp('Reading audio files from .raw and wav files: ');
    for i = 1:length(D); for j = 1:length(F{i}); D1{end+1} = [D{i} '\' F{i}{j}]; if strcmp(F{i}{j}(end-3:end),'.wav'); wavfiles{end+1} = F{i}{j}; wavstr{end+1} = D1{end}; disp(wavstr{end}); end; end; end
    for i = 1:length(movies); if strcmp(movies{i}(end-2:end),'raw'); disp([movieloc movies{i}]); end; end
    readaudiofiles2;
    disp(['Check that audio files were successfully written, then if you''re sure, you can delete files after ' LastOnMovie]);
    disp('Any premade wav files were read and left in their original directory');
else
    disp('no audio was extracted from vids');
end
warning('on','all');
%
D = dir([movieloc '*.mp4']); % get the total lengths of the other files as well as the stored file time
if isempty(D); D = dir([movieloc '*.mov']); end
if isempty(D) && audioonly; D = dir([movieloc '*.wav']); 
   if isempty(D) && audioonly; D = dir([DIR '*.wav']); end
end
% some other code if you need to restrict what movies it finds
% todel = [];
% for i = 1:length(D); if length(D(i).name) ~= 10; todel = [todel; i]; end; end % in case you've made any SxS files already or something silly
% D(todel) = [];
% D2 = dir([movieloc 'CATS*.mp4']);  todel = [];
% for i = 1:length(D2); if length(D2(i).name) ~= 12; todel = [todel; i]; end; end % in case you've made any SxS files already or something silly
% D2(todel) = [];
% D = [D; D2]; % new 2ks have longer file names and start with CATS
vidDN = vertcat(D.datenum); vidDN = vidDN - floor(vidDN(1)); % dates don't match, just worry about the time (but this way keeps track of if the video goes to the next day)


k = 1; vidNam = {}; vidNum = [];
for i = 1:length(movies)
    if ~strcmp(movies{i}(end-2:end),'raw') && ~strcmp(movies{i}(end-2:end),'wav')
%         vid = mmread([movieloc movies{i}],[1 10]); % just read a few frames to get the total duration
        try vid = VideoReader([movieloc movies{i}]);
            
            if isnan(vidDurs(movN(i)))
                vidDurs(movN(i)) = vid.Duration;
            end
        catch;  vid = mmread([movieloc movies{i}],[1 10]); % just read a few frames to get the total duration, this was for old style videos, and seems to work for some 4ks
            if isnan(vidDurs(movN(i)))
                vidDurs(movN(i)) = vid.totalDuration;
            end
        end
        vidNam(k,1) = movies(i); 
        vidNum(k,1) = str2num(vidNam{k}(end-7:end-4)); k = k+1;
    else
        load([DIR movies{i}(1:end-4) 'audio.mat'],'totalDuration');
        vidDurs(movN(i)) = totalDuration;
        if audioonly
            vidNam(k,1) = movies(i);
            vidNum(k,1) = str2num(vidNam{k}(end-7:end-4)); k = k+1;
        end
    end
end
oi = nan(size(vidDurs));
% if ~timestamps; oi(vidNum) = vidDN(ismember(cellstr(vertcat(D.name)),movies)); end
vidDN = oi;
oi = cell(size(vidDurs));
oi(vidNum) = vidNam;
vidNam = oi;

% may need to restart here if there was an error somewhere along the way.
% Next line recreates an empty frameTimes and original frameTimes
% (oframeTimes).  oframeTimes is the time embedded in the frame metadata (the output from mmread,
% but if there were any skips in data this can be off, so frameTimes reads
% the timestamp written on the frame.
% frameTimes = cell(max(movN),1); oframeTimes = frameTimes;
if exist('vidNums','var') && ~isempty(vidNums) % if you signaled to only read a couple of the videos, load all the videos first
    try load([dataloc datafile(1:end-4) 'movieTimes.mat'],'frameTimes','oframeTimes','vidDN','vidDurs'); disp('existing movietimes file loaded');
    catch
        try
            load([dataloc datafile(1:end-4) 'movieTimesTEMP.mat'],'frameTimes','oframeTimes','vidDN','vidDurs');
            disp('movieTimesTEMP file loaded');
        catch
            disp('WARNING: No old frameTimes found, resulting frameTimes file will only be for indicated videos'); 
        end
    end
    for n = 1:length(vidNums)
        frameTimes{vidNums(n)} = []; oframeTimes{vidNums(n)} = [];
    end
else vidNums = movN;
end
% 
% start reading frames here:
 badvidDN = false(size(vidDN)); viddifs = zeros(size(badvidDN));
for n = 1:length(movies) 
    if isempty(intersect(movN(n),vidNums)); continue; end
    if audioonly; vidDN(movN(n)) = datenum(movies{n}(min(regexp(movies{n},'-'))+1:max(regexp(movies{n},'-'))-1),'yyyymmdd-HHMMSS-fff'); 
        try if abs(vidDN(movN(n)) - (vidDN(movN(n)-1) + vidDurs(movN(n)-1)/24/60/60))>60/24/60/60; disp(['WARNING!: audio ' num2str(movN(n)) ' appears to start ' num2str(vidDN(movN(n)) - (vidDN(movN(n)-1) + vidDurs(movN(n)-1)/24/60/60)) ' s from the end of the previous file, perhaps check your audio rate?']); end; catch; end
        continue; 
    end
%     try 
        vidDN(movN(n)) = datenum(movies{n}(min(regexp(movies{n},'-'))+1:max(regexp(movies{n},'-'))-1),'yyyymmdd-HHMMSS-fff');
%     catch; warning(['Cannot read precise video start time from video ' num2str(movN(n)) ', will try to read video from timestamps on video, else may need to adjust manually'])
%         badvidDN(movN(n)) = true;
%          d = regexp(movies{n},'-');
%         if ~strcmp(movies{n}(end-2:end),'raw')
%             day = datenum(movies{n}(d(1)+1:d(2)-1),'yyyymmdd');
%             starttime = 0; endtime = dur; videoL = dur+1; DAY = 0; badmovie = false;
%             readwirelessvideo2;
%             oi = checkbadframes(vid.times);
%             if any(diff(oi)<0 | diff(oi)>2.5*median(diff(oi))); dur0 = dur; dur = round(dur/3); readwirelessvideo2; oi = checkbadframes(vid.times); dur = dur0; end
%             if ~any(diff(oi)<0 | diff(oi)>1/24/60/60) % if a gap is bigger than 1 second, try to fix it later2.5*median(diff(oi)));
%                 vidDN(movN(n)) = day + vid.times(1);
%                 disp(['Video ' num2str(movN(n)) ' calculated to start at ' datestr(vidDN(movN(n)),'yyyy-mmm-dd HH:MM:SS.fff')]);
%                 badvidDN(movN(n)) = false;
%             else
%                  try vidDN(movN(n)) = datenum(movies{n}(d(1)+1:d(3)-1),'yyyymmdd-HHMMSS');
%                      warning(['movie file ' num2str(movN(n)) ' start time could only be read to the nearest second: ' datestr(vidDN(movN(n)))]);
%                      disp('RECOMMEND synching movie times via surfacings');
%                      badvidDN(movN(n)) = false;
%                  catch; error('No video start time, recommend avoiding simpleread')
%                  end
%             end
%         else
%             try vidDN(movN(n)) = datenum(movies{n}(d(1)+1:d(3)-1),'yyyymmdd-HHMMSS');
%             disp(['Audio file ' num2str(movN(n)) ' start time could only be read to the nearest second: ' datestr(vidDN(movN(n)))]);
%             badvidDN(movN(n)) = false;
%             catch; warning('Could not read audio time stamp');
%             end
%         end
%     end
%     if strcmp(movies{n}(end-2:end),'raw') || n >  mlast-m1+1
%         vidNam{movN(n)} = movies{n};
%        if movN(n)>movN(mlast-m1+1); vidNam{movN(n)} = []; vidDurs(movN(n)) = nan; end
%        continue; 
%     end
    starttime = 0;
    endtime = dur;
    videoL = dur+1;
    DAY = 0;
    badmovie = false;
    try
    OBJ =  VideoReader([movieloc movies{n}]);
    fr = OBJ.FrameRate;
     catch;  vid = mmread([movieloc movies{i}],[1 10]); % just read a 
         fr = vid.rate;
    end
frameTimes{vidNum(n)} = 1/fr:1/fr:OBJ.Duration;
   
%     while endtime<videoL+dur;
%         clear M vid aud;
%          vid = struct();
%         vid.times = starttime+1/fr:1/fr:endtime;
%         vid = read(OBJ,[starttime*fr+1 endtime*fr]);
        
%         if ~timestamps
%             [vid,~] = mmreadFT([movieloc movies{n}], [],[starttime endtime],false,true); %FT saves some time, but still not as efficient as it should be since it reads the video and not just the time
%         elseif ~simpleread
%             readwirelessvideo2; 
%             if badmovie; if ~exist([movieloc 'bad movies\'], 'dir'); mkdir([movieloc 'bad movies\']); end
%                 movefile([movieloc movies{n}], [movieloc 'bad movies\' movies{n}]);
%                 frameTimes{vidNum(n)} = []; oframeTimes{vidNum(n)} = []; vidDN(vidNum(n)) = nan;
%             end
%         else %simpleread and timestamps.  This only reads some of the timecodes from the video frame, and compares it to what is read from mmread
%             vid = mmread([movieloc movies{n}], [],[starttime endtime],false,true);
%             flag = true; iii = -1;
%             while flag && iii<100 % try 100 frames until you get a good read
%                 iii = iii+1;
%                 [newf, flag] = gettimestamp(vid.frames(end-iii).cdata,false);
%             end
%             if ~flag; viddif = newf*24*60*60-vid.times(end-iii)-(vidDN(movN(n))-floor(vidDN(movN(n))))*24*60*60; else viddif = nan; end
%             if abs(viddif)>timewarn || isnan(viddif);
%                 warning(['end frame of this snippet appears to be ' num2str(viddif) ' s off from video''s time stamp']);
%             end
%             % this resets oframeTimes in case some vid.times had been added
%             % during readwireless2 (which was run to check the initial time
%             % stamp)
%             if starttime == 0; oframeTimes{movN(n)} = frameTimes{movN(n)}; end
%             oframeTimes{movN(n)} = [oframeTimes{movN(n)} vid.times];
%         end
    
%         if ~badmovie; 
%             frameTimes{movN(n)} = [frameTimes{movN(n)} vid.times];
% %     end
% %         if sum(isnan(vid.times)) > 0; lbl = [' with ' num2str(sum(isnan(vid.times))) ' bad data points (will be fixed at movie end)']; else lbl = ''; end
% %         disp([num2str(endtime) ' sec completed' lbl]);
% %         videoL = OBJ.Duration;
%         if timestamps && ~simpleread; starttime = endtime-0.1; else starttime = endtime; end
%         endtime = endtime+dur;
%         if endtime+3>videoL; endtime = endtime+3;
%             disp('Added 3 seconds to second to last partial read (last partial read would have been very short');
%         end % if the final clip would be < 3 seconds, just incorporate it into the 2nd to last clip
%     end
%     if ~simpleread
%         [frameTimes{movN(n)},numbad,numbadsect] = checkbadframes(frameTimes{movN(n)});
%         lbl = [' with ' num2str(numbad) ' frames that could not be read and were interpolated, and ' num2str(numbadsect) ' frames that were shifted and corrected'];
%     else if abs(viddif)>timewarn; lbl = [', end frame of movie appears to be ' num2str(viddif) ' s off from video''s time stamp- can choose to offset movie start time at end of process']; else lbl = [', and end of video is within time differential threshold (' num2str(viddif) 's)']; end
%         viddifs(movN(n)) = viddif;
%     end
%     
%     if timestamps && ~badmovie && ~simpleread
%         if badvidDN(movN(n))
%             vidDN(movN(n)) = frameTimes{movN(n)}(1);
%         end
%             frameTimes{movN(n)} = (frameTimes{movN(n)} - frameTimes{movN(n)}(1))*24*60*60;
%              if n>1 && vidDN(movN(n)) < vidDN(movN(n-1)); vidDN(movN(n)) = vidDN(movN(n))+1; DAY = DAY+1; end
% 
%         vidDurs(movN(n)) = frameTimes{movN(n)}(end)+median(diff(frameTimes{movN(n)}));
%     end
%     disp([movies{n} ' completed' lbl]);
%     save([dataloc datafile(1:end-4) 'movieTimesTEMP.mat'],'frameTimes','vidDN','videoL','vidDurs');
%     if timestamps; save([dataloc datafile(1:end-4) 'movieTimesTEMP.mat'],'oframeTimes','-append'); end
end
% 
% if ~timestamps && sum(badvidDN)>= length(movN)
%     disp('Using metadata timestamps from video files (times could not be read from video name)')
%     disp('Do these timestamps come from the beginning of the video (e.g. gopro?) or the end of the video (e.g. 2k)?');
%     s = input('1 = beginning, 2 = end');
%     if s == 2;
%         vidDN = vidDN - vidDurs/60/60/24;  % new CATS videos from 2k save date stamp at the END of recording rather than the beginning
%     end
% end
% 
% if timestamps && ~simpleread  && ~audioonly% if vid start time is read from the time stamps.
%     for n = 1:length(movies)
%         if isempty(intersect(movN(n),vidNums)); continue; end
%         if ~strcmp(movies{n}(end-2:end),'raw') && n<=  mlast-m1+1;
%             DAY = floor(datenum(movies{n}(min(regexp(movies{n},'-'))+1:max(regexp(movies{n},'-'))-1),'yyyymmdd-HHMMSS'));
%             vidDN(movN(n)) = DAY+vidDN(movN(n))-floor(vidDN(movN(n)));
%             continue;
%         end
%     end
% end
% 
% % make some plots to check if there are any potential errors in the frame
% % reads
if ~audioonly
% figure; hold on; I = 1; numrepeats = zeros(size(vidDN));
% title('Plot of frameTimes.  Red stars indicate a jump in time (if no red stars, you should be good)');
% for i = 1:length(frameTimes); 
%     if isempty(frameTimes{i}); continue; end
%     plot(I:I+length(frameTimes{i})-1,vidDN(i)+frameTimes{i}/24/60/60); 
%     if ~isempty(frameTimes{i}); text(I,frameTimes{i}(1)/24/60/60+vidDN(i),num2str(i)); end; 
%     oi = find(diff(frameTimes{i})>5*median(diff(frameTimes{i})) | diff(frameTimes{i}) < -median(diff(frameTimes{i})));% mark if it's a whole frame backward
%     oi2 = find(abs(diff(frameTimes{i})) < median(diff(frameTimes{i}))/2 | (diff(frameTimes{i}) < 0 & diff(frameTimes{i}) > -median(diff(frameTimes{i}))));
%     numrepeats(i) = length(oi2); 
%     % fix the very close frames that are less than 0 to be just bigger than 0 so that time
%     % never decreases.
%     oi2 = find((diff(frameTimes{i}) < 0 & diff(frameTimes{i}) > -median(diff(frameTimes{i}))));
%     frameTimes{i}(oi2+1) = frameTimes{i}(oi2)+.001;
%     if ~isempty(oi); 
%         plot(I+oi-1,vidDN(i)+frameTimes{i}(oi)/24/60/60,'r*','markersize',8); hold on;
%         df = round(1000*diff(frameTimes{i}))/1000;
%         for ii = 1:length(oi); text(I+oi(ii)-1,vidDN(i)+frameTimes{i}(oi(ii))/24/60/60,[num2str(df(oi(ii))) ' s jump at vidFrame ' num2str(oi(ii))],'rotation',90,'verticalalignment','bottom'); end
%         title('zoom in and check red stars, may have to check videos files and adjust frameTimes and vidDN if there is an error');
%     end
%     I = I+length(frameTimes{i});
%     if (~simpleread || badvidDN(i))&&~isnan(vidDN(i)); disp(['video # ' num2str(i) ' calculated to start at ' datestr(vidDN(i),'HH:MM:SS.fff')]); end
%     if simpleread && abs(viddifs(i))>timewarn; disp(['end frame of movie #' num2str(i) ' appears to be ' num2str(viddifs(i)) ' s off from video''s time stamp']); end
%     if simpleread && numrepeats(i)~=0; disp(['movie #' num2str(i) ' appears to have ' num2str(numrepeats(i)) ' repeated frames']); end
% end; xlabel('Frame #'); 
% set(gca,'yticklabel',datestr(get(gca,'ytick'),'HH:MM:SS.fff'));  
% % these lines examine any repeated time stamps and offer a solution if
% % appropriate
% % i = 25;  oi2 = find(abs(diff(frameTimes{i})) < median(diff(frameTimes{i}))/2 | (diff(frameTimes{i}) < 0 & diff(frameTimes{i}) > -median(diff(frameTimes{i}))));
% % for ii = 1:length(oi2); disp(frameTimes{i}(oi2(ii)-3:oi2(ii)+3)); end
% % for ii = 1:length(oi2); frameTimes{i}(oi2(ii)) = mean(frameTimes{i}([oi2(ii)-1 oi2(ii)+1])); end
% % for ii = 1:length(oi2); disp(frameTimes{i}(oi2(ii)-3:oi2(ii)+3)); end
% if any(abs(viddifs)>timewarn)
%     warning('Some movies had initital times (vidDN) + durations that were offset from the final frame embedded in the movie by the above listed amounts.  Recommend adjusting these times if they are small (< 1 s) and consistent throughout the movie times read process (see notes above for the specific video).  If they are big, recommend rerunning cell 1 with simpleread set to off for the movies that need it (this reads the embedded timestamp on each frame)');
%     
%     pp = input('adjust vidDN by the above listed amounts? 1 = yes, 2 = no ');
%     if pp == 1
%         for i = 1:length(frameTimes); 
%             if isempty(frameTimes{i}); continue; end
%             if abs(viddifs(i))>timewarn
%                 vidDN(i) = vidDN(i) + viddifs(i)/24/60/60;
%                 disp(['Movie ' num2str(i) ' start time adjusted by ' num2str(viddifs(i)) 's ']);
%             end
%         end
%     end
% end

% frameSize = [vid.width vid.height];
frameSize = [OBJ.Width OBJ.Height];
else
    frameSize = [nan nan];
end
vidDurs(movN(mlast-mfirst+1)+1:end) = []; frameTimes(movN( mlast-mfirst+1)+1:end) = []; vidDN(movN( mlast-mfirst+1)+1:end) = []; vidNam(movN( mlast-mfirst+1)+1:end) = []; try oframeTimes(movN( mlast-mfirst+1)+1:end) = []; catch; end
vid4k = true;
save([dataloc datafile(1:end-4) 'movieTimes.mat'],'vidDurs','frameTimes','movies','vidDN','vidNam','frameSize','vid4k');
if timestamps; save([dataloc datafile(1:end-4) 'movieTimes.mat'],'oframeTimes','-append'); end
disp('movieTimes file saved successfully, truncated to only include on whale files.');

warning('on','all');