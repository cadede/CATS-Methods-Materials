%
% David Cade
% version 9.29.15
% Goldbogen Lab
% Stanford University

%64-bit full version (with graph)
% rearranges tag videos into side x side format.  Cuts the video into
% chunks of length "dur" (in seconds).  Makes a new folder to put the cut pieces in.
% includes audio and video with the same frame rate as the original.  File
% sizes are much larger than the original.  If you want to play around with different
% compressors that may give lower file sizes or allow for the whole video
% to be written as one file, type "list = mmwrite('','ListAviVideoEncoders')" at a >> prompt to see your
% compressor options.  Then change "comp" to a number corresponding to the
% compressor you want to try.  Or download "freemake" video software and
% just splice together the resulting clips
tic;

firststart =000; % enter starttime of the first video you want to render (in seconds).  0 if the whole video, something else if you want to cut some off
lastend = 30000;% 78*60*60 + 3*60 + 43 - 100*60*60;%40; %enter a big number if you want to do the whole last video (anything longer than the length of the file- generally 100 hours in seconds to be bigger than the weird reading wireless videos)
swapRL = false; % set to true if you want R and L cameras to be reversed (good for fro-back deployments on the right side of the animal).

justflownoise = true; % set to true if you want to display one speed metric instead of the composite speed
justjiggle = true; % set to true if you want to use jiggle instead of flow noise
combinesmall = true; % set to 0 if you do NOT want files to be combined into the previous video.  Default settings are 30 minute videos as long as gaps are less than 3 minutes
smallfont = false; %processors seem to handle font size differently.  For some newer videos need to set this to true to keep the font reasonable in the data box
dur = 5; % how much of the video (s) to read each iteration.  inf = make the whole video.  if you're having memory problems, I feel bad for you, son.  Try lowering the duration
comp = 5; % 0 means, use the last listed video compressor.  Change this number to use a different one.  run the "list" line later to see all compresors, MJPEG seems to work well
% boxP = 300/2560;  % proportion of the bottom graph you want to be the data box
boxsize = 0; %boxP*2560; %size of data box in pixels.  300 works for current font size and height
filtspeed = true; % if you want to smooth the speed
% withgraph = 1;  % 0 = no graph, 1 = platypus data, 2 = dtag data
% withprh = 1; % adds pitch, roll, heading, but you need the file
% threegraphs = 1; % when you have prh, also includes the whole dive profile.

Hz = 10; %refresh rate of graph
gap = 2; %pixels between video and graph
% calibrate = 0; %if you want to automatically calibrate the timcal (only works if you have the first video file recorded)
lowrat = 1; % the ratio of the size of the lower graph to the upper one
% leftres = 1280; % if your left side monitor is lower res than your main monitor (assuming your main monitor isn't ~=2560).  Put nan if it's the same.
% makenew = true;

% a = getdrives;
% for i = 1:length(a)
%     [~,vol]=system(['vol ' a{i}(1) ':']);
%     if strfind(vol,'CADE2'); vol = a{i}(1); break; end
% end
co = [0 0 1;
      0 0.5 0;
      1 0 0;
      0 0.75 0.75;
      0.75 0 0.75;
      0.75 0.75 0;
      0.25 0.25 0.25];
set(groot,'defaultAxesColorOrder',co);

cf = pwd; %cd([vol ':\CATS']);
[filename,fileloc]=uigetfile('*.*', 'select video files from one deployment','multiselect','on'); % can do multiple files successively
if ischar(filename); filename = {filename}; end
cd(fileloc); [prhfile,prhloc] = uigetfile('*.mat','select prh file');

try cd('Z:\New stuff'); catch; end
[~,filedest] = uigetfile('*.*','choose any file in the directory you want to put the partial files in, press cancel to use the same folder as the videos');
cd(cf);
if sum(filedest==0) || isempty(filedest); filedest = fileloc; end
%
load([prhloc prhfile]); %viddeploy(1) = [];
starttime = [firststart zeros(1,length(filename)-1)]; %zeros(size(filename)); % can adjust these if you don't want to render the whole videos
et = [100*60*60*ones(1,length(filename)-1) lastend];
if INFO.tagnum>=40 && ~ismember(INFO.tagnum,[50 51])
    kitten = true; else kitten = false;
end


if kitten
    try D = dir(fileloc(1:end-4)); D = {D.name}'; load([fileloc(1:end-4) D{~cellfun(@isempty,cellfun(@(x) strfind(x,'movieTimes'),D,'uniformoutput',false))}],'frameTimes','oframeTimes','frameSize');
    catch; try D = dir(fileloc(1:end)); D = {D.name}'; load([fileloc(1:end) D{~cellfun(@isempty,cellfun(@(x) strfind(x,'movieTimes'),D,'uniformoutput',false))}],'frameTimes','oframeTimes','frameSize');
        catch; [framename,frameloc]=uigetfile('*.*', 'select movieTimes file for wireless videos'); load([frameloc framename],'frameTimes','oframeTimes','frameSize');
        end
    end
end

if ~exist('speed','var');
    speed = table(speedFN, cell(size(speedFN)),'VariableNames',{'comp','type'});  speed.type(:) = {'FN'};
end

if justflownoise
    if ~exist('speedFN','var'); speedFN = speed.FN; end
    if justjiggle; speedFN = speed.JJ; end
end

if filtspeed;
    try speed.comp = runmean(speed.comp,round(fs)); catch; end
    try speedFN(isnan(speedFN)) = 0; speedFN = runmean(speedFN,round(fs)); catch; end
end
%
% some info
Ag = nan(size(Aw));
if ~exist('pitchgy','var') || sum(isnan(pitchgy))== length(pitchgy); pitchgy = pitch; rollgy = roll; headgy = head; end
Ag(:,1) = sin(pitchgy); Ag(:,2) = -sin(rollgy).*cos(pitchgy);  Ag(:,3) = -cos(rollgy).*cos(pitchgy);
%static acceleration of the animal
Sa = Aw - Ag;
njerkSa = (9.81*fs)* sqrt(diff(Sa).^2*ones(3,1)) ; njerkSa(end+1) = njerkSa(end);
njerk = (9.81*fs)*sqrt(diff(Aw).^2*ones(3,1)) ; njerk(end+1) = njerk(end);

% size of single video
tic
[vid,aud] = mmread([fileloc filename{1}], [1 2]);
toc
if abs(vid.width/vid.height-16/9)>.1; warning ('Single video is not 16 x 9. '); end
vidW = vid.width; vidH = vid.height;
ysize = round(2.15/14.3*vidH*14.3/10);%260;%215; % vertical pixel size of graph % assumes a screen area of a media player that is 14.3 cm tall by 30.8 wide
xsize = round(21.5/17.7*vidW);
if vidW > 2500; adF = 14; elseif vidW > 1500; adF = 2; else adF = 0; end % adjust fontsize slightly
scrn = get(0,'screensize');
if scrn(3)>1500 && vidW<2500; adF = adF - 3; if vidW > 1500; adF = adF + 1; end; end

% basic top overall graph
DN2= 1:length(DN); %changes the DN graphs back to regular indices the lazy way (instead of deleting the DN wrapper)
clear pat;
fig = figure(2); clf; %set(fig2,'color',[.8 .8 .8]);
a = find(p>1,1,'first');
%     a = I29 - 400; % I set the I29 to be the index of the first video I wanted to see (i.e. I29 = find(T.DV(:,4)>=16&T.DV(:,5)>=24&T.DV(:,6)>=04,1,'first')), then it only graphs from there on out
b = find(p>1,1,'last');
nonpat = nan(size(viddeploy));
labs = num2str(viddeploy); %labs = labs(:,5:6);
vidDurs2 = vidDurs; vidDurs2(viddeploy(1)) = vidDurs2(viddeploy(1)) - starttime(1);
vidDN2 = vidDN; vidDN2(viddeploy(1)) = vidDN2(viddeploy(1))+starttime(1)/24/60/60;
if combinesmall; combos = getcombos(vidDN2,vidDurs,viddeploy,1800,1200); else combos = num2cell(viddeploy); end
for i = 1:length(combos); for j = 2:length(combos{i}); labs(viddeploy==combos{i}(j),:) = ' '; end; end
for i = 1:length(viddeploy)
    [~,a1] = min(abs(DN-(vidDN(viddeploy(i))+starttime(i)/24/60/60)));
    %         a1 = round(a1+starttime(i)*fs);
    plot(DN2([a1 a1]),[-10 1000],'r','linewidth',2); hold on;
    ys = get(gca,'ylim');
    b1 = a1 + min(round((et(i)-starttime(i))*fs),round((vidDurs(viddeploy(i))-starttime(i)) * fs)) -1;
    nonpat(i) = rectangle('position',[DN2(a1) ys(1) DN2(b1)-DN2(a1) ys(2)],'facecolor',[255 192 203]/255);
    oi = get(gca,'children');
    set(gca,'children',[oi(2:end); oi(1)]);
    text(DN2(a1),0,labs(i,:),'verticalalignment','bottom','fontsize',10+adF,'color','r')
end
plot(DN2(a:b),p(a:b),'linewidth',1.5);
set(fig,'units','pixels','Position',[1,50,xsize,round(ysize*lowrat)]); % assumes 2560 frame size, which it should be
set(gca,'xlim',[DN2(a) DN2(b)],'ylim',[0 max(p(a:b))],'fontsize',16 + adF);
oi = datestr(DN(get(gca,'xtick')),'HH:MM:SS');
set(gca,'xticklabel',oi,'ydir','rev');
xlabel('Local Time'); ylabel('Depth');
set(gca,'units','points');
pos = get(gca,'position'); pos(1) = 50;
set(gca,'position',pos); set(gca,'units','normalized'); pos = get(gca,'position');
set(gca,'position',[pos(1:2) 1-pos(1)-.005 pos(4)]);
xs = get(gca,'xlim');
ys = get(gca,'ylim');
text(xs(1),ys(2),'Full Deployment Record','fontsize',20 + adF,'verticalalignment','bottom');
prhN = regexp(prhfile,' ')-1;
saveas(fig,[prhloc prhfile(1:prhN) ' TDR' '.bmp']);

% make long prh graph
fig5 = figure(5); clf; set(fig5,'color','w');
set(gca,'units','points');
pos = get(gca,'position'); pos(1) = 50;
set(gca,'position',pos); set(gca,'units','normalized'); pos = get(gca,'position');
set(gca,'position',[pos(1:2) 1-2.3*pos(1) pos(4)]);
xs = xs(1):xs(end); %xs from above
[axp5, hp51,hp52] = plotyy(xs,pitch(xs)*180/pi,xs,head(xs)*180/pi);
hold(axp5(1),'on'); hold(axp5(2),'on');
hp53 = plot(xs,roll(xs)/pi*180,'r.','parent',axp5(2),'markersize',4);
hp54 = plot(-10,5,'r-','linewidth',3); hp55 = plot(-10,5,'b-','linewidth',3); hp56 = plot(-10,5,'k--','linewidth',2);
set(hp51,'color','g','linewidth',2); set(hp52,'color','b','marker','.','markersize',4,'linestyle','none');
set(axp5,'xlim',[xs(1) xs(end)]);
set(axp5(1),'ylim',[-90 90],'ytick',-90:30:90,'box','off','ycolor','g'); set(axp5(2),'ycolor','k','ylim',[-181 181],'ytick',-120:60:120);
ylabel('degrees (pitch)');
%         if sum([headgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))>60/180*pi; pitchgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))>30/180*pi]) > sum([headgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))<-60/180*pi; pitchgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))<-30/180*pi])
%             legloc = 'southwest'; else legloc = 'northwest'; end
%         leg = legend(hp6,'Gyros','location',legloc,'orientation','horizontal');
set(get(axp5(2),'ylabel'),'string',{'degrees' '(roll and head)  '},'fontsize',16 + adF);
set(get(axp5(1),'ylabel'),'fontsize',16 + adF);
set(fig5,'units','pixels','Position',[1,70,xsize,round(ysize*lowrat)]);%fix this
set(gca,'units','points'); % these four lines repeated to ensure it stretches out right
pos = get(gca,'position'); pos(1) = 50;
set(gca,'position',pos); set(gca,'units','normalized'); pos = get(gca,'position');
set(gca,'position',[pos(1) pos(2)+.03 1-2.5*pos(1)-.005 pos(4)]);
set(axp5,'xtick',round(xs(1)+(xs(end)-xs(1))/10:(xs(end)-xs(1))/10:xs(end)));
oi = datestr(DN(get(axp5(2),'xtick')),'HH:MM:SS');
set(axp5,'xticklabel',oi,'fontsize',16 + adF);

%     plot(axp5(1),xs,skipgaps(pitchgy(xs)*180/pi,100),'g--','linewidth',2);
%     plot(axp5(2),xs,skipgaps(rollgy(xs)*180/pi,100),'r--','linewidth',2);
%     plot(axp5(2),xs,skipgaps(headgy(xs)*180/pi,100),'b--','linewidth',2);
ys = get(axp5(1),'ylim');
xs = get(axp5(1),'xlim');
ltext = text(xs(1)-(xs(2)-xs(1))*.028,ys(1)-(ys(2)-ys(1))/40,'Local Time: ','parent',axp5(1),'verticalalignment','top','fontname',get(axp5(1),'fontname'),'fontsize',get(axp5(1),'fontsize'),'horizontalalignment','left');
saveas(fig5,[prhloc prhfile(1:prhN) ' prh' '.bmp']);
CONTINUE = false; % Should be false, but set to true if you want to continue a previously interrupted cycle.  FYI, best time to interrupt is when the graphs are blazing.
startref = starttime; % for adjusting starttime for wireless videos

% make dolphin figure
 boxsize = xsize-vidW -gap;
 boxH = round(457/300*boxsize); %457/30
 fig6 = figure(6); clf;
 dH = vidH - boxH - gap;
 dW = boxsize;
 set(fig6,'units','pixels','Position',[scrn(3)-2*dW,0,2*dW,min(2*dH,0.9*scrn(4))],'color','w');
 set(gca,'position',[.01 .01 .962 .98],'xtick',[],'ytick',[]);
 F = plot_3d_model;
 [fpk,q] = dsf(Aw(tagon,:),fs,fs); % determine dominant stroke frequency;
 disp(['dominant stroke frequency: ' num2str(fpk) ' quality: ' num2str(q)]);
 [bodypitch,bodyroll] = a2pr([Aw(:,1:2) -Aw(:,3)],fs,fpk/5); bodyroll = -bodyroll; %uses Johnson method and then rotates back to normal axis orientation.
 bodyhead = m2h([Mw(:,1:2) -Mw(:,3)],[Aw(:,1:2) -Aw(:,3)],fs,fpk/5);
 prh = [bodypitch -bodyroll bodyhead]; %negative roll 'cause that's how tagtools rolls.
 rot_3d_model_NED(F,prh(find(tagon,1),:));
 prh(~tagon,:) = nan;
 
 
%% 
tagnum = gettagnum(prhfile);
% if tagnum>4 && tagnum < 40; T = nan(size(T)); Light = nan(size(Light)); end
if CONTINUE; startn = n; j = max(1,j-1); else startn = 1; end % lastVideo = length of last video % lastVideo = 0;

for n = startn:length(filename) %24:26 21:23 18:20]
    if ~CONTINUE || n>startn || j == 1
        vidN = viddeploy(strcmp(filename{n},vidNam(viddeploy)));
        if isempty(vidN); continue; end
        %     pr = p;
        endtime = et(n);
        if kitten; endtime = oframeTimes{vidN}(find(frameTimes{vidN}<=et(n),1,'last')); end
        totaltime = endtime - starttime(n);
%         stopj = floor((totaltime-1-.001)/dur)+1; %-1 is a correction factor for videos that go just one second over the threshold, 'cause then what's the point?
        stopj = floor((totaltime-.001)/dur)+1; %got rid of correction factor from above, we'll include that last second thank you very much.
        if kitten && starttime(n)>0; stopj = floor(((endtime-oframeTimes{vidN}(find(frameTimes{vidN}>=starttime(n),1,'first')))-.001)/dur)+1; end
        startj = 1;
    else startj = j;
    end
    clear audAdj;
    try
        audAdj = load([prhloc 'VideoData\' filename{n}(1:end-4) 'audio.mat']);
        audAdj = audAdj.aud;
    catch
    end
    %     if kitten; lastframe = 0; end
    c = find(cellfun(@(x) ismember(vidN,x), combos));
    ci = find(combos{c}==vidN);
    if ci == 1; combos{c}(ismember(combos{c}, find(isnan(vidDN)))) = []; end
%     if isnan(vidDN(vidN)); disp(['Skipped video num ' num2str(vidN)]); continue; % if the video was so short you couldn't even place it, just don't do it.
%     else combos{c}(ismember(combos{c}, find(isnan(vidDN)))) = [];
%     end
    j = startj;
    if j == 1; endf = 0; end
    while j<=stopj
        %
        oj = j;
        if ci~=1; oj = j+1; end
        if j>stopj  %j*dur>endtime || dur*(j-1)+starttime(n) > endtime;
            break; end
        clear vid;
        if isinf(dur); [vid,aud] = mmread([fileloc filename{n}]); else
            ST0 = dur*(j-1)+starttime(n); if kitten; ST0 =  oframeTimes{vidN}(find(frameTimes{vidN}>=starttime(n),1,'first')) + dur*(j-1); end
            if j>1; ST0 = ST0-kitten*0.1; end
%             if j>1; error('fix me'); end
            [vid,aud] = mmread([fileloc filename{n}], [],[ST0 min(dur+ST0,endtime)]); end
        if endtime>vid.totalDuration; endtime = vid.totalDuration; totaltime = endtime - starttime(n);
            if stopj>ceil((totaltime-1)/dur); stopj = ceil((totaltime-1)/dur); end; end %has to be before the kitten section since you are reading video from the raw (wrong) time
        startref(n) = starttime(n);
        if kitten
            readingweird = false;
            try TD = min(vid.totalDuration, endtime); catch; TD = 100*60*60; end
            if length(vid.frames)<90
                disp(['Video ' filename{n} ' is reading weirdly (or the last partial is very short), checking all frames for information']);
                readingweird = true;
%                 kk = 0;
                ST =ST0; %dur*(j-1)+starttime(n);
                ET = min([dur*j+ST0,endtime, 100*60*60]); % endtime
                while length(vid.frames)<150 && ET+.01<TD && ET < 100*60*60; % got rid of 30+60 since ET will always be shorter
                    clear vid;
                    j = j+1;
                    ET = min([dur*j+ST0,endtime, 100*60*60]); % endtime
                    [vid,aud] = mmread([fileloc filename{n}], [],[ST ET]);
                    if mod(round((ET-ST)*2/60)/2,300) == 0; kk = round((ET-ST)/60/60); disp([num2str(kk) ' "hours" read']);
                    end
                end
            end
            sf = find(oframeTimes{vidN} == vid.times(1));
            df = endf-sf+1;
            if sf<endf; vid.times = vid.times (df+1:end); vid.frames = vid.frames(df+1:end); sf = sf+df; end % shorten the videos by as much earlier you started the vdieo (accounts for the propensity of these videos to take a random frame as the first frame, messing up everything
            endf = min(find(oframeTimes{vidN} == vid.times(end)));
            if isempty(sf); error('vid.times(1) not fround in oframeTimes'); end
            if isempty(endf); %% allows for up to 0.1 seconds of extra frames in the movie
                [ed,endf] = min(abs(oframeTimes{vidN}-vid.times(end))); 
                if ed>0.1;  error('vid.times(end) not fround in oframeTimes'); else disp(['End Frame was ' num2str(ed) 'seconds off, but proceeding']); end
%                 
%                 if j>=stopj; 
%                     
%                 else
%                     error('vid.times(end) not fround in oframeTimes');
%                 end
            end
            lf = sf+length(vid.times)-1; ef = 0;
            if lf>length(frameTimes{vidN}); lf = length(frameTimes{vidN}); ef = sf+length(vid.times)-1 - lf; disp(['Warning: video ' filename{n} ' length had ' num2str(ef) ' extra frame(s) read, adding time for those frames']); end
            vid.times(1:end-ef) = frameTimes{vidN}(sf:lf); if ef~=0; vid.times(end-ef+1:end) = vid.times(end-ef)+1/30:1/30:vid.times(end-ef)+ef/30; end
            %             if vid.totalDuration>60*60;
            vid.totalDuration = frameTimes{vidN}(end)+1/30; %end
            if starttime(n) ~= 0 && j == 1; startref(n) = vid.times(1); end
            
        end
        j = j+1;
        
        %         if vid.width == 1280; vidW = 2560; else vidW = vid.width; end
        %         boxsize = vidW*boxP;
        
        % make figure
        fig1 = figure(1); clf; set(fig1,'color','w')
        vidstampend = vidDN(combos{c}(end))+min([vidDurs(combos{c}(end)) et(n+length(combos{c})-ci)])/24/60/60;
        %         sp1 = subplot(4,1,1:3);
        vidstamp1 = vidDN(combos{c}(1));
        vidstamp1 = vidstamp1 + startref(viddeploy == combos{c}(1))/24/60/60;
        vidstamp = vidDN(vidN);
        vidstamp = vidstamp + startref(n)/24/60/60;
        dataendtime = (vidstampend -vidstamp1)*24*60*60;
        [~,b] = min(abs(DN-vidstamp1));
        if justflownoise
            [ax1, h1, h2]= plotyy(DN2(b:round(fs*(dataendtime)+b)),p(b:round(fs*(dataendtime)+b)),DN2(b:round(fs*(dataendtime)+b)),speedFN(b:round(fs*(dataendtime)+b),:)); hold on;
        else
            [ax1, h1, h2]= plotyy(DN2(b:round(fs*(dataendtime)+b)),p(b:round(fs*(dataendtime)+b)),DN2(b:round(fs*(dataendtime)+b)),speed(b:round(fs*(dataendtime)+b),:),@plot,@plotspeed); hold on;
        end
        set([h1; h2],'linewidth',2);
        if justflownoise
            ymax = max(speedFN(b:round(fs*(dataendtime)+b)));
        else
            ymax = max(speed.comp(b:round(fs*(dataendtime)+b)));
        end
        if ymax == 0 || isnan(ymax); ymax = 1; end
        set(ax1(2),'ylim',[0 ymax],'nextplot','add');
        set(fig1,'units','pixels','Position',[1,400,xsize,ysize]);
        %             set(fig,'units','normalized');
        %             p = get(fig,'position');
        %             set(fig,'units','normalized','Position',[0,.1,1,p(4)]);
        set(ax1,'xlim',[DN2(b) DN2(ceil(fs*(dataendtime)+b))]); % changed round to ceil for a unique scenario when you need one more hit
        ylim([-2 max(p(b:round(fs*(dataendtime)+b)))]);
        %         set(ax1,'xtick',b+fs*60:fs*60:b+20*fs*60)
        set(ax1(1),'ydir','rev','box','off');
        tagnan = double(tagon); tagnan(~tagnan) = nan;
        oi = max(njerk(b:round(fs*(dataendtime)+b)).*tagnan(b:round(fs*(dataendtime)+b)));
        oi = ceil(oi/ymax*2)/2;
        h3 = plot(ax1(2),b:round(fs*(dataendtime)+b),njerk(b:round(fs*(dataendtime)+b)).*tagnan(b:round(fs*(dataendtime)+b))/oi,'m');
        try usePaddles = INFO.usePaddles; catch; end;
        if usePaddles; speedoi = 'SpeedP'; else speedoi = 'Speed'; end
        set(get(ax1(2),'ylabel'),'string',[speedoi ' (m/s)'],'fontsize',16 + adF);
        %         if vidW == 2560; JP = 1.037; else JP = 1.062; end
        if vidW>2500; JP = 1.039; SPACE = '      '; elseif vidW > 1500; JP = 1.044; SPACE = '    '; else JP = 1.057; SPACE = ''; end; if ~smallfont; JP = JP - 0.01; end
        jtext = text(DN2(round(fs*(dataendtime)*JP+b)),max(get(gca,'ylim')),[SPACE 'Jerk/' num2str(oi) ' (m/s^{3})'],'rotation',90,'fontsize',16 + adF,'color','m');
        %             xlabel('Local Time');
        ylabel('Depth'); set(get(gca,'ylabel'),'fontsize',16 + adF);
        set(gca,'units','points');
        pos = get(gca,'position'); pos(1) = 50;
        set(gca,'position',pos); set(gca,'units','normalized'); pos = get(gca,'position');
        if vidW>2500
             pos = [0.024 0.1539 1-4.2*0.0143 .97-0.1539];
        elseif vidW>1500
            if scrn(3)<1500; pos = [pos(1) .95*pos(2) 1-2.5*pos(1) .97-.95*pos(2)];
            else pos = [pos(1) 0.1539 1-2.5*pos(1) .97-0.1539]; end
        elseif scrn(3)<1500
            pos = [pos(1) 1.55*pos(2) 1-2.5*pos(1) .97-1.55*pos(2)];
        else
            pos = [pos(1) .1974 1-2.5*pos(1) .97-.1974];
        end
        
        set([gca ax1],'position',pos,'xticklabel',[]);
        
        %         sp2 = subplot(4,1,4);
        %         h4 = plot(sp2,DN2(b:round(fs*(endtime-starttime(n))+b)),Sa(b:round(fs*(endtime-starttime(n))+b),1),'r');
        MSA = msa(Aw);
        xs = b:round(fs*(dataendtime)+b);
        %         set(sp2,'position',[pos(1) .11 1-2.5*pos(1) pos(2)-.115]);
        %         set(sp2,'xlim',[DN2(b) DN2(round(fs*(endtime-starttime(n))+b))]);
        set(ax1,'xtick',round(xs(1)+(xs(end)-xs(1))/10:(xs(end)-xs(1))/10:xs(end)));
        
        %         if oj == 1 && max(DN(get(ax1(1),'xlim'))-DN(b))*24*60*60>combinesmall %new
        %             lastVideo = 0;
        %         end
        
        xlab = datestr(DN(get(ax1(1),'xtick'))-DN(b),'MM:SS'); %datestr(DN(get(gca,'xtick')),'HH:MM:SS'); %+lastVideo/24/60/60
        set(ax1,'xticklabel',xlab,'fontsize',16 + adF);
        %         set(sp2,'xticklabel',xlab,'fontsize',16,'ylim',[min(Sa(xs,1)) max(Sa(xs,1))],'ycolor','r');
        %         ylabel('Sa_{x}','parent',sp2);  pos = get(get(sp2,'ylabel'),'position');
        xs = get(ax1(1),'xlim');
        %         set(get(sp2,'ylabel'),'position',[xs(1)-(xs(2)-xs(1))*.017 pos(2) pos(3)]);
        ys = get(ax1(1),'ylim');
        try if ~any(get(ax1(1),'ytick')>0 & get(ax1(1),'ytick')<max(get(ax1(1),'ylim'))); set(ax1(1),'ytick',[0 max(floor(max(ys)), 10*floor(max(ys-1)/10))]); end
        catch
            if min(p(xs(1):xs(2)))+1>max(p(xs(1):xs(2))); set(ax1(1),'ytick',-1:1:1);
            else
                error('Check Depth');
            end
        end
        vtext = text(xs(1)-(xs(2)-xs(1))*.028,ys(2),'Time: ','parent',ax1(1),'verticalalignment','top','fontname',get(ax1(1),'fontname'),'fontsize',get(ax1(1),'fontsize'),'horizontalalignment','left');
        %
        if vidN < 10; vN = ['0' num2str(vidN)]; else vN = num2str(vidN); end
        vNs = num2str(combos{c}');
        vidts = nan(size(vNs,1),1);
        for i = 1:size(vNs,1)
            [~,a1] = min(abs(DN-(vidDN(combos{c}(i))+startref(viddeploy == combos{c}(i))/24/60/60)));
            %         a1 = round(a1+starttime(i)*fs);
            if i ==1; a1 = a1 + 10; end
            plot(DN2([a1 a1]),[-10 1000],'r','linewidth',2); hold on;
            ys = get(gca,'ylim');
            b1 = a1 + min(round((et(viddeploy == combos{c}(i))-startref(viddeploy == combos{c}(i)))*fs),round((vidDurs(combos{c}(i))-startref(viddeploy == combos{c}(i))) * fs)) -1;
            if i == ci; fc = [255 250 205]; else fc = [255 240 245]; end
            rectangle('position',[DN2(a1) -10 DN2(b1)-DN2(a1) 500],'facecolor',fc/255);
            oi = get(gca,'children');
            set(gca,'children',[oi(2:end); oi(1)]);
            if i == 1; y2 = mean(ys); else y2 = 0; end
            vidts(i) = text(DN2(a1+20),y2,vNs(i,:),'verticalalignment','top','horizontalalignment','left','fontsize',10,'color','r');
        end
        if ~exist ([prhloc 'graphs'],'dir'); mkdir ([prhloc 'graphs']); end
        if oj == 1; savefig(fig1,[prhloc 'graphs\' prhfile(1:prhN) ' (' vN ') speedgraph' num2str(fs) 'Hz.fig']); end
        %
        
        
        %
        % set up basic prh graph
        fig3 = figure(3); clf; set(fig3,'color','w');
        set(gca,'units','points');
        pos = get(gca,'position'); pos(1) = 50;
        set(gca,'position',pos); set(gca,'units','normalized'); pos = get(gca,'position');
        set(gca,'position',[pos(1:2) 1-2.3*pos(1) pos(4)]);
        %         vidstamp = vidDN(strcmp(filename{n},vidNam))-timcal;
        %         vidstamp = vidstamp + starttime(n)/24/60/60;
        %         [~,b] = min(abs(T.DN-vidstamp));
        %         c = round((b-I29+10.975*40)/4);
        xs = b:round(fs*(dataendtime)+b);
        [axp, hp1,hp2] = plotyy(DN2(xs),pitch(xs)*180/pi,DN2(xs),head(xs)*180/pi);
        hold(axp(1),'on'); hold(axp(2),'on');
        hp3 = plot(DN2(xs),roll(xs)/pi*180,'r.','parent',axp(2),'markersize',5);
        hp4 = plot(-10,5,'r-','linewidth',3); hp5 = plot(-10,5,'b-','linewidth',3); hp6 = plot(-10,5,'k--','linewidth',2);
        set(hp1,'color','g','linewidth',3); set(hp2,'color','b','marker','.','markersize',5,'linestyle','none');
        set(axp,'xlim',DN2([xs(1) xs(end)])); set(axp(1),'ylim',[-90 90],'ytick',-90:30:90,'box','off','ycolor','g'); set(axp(2),'ycolor','k','ylim',[-181 181],'ytick',-120:60:120);
        ylabel('degrees (pitch)');
%         if sum([headgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))>60/180*pi; pitchgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))>30/180*pi]) > sum([headgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))<-60/180*pi; pitchgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))<-30/180*pi])
%             legloc = 'southwest'; else legloc = 'northwest'; end
%         leg = legend(hp6,'Gyros','location',legloc,'orientation','horizontal');
        set(get(axp(2),'ylabel'),'string',{'degrees' '(roll and head)  '},'fontsize',16 + adF);
        set(get(axp(1),'ylabel'),'fontsize',16 + adF);
        set(fig3,'units','pixels','Position',[1,70,xsize,round(ysize*lowrat)]);%fix this
        set(gca,'units','points'); % these four lines repeated to ensure it stretches out right
        pos = get(gca,'position'); pos(1) = 50;
        set(gca,'position',pos); set(gca,'units','normalized'); pos = get(gca,'position');
        if vidW > 1500
            set(gca,'position',[pos(1) pos(2)+.025 1-2.5*pos(1)-.005 pos(4)]);
        else
            set(gca,'position',[pos(1) pos(2)+.045 1-2.5*pos(1)-.005 pos(4)]);
        end
        set(axp,'xtick',round(xs(1)+(xs(end)-xs(1))/10:(xs(end)-xs(1))/10:xs(end)));
        oi = datestr(DN(get(axp(2),'xtick')),'HH:MM:SS');
        set(axp,'xticklabel',oi,'fontsize',16 + adF);
        
%         plot(axp(1),DN2(xs),skipgaps(pitchgy(xs)*180/pi,100),'g--','linewidth',2);
%         plot(axp(2),DN2(xs),skipgaps(rollgy(xs)*180/pi,100),'r--','linewidth',2);
%         plot(axp(2),DN2(xs),skipgaps(headgy(xs)*180/pi,100),'b--','linewidth',2);
        ys = get(axp(1),'ylim');
        xs = get(axp(1),'xlim');
        ltext = text(xs(1)-(xs(2)-xs(1))*.028,ys(1)-(ys(2)-ys(1))/40,'Local Time: ','parent',axp(1),'verticalalignment','top','fontname',get(axp(1),'fontname'),'fontsize',get(axp(1),'fontsize'),'horizontalalignment','left');
        %
        if vidN < 10; vN = ['0' num2str(vidN)]; else vN = num2str(vidN); end
        if oj == 1; savefig(fig3,[prhloc 'graphs\' prhfile(1:prhN) ' (' vN ') prhgraph' num2str(fs) 'Hz.fig']); end
        
        %            +(j-1)*dur/24/60/60; %time stamp of the video fragment
        
        if oj == 1 % top overall graph
            %             try delete(pat); catch; end
            fig2 = figure(2); hold on;
            [~,a] = min(abs(DN-vidstamp));
            b = a + round(dataendtime*fs);%round(vid.totalDuration * fs) -1;
            %             pat = patch([T.DN(a) T.DN(b) T.DN(b) T.DN(a)],[1000 1000 -10 -10],'y','facealpha',.2);
            ys = get(gca,'ylim');
            if ~exist('pat','var')
                pat = rectangle('position',[DN2(a) ys(1) DN2(b)-DN2(a) ys(2)],'facecolor','y');
                oi = get(gca,'children');
                oldpat = find(ismember(oi,nonpat),1,'first');
                set(gca,'children',[oi(2:oldpat-1); oi(1); oi(oldpat:end)]);
            else
                set(pat,'position',[DN2(a) ys(1) DN2(b)-DN2(a) ys(2)]);
            end
            set(fig2,'units','pixels','Position',[1,50,xsize,ysize]);
            topM = getframe(fig2);
            scrn = get(0,'screensize');
            pos = get(gcf,'position');
            while pos(3)+pos(1)-1>scrn(3)
                pos(1) = pos(1) - scrn(3) ;
                set(fig2,'position',pos);
                topM2 = getframe(fig2);
                topM.cdata = [topM.cdata topM2.cdata];
            end
            %             if size(topM.cdata,2)< 2560
            %                 topM.cdata= [zeros(size(topM.cdata,1),floor((2560-size(topM.cdata,2))/2),size(topM.cdata,3)) topM.cdata zeros(size(topM.cdata,1),ceil((2560-size(topM.cdata,2))/2),size(topM.cdata,3))];
            %             end
        end
        %         if vidW == 2560; vidH = vid.height/2; else vidH = vid.height; end
        oi = uint8(zeros(vidH+2*ysize + 2*gap,xsize,3)); % a matrix for each recreated frame
        
        % info pane
        if smallfont; adj = 2; else adj = 0; end
        
        
        istamp = vidstamp+vid.times(1)/24/60/60-startref(n)/24/60/60;
        [~,b] = min(abs(DN-istamp));
        DNclose0 = round(DN(b)*24*60*60*Hz)/24/60/60/Hz;
        %
        boxsize = xsize-vidW -gap;
        boxH = round(457/300*boxsize); %457/300 is box ratio from sxs style
        if vidW > 3000
            adj = -28; %makes the font bigger for a bigger box
        elseif vidW > 1500
            adj = -6;
        else
            adj = 1;
        end
        if scrn(3)>1500; adj = adj + 2; if vidW > 1500; adj = adj +1; end; end
        for i = 1:length(vid.frames)
            if tagnum == 40 && DN(1)<datenum([2017 08 12 0 0 0])
                vid.frames(i).cdata = rot90(vid.frames(i).cdata,2);
            end
            %             if i == 299; error ('fix me'); end
            istamp = vidstamp+vid.times(i)/24/60/60-startref(n)/24/60/60;
            [~,b] = min(abs(DN-istamp));
            DNclose = round(DN(b)*24*60*60*Hz)/24/60/60/Hz;
            %             disp(i)
            %             disp(datestr(istamp,'HH:MM:SS.fff'));
            %             disp(datestr(DNclose,'HH:MM:SS.fff'));
            if i == 1 || DNclose>=DNclose0+.99/24/60/60/Hz% .99 is the wiggle room for rounding errors with doubles
                DNclose0 = DNclose;
                fig4 = figure(4); clf;
                set(fig4,'units','pixels','Position',[1,50,boxsize,boxH],'color','w');
                set(gca,'position',[.01 .01 .962 .98],'xtick',[],'ytick',[]);
                box on;
                [~,b2] = min(abs(DN-(istamp+1/24/60/60/Hz))); b2 = b2 -1;
                iD = datestr(DNclose,'HH:MM:SS.fff');
                %                 disp(iD);
                text(.05,.95,['Time = ' iD(1:end-2)],'color','k','fontsize',16-adj);
                text(.05,.88,['Depth = ' sprintf('%.1f', mean(p(b:b2))) ' m'],'color','b','fontsize',16-adj);
                if usePaddles; speedoi = mean(speedP(b:b2)); elseif justflownoise; speedoi = nanmean(speedFN(b:b2)); else speedoi = mean(speed.comp(b:b2)); end
                if justjiggle; spoi = 'JJ'; elseif justflownoise; spoi = 'FN'; elseif usePaddles; spoi = 'P'; else spoi = ''; end
                text(.05,.81,['Speed' spoi ' = ' sprintf('%.1f', speedoi) ' m/s'],'color','g','fontsize',16-adj,'fontweight','bold');
                text(.05,.74,['Jerk = ' sprintf('%.1f', mean(njerk(b:b2))) ' m/s^{3}'],'color','m','fontsize',16-adj);
                %                 text(.05,.67,['Pitch (Gyro)= ' sprintf('%.0f',circ_mean(pitch(b:b2))*180/pi) '{\circ} (' sprintf('%.0f',circ_mean(pitchgy(b:b2))*180/pi) '{\circ})'],'color','g','fontsize',16-adj,'fontweight','bold');
                %                 text(.05,.60,['Roll (Gyro)= ' sprintf('%.0f',circ_mean(roll(b:b2))*180/pi) '{\circ} (' sprintf('%.0f',circ_mean(rollgy(b:b2))*180/pi) '{\circ})'],'color','r','fontsize',16-adj);
                %                 text(.05,.53,['Head (Gyro)= ' sprintf('%.0f',circ_mean(head(b:b2))*180/pi) '{\circ} (' sprintf('%.0f',circ_mean(headgy(b:b2))*180/pi) '{\circ})'],'color','b','fontsize',16-adj)
                %                 text(.05,.53,['|Sa| = ' sprintf('%.1f', sqrt(sum(Sa(b,:).^2))) ' m/s^{2}'],'color','k','fontsize',16);
                text(.05,.67,['Pitch = ' sprintf('%.0f',circ_mean(pitch(b:b2))*180/pi) '{\circ}'],'color','g','fontsize',16-adj,'fontweight','bold');
                text(.05,.60,['Roll = ' sprintf('%.0f',circ_mean(roll(b:b2))*180/pi) '{\circ}'],'color','r','fontsize',16-adj);
                text(.05,.53,['Head = ' sprintf('%.0f',circ_mean(head(b:b2))*180/pi) '{\circ}'],'color','b','fontsize',16-adj)
                text(.05,.46,['|Aw| (MSA) =' sprintf('%.2f', mean(sqrt(sum(Aw(b:b2,:).^2,2)))) ' g (' sprintf('%.2f', mean(sqrt(sum(MSA(b:b2,:).^2,2)))) ')'],'color','k','fontsize',16-adj);
                %                 text(.05,.46,['|Aw| (|Sa|) =' sprintf('%.2f', mean(sqrt(sum(Aw(b:b2,:).^2,2)))) ' g (' sprintf('%.2f', mean(sqrt(sum(Sa(b:b2,:).^2,2)))) ')'],'color','k','fontsize',16-adj);
                text(.05,.38,['Aw_{x} (Sa_{x}) = ' sprintf('%.2f', mean(Aw(b:b2,1))) ' g (' sprintf('%.2f',mean(Sa(b:b2,1))) ')'],'color','k','fontsize',16-adj);
                text(.05,.30,['Aw_{y} (Sa_{y}) = ' sprintf('%.2f', mean(Aw(b:b2,2))) ' g (' sprintf('%.2f',mean(Sa(b:b2,2))) ')'],'color','k','fontsize',16-adj);
                text(.05,.22,['Aw_{z} (Sa_{z}) = ' sprintf('%.2f', mean(Aw(b:b2,3))) ' g (' sprintf('%.2f',mean(Sa(b:b2,3))) ')'],'color','k','fontsize',16-adj);
                %                 if mean(Light(b:b2))>3499.8; gtl = ' {\geq} '; else gtl = ' = '; end
                %                 text(.05,.15,['Light' gtl sprintf('%.0f', mean(Light(b:b2))) ' '],'color','k','fontsize',16-adj);
                text(.05,.15,['Light = ' sprintf('%.0f', mean(Light(b:b2))) ' '],'color','k','fontsize',16-adj);
                text(.05,.085,['Temperature = ' sprintf('%.1f', mean(T(b:b2))) '{\circ}'],'color','k','fontsize',16-adj);
                text(.05,.025,[num2str(fs) ' Hz datafile index: ' num2str(b) ':' num2str(b2) ],'fontsize',12-adj);
                Mbox = getframe(fig4);
                
                % dolphin figure
                rot_3d_model_NED(F,circ_mean(prh(b:b2,:)));
                                
                try delete(h);catch; end; % add lines
                figure (1);
                h = plot(ax1(1),[DN2(b) DN2(b)],[-5 1000],'r','linewidth',2);
                if ci>1 && i == 1; x2 = x1; end
                if (vid.times(i)-startref(n)>3 && vid.times(i)-startref(n)<10 && ci == 1) || (ci == 1 && vid.times(end)<=3 && i>length(vid.times)/2)
                    x1 = min(get(ax1(1),'xlim')); x2 = x1;
                    %                     xend = max(get(sp1,'xlim'));
                    if DN(end)<=vidstamp1+dataendtime/24/60/60; % this if checks if the prh file goes to the end of the video.  Else if just graphs to the end of the prh file.
                        xend = length(DN);
                    else
                        xend = find(DN - (vidstamp1+dataendtime/24/60/60)>0,1,'first'); % finds the index just barely bigger then the end of the video
                    end
                    try delete(fullt);  catch; end
                    if length(combos{c}) == 1; fullt = text(x1,0,['Full Video ' num2str(vidN) ' Record'],'verticalalignment','top','horizontalalignment','left','fontsize',20+adF,'fontweight','bold','parent',ax1(1));
                    else fullt = text(x1,0,['Full Video ' num2str(vidN) ' - ' num2str(combos{c}(end)) ' Record'],'verticalalignment','top','horizontalalignment','left','fontsize',20+adF,'fontweight','bold','parent',ax1(1));
                    end
                elseif vid.times(i)-startref(n)>3 || ci~=1
                    % add in here adjusting y values of speedFN and jerk
                    try delete(fullt);  catch; end
                    if (istamp>DN(round(x2)) || (i == 1 && ci>1))&&x2<xend-2*fs %add in some wiggle room, it's okay if less than 2 seconds go past the end time (added since the new videos were reading weirdly)
                        if x2>x1; [~,x1] = min(abs(DN-istamp)); try delete(leg); catch; end; end
                        x2 = min([x1 + 10*fs*60; xend-4*fs*60]); % if the end is less than 14 minutes away, just go there
                        if x2<x1 + 10*fs*60; x2 = xend; end
                    end
                    set(vidts,'visible','on');
                    tx = get(vidts,'Position'); if iscell(tx); tx = vertcat(tx{:}); end; tx = tx(:,1);
                    set(vidts(tx<x1|tx>x2),'Visible','off');
                    
                    set([ax1 axp],'xlim',[x1 x2]);
                    %                         [~,onminute] = min(abs((DN(x1+(58*fs:62*fs))-vidstamp1)*24*60-1)); %should account for rounding errors so that ticks are on the minute?
                    onminute = 0;
                    ticks = x1+60*fs+onminute:60*fs:x2;
                    if (x2-x1)/fs/60<6; ticks = x1+60*fs+onminute-30*fs:30*fs:x2; end
                    if (x2-x1)/fs/60<3; ticks = x1+60*fs+onminute-30*fs-10*fs:20*fs:x2; end
                    if isempty(ticks); ticks = round(mean(x1,x2)); end
                    set([ax1 axp],'xtick',ticks);
                    xlab = datestr(DN(get(gca,'xtick')),'HH:MM:SS');
                    %                         [~,startb] = min(abs(DN-vidstamp));
                    xlab2 = datestr(DN(get(gca,'xtick'))-vidstamp1,'MM:SS'); %+lastVideo/24/60/60,
                    set(axp,'xticklabel',xlab,'fontsize',16 + adF);
                    set(ax1,'xticklabel',xlab2,'fontsize',16 + adF);
                    %                         [~,b1] = min(abs(DN-DN(x1))); [~,b2] = min(abs(DN-DN(x2)));
                    if justflownoise; ymax = max(speedFN(x1:x2));
                    else ymax = max(speed.comp(x1:x2)); end
                    if ymax == 0 || isnan(ymax); ymax = 1; end
                    set(ax1(2),'ylim',[0 ymax]);
                    pos = get(jtext,'position');  pos(1) = (x2-x1)*JP+x1; set(jtext,'position',pos);
                    pos = get(vtext,'position'); pos(1) = x1-(x2-x1)*.028; set(vtext,'position',pos);
                    pos = get(ltext,'position'); pos(1) = x1-(x2-x1)*.028; set(ltext,'position',pos);
                    %                         set(sp2,'ylim',[min(Sa(x1:x2,1)) max(Sa(x1:x2,1))]); %untested order
                    %                         pos = get(get(sp2,'ylabel'),'position');
                    %                         set(get(sp2,'ylabel'),'position',[x1-(x2-x1)*.017 pos(2) pos(3)]);
                end
                fig = fig1;
                if vid.times(i)-vid.times(1) <= 3 && oj == 1 && ci == 1 && ~(ci == 1 && vid.times(end)<=3 && i>length(vid.times)/2)% first three seconds
                    fig = fig2;
                end
                set(fig,'units','pixels','Position',[1,50,xsize,ysize]);
                M = getframe(fig);
                scrn = get(0,'screensize');
                pos = get(gcf,'position');
                if size(M.cdata,2)< xsize
                    while pos(3)+pos(1)-1>scrn(3)
                        pos(1) = pos(1) - scrn(3) ;
                        set(fig,'position',pos);
                        M2 = getframe(fig);
                        M.cdata = [M.cdata M2.cdata];
                    end
                end
                if size(M.cdata,2)< xsize
                    M.cdata= [zeros(size(M.cdata,1),floor((xsize-size(M.cdata,2))/2),size(M.cdata,3)) M.cdata zeros(size(M.cdata,1),ceil((xsize-size(M.cdata,2))/2),size(M.cdata,3))];
                end
                figure (3);
                try delete(hp);catch; end;
                hp = plot([DN2(b) DN2(b)],[-100 100],'k','linewidth',3);
                %                 c = round((b-I29+10.975*40)/4);
                
                figb = fig3;
                if vid.times(i)-vid.times(1) <= 3 && oj == 1 && ci == 1  && ~(ci == 1 && vid.times(end)<=3 && i>length(vid.times)/2)% first three seconds
                    figb = fig5;
                end
                set(figb,'units','pixels','Position',[1,50,xsize,round(lowrat*(ysize))]);
                Mp = getframe(figb);
                scrn = get(0,'screensize');
                pos = get(gcf,'position');
                if size(Mp.cdata,2)< xsize
                    while pos(3)+pos(1)-1>scrn(3)
                        pos(1) = pos(1) - scrn(3) ;
                        set(figb,'position',pos);
                        Mp2 = getframe(figb);
                        Mp.cdata = [Mp.cdata Mp2.cdata];
                    end
                end
                
                if size(Mp.cdata,2)< xsize
                    Mp.cdata= [zeros(size(Mp.cdata,1),floor((xsize-size(Mp.cdata,2))/2),size(Mp.cdata,3)) Mp.cdata zeros(size(Mp.cdata,1),ceil((xsize-size(Mp.cdata,2))/2),size(Mp.cdata,3))];
                end
                % corner dolphin
                figure(6);
                Df = getframe(fig6);
                dsize = size(Df.cdata);
                if mod(dsize(1),2) == 1; Df.cdata(end+1,:,:) = 1; end
                Df.cdata = Df.cdata(ceil(dsize(1)/2-dH/2 + 1):floor(dsize(1)/2+dH/2),ceil(dsize(2)/2-dW/2 + 1):floor(dsize(2)/2+dW/2),:);
                dsize = size(Df.cdata); if dsize(1)~=dH; Df.cdata(end+1,:,:) = 1; end; if dsize(2)~=dW; Df.cdata(:,end+1,:) = 1; end
                
            end
            %             if vidW == 2560; vidH = vid.height/2; else vidH = vid.height; end
            %             oi = uint8(zeros(vidH+round((1+lowrat)*ysize) + 2*gap,vidW,3));
            oi = uint8(zeros(vidH+2*ysize + 2*gap,xsize,3)); %
            %             oi(1:vidH,1:vid.width,:) = vid.frames(i).cdata;
            oi(1:vidH,xsize-vid.width+1:xsize,:) = vid.frames(i).cdata;
            %             if vidW == 2560
            %             oi(1:vidH,vid.width+1:end,:) = vid.frames(i).cdata(vid.height/2+1:end,:,:);
            %             end
            %             if swapRL; OI = oi(1:vidH,1:size(oi,2)/2,:); oi(1:vidH,1:size(oi,2)/2,:) = oi(1:vidH,size(oi,2)/2+1:size(oi,2),:); oi(1:vidH,size(oi,2)/2+1:size(oi,2),:) = OI; end %swaps L and R frames
            oi(vidH+gap+1:vidH+gap+ysize,:,:) = M.cdata;
            oi(vidH+2*gap+1+ysize:end,:,:) = Mp.cdata;
            %             oi(vidH-boxH+1:vidH,xsize-boxsize+1:end,:) = Mbox.cdata;
            oi(vidH-boxH+1:vidH,1:boxsize,:) = Mbox.cdata;
            if size(oi,1)/2 ~= round(size(oi,1)/2); oi = [oi; zeros(1,size(oi,2),3)]; end % if the video height isn't even
            % add dolphin
            oi(1:dH,1:dW,:) = Df.cdata;            
            vid.frames(i).cdata = oi;
        end
        if ci == 1; % if you are not adding on to a last video, make a new directory
            if ~exist ([filedest 'partial\' prhfile(1:regexp(prhfile,' ')-1) '\' filename{n}(1:end-4)],'dir'); mkdir ([filedest 'partial\' prhfile(1:regexp(prhfile,' ')-1) '\' filename{n}(1:end-4)]); end
            dirN = n; %else dirN stays the same
        end
        vid.width = xsize;
        vid.height = vidH+2*ysize + 2*gap;
        if vid.height/2 ~= round(vid.height/2); vid.height = vid.height+1; end
        list = mmwrite('','ListAviVideoEncoders');
        if ~comp; conf.videoCompressor = char(list(end)); else conf.videoCompressor = char(list(comp)); end
        %         list = mmwrite('','ListAviAudioEncoders');
        %         conf.audioCompressor = char(list(end));
        %         try
        if exist('audAdj','var')
            audA = audAdj;
            audI = [round(vid.times(1)*audAdj.rate) + 1, min(round(vid.times(end)*audAdj.rate)+1,size(audA.data,1))];
            audA.data = audA.data(audI(1):audI(2),:);
            audI = find(audAdj.times>=vid.times(1) & audAdj.times<vid.times(end));
            audA.times = vid.times(1); %audA.times(audI); audA.frames = audA.frames(audI);
        else audA = aud;
        end
%         if n~=15 || j>331
        if isinf(dur); mmwrite([filedest 'partial\' filename{dirN}(1:end-4) '\' filename{n}(1:end-4) 'SxS.avi'],vid,aud,conf); else
            if oj ==1; tail = []; else tail = [num2str(dur*(j-2)+starttime(n)) 'sec']; end
            if vidN < 10; vN = ['0' num2str(vidN)]; else vN = num2str(vidN); end
            if combos{c}(end)<10;  vN2 = ['0' num2str(combos{c}(end))]; else vN2 = num2str(combos{c}(end)); end
            if ci == 1 && oj == 1 && length(combos{c})>1; vN = [vN '-' vN2]; es = ' '; else es = ''; end
            mmwrite([filedest 'partial\' prhfile(1:regexp(prhfile,' ')-1) '\' filename{dirN}(1:end-4) '\' es prhfile(1:regexp(prhfile,' ')-1) ' (' vN ')' tail '.avi'],vid,audA,conf);
        end
%         end
        %         catch err
        %             if isinf(dur); mmwrite([filedest 'partial\' filename{n}(1:end-4) '\' filename{n}(1:end-4) 'SxS.avi'],vid,aud,conf); else
        %                 if j ==1; tail = []; else tail = [num2str(dur*(j-1)+starttime(n)) 'sec']; end
        %                 if vidN < 10; vN = ['0' num2str(vidN)]; else vN = num2str(vidN); end
        %                 mmwrite([filedest 'partial\' prhfile(1:10) '\' filename{n}(1:end-4) '\' prhfile(1:10) ' (' vN ')' tail '.avi'],vid,aud);
        %             end
        %         end
        % these lines should have allowed the continual writing of a file to
        % put all the pieces into one video, but it always crashed matlab for me
        %         if j == 1; mmwrite([fileloc filename{n}(1:end-4) '\' filename{n}(1:end-4) 'SxS.avi'],vid,aud,conf,'Continue'); end
        %         if j == stopj; mmwrite([fileloc filename{n}(1:end-4) '\' filename{n}(1:end-4) 'SxS.avi'],vid,aud,conf,'Initialized');
        %         else mmwrite([fileloc filename{n}(1:end-4) '\' filename{n}(1:end-4) 'SxS.avi'],vid,aud,conf,'Continue','Initialized'); end
        disp([num2str((j-1)*dur) ' sec written']);
    end
    %     if n == length(filename); curVideo = min(vid.totalDuration,endtime)-starttime(n); else curVideo = vid.totalDuration-starttime(n); end
    %     if curVideo<combinesmall; lastVideo = lastVideo+curVideo; elseif combinesmall; lastVideo = curVideo; end
    disp ([filename{n} ' complete']);
    if ci == length(combos{c})
        continue
    else
        tail = [num2str(dur*(j-2)+starttime(n)) 'sec'];
        if ~exist('frameSize','var'); vid2 = mmread([fileloc filename{n}],[1 5]); frameSize = [vid2.width vid2.height]; end
        blankframe = vid.frames(end).cdata;
        blankframe(1:frameSize(2),end-frameSize(1)+1:end,:) = 0;
        for ii = 1:90; vid.frames(ii).cdata = blankframe; end
        vid.frames(91:end) = [];
        vid.times(1:90) = 1/30:1/30:3; vid.times(91:end) = [];
        vid.rate = 30;
        vid.totalDuration = 3;
        tail = [tail 'gap'];
        mmwrite([filedest 'partial\' prhfile(1:regexp(prhfile,' ')-1) '\' filename{dirN}(1:end-4) '\' prhfile(1:regexp(prhfile,' ')-1) ' (' vN ')' tail '.avi'],vid,conf);
    end
    
end
clear vid
aoi = find(~isnan(vidDN),1,'first');
oi = find(isnan(vidDN)); oi(oi<aoi) = [];
if ~isempty(oi)
    badmoviesloc = [filedest 'partial\' prhfile(1:regexp(prhfile,' ')-1) '\badmovies\'];
    D = dir([fileloc 'bad movies\']); D = {D(~vertcat(D.isdir)).name};
    if ~exist(badmoviesloc,'dir'); mkdir(badmoviesloc); end
    for i = 1:length(D)
        clear vid aud
        [vid,aud] = mmread([fileloc 'bad movies\' D{i}]);
        mmwrite([badmoviesloc D{i}(1:end-3) 'avi'],vid,aud,conf);
    end
end
    
