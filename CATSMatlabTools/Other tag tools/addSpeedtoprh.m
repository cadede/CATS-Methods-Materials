% First run importd3tagdata to get a mat file with Adata, and also ensure you have the wav files available
% 

% make flownoise and jiggle variables
cf = pwd; %cd([vol ':\Tag data\Dtags\SOCAL']);
[prhfile,prhloc]=uigetfile('*.mat', 'select prh file');
cd(prhloc)
[filename,fileloc]=uigetfile('*.mat', 'select mat file with Adata');
[~,audiodir]=uigetfile('*.wav', 'select a wav file in folder with wav data');
cd(cf);

vars = load([prhloc prhfile],'vidDN','timedif','audstart','vidDurs','vidNum','fs','ofs','camon','tagon','nocam','noaud','nopress','df','CAL','Hzs','DN','whaleName'...
    ,'p','pitch','A','At','roll');
names = fieldnames(vars);
for i = 1:length(names)
    eval([names{i} ' = vars.' names{i} ';']);
end
vars.Depth = p;
load([fileloc filename])
%
if ~exist('vidDN','var') || isempty(vidDN) 
    vidDN = []; vars.vidDN = vidDN; 
  warning('vidDN variable is empty (no timestamps from video/audio files), will assume audio files start at start of data files'); 
end
if ~exist('DN','var'); DN = (0:length(vars.p)-1)'; DN = DN/fs/24/60/60 + Atime(1); vars.DN = DN; end
if ~exist('ODN','var'); ODN = Atime(1); vars.ODN = ODN; end
if ~exist('nocam','var'); nocam = true; vars.nocam = nocam; end
if ~exist('audstart','var'); audstart = ODN; warning('no audstart variable found, assuming audio starts at original data start time'); end
if ~exist('At','var'); At = A; vars.At = At; end
if ~exist('tagon','var'); tagon = gettagon(p,fs,DN(1),At); vars.tagon = tagon; end
if ~exist('timedif','var'); timedif = 0; end
vars.audstart = audstart;
if ~exist('noaud','var') || ~noaud
[flownoise,AUD] = getflownoise(audiodir,vars); noaud = false;
else; noaud = true;
end

tag1 = find(vars.tagon,1);
tag2 = find(vars.tagon,1,'last');
    disp('Done importing, check out figure 300 to examine data for outliers');
% plot data.  Look for outliers, may have to remove data above a threshold
% if there are spikes (sometimes happens at start and end of recordings)
if sum(isnan(flownoise)) ~= length(flownoise)
    figure(300); clf; set(300,'windowstyle','docked');
    ax = plotyy(DN(tag1:tag2),p(tag1:tag2),DN(tag1:tag2),flownoise(tag1:tag2));
    set(ax(1),'ydir','rev');
    legend('Depth','Flow noise (dB)');
    ylabel('Depth (m)','parent',ax(1));
    ylabel('Flow noise (dB)','parent',ax(2));
    %         text(min(get(gca,'xlim')),max(get(gca,'ylim')),'Press Enter if okay, or click on the threshold above which points are considered outliers','verticalalignment','top','fontsize',16,'parent',ax300);
    %         [~,y,button] = ginput(1);
    
    %         if ~isempty(button)
    %             flownoise(flownoise>y) = nan;
    %             plot(tag1:tag2,DB(tag1:tag2),'s');
    %         end
end
% clear vars
if ~exist('flownoise','var') || isempty(flownoise); warning('no audio files detected, flownoise is all nans'); flownoise = nan(size(Depth)); end
try Afs = Hzs.Afs; catch; Afs = round(1/mean(diff(Atime(50:60))*24*60*60)); end
disp(['Accelerometer sample rate is ' num2str(Afs) ' Hz'])

if Afs>180; maxfilt = 90; else maxfilt = round(.9*Afs/2); warning('Acc sample rate is less than 200 Hz, results of speed calibraiton may be unreliable'); end
if Afs>100; minfilt = 10; else minfilt = ceil(.1*Afs/2); end
try 
JX = TagJiggle(Adata(:,1),Afs,fs,[minfilt maxfilt],.5,Atime+timedif/24,DN); % 10 and 90 are the high-pass and low-pass filter frequencies. The higher number will have to be < .5* Afs.
JY = TagJiggle(Adata(:,2),Afs,fs,[minfilt maxfilt],.5,Atime+timedif/24,DN);
JZ = TagJiggle(Adata(:,3),Afs,fs,[minfilt maxfilt],.5,Atime+timedif/24,DN);
J = TagJiggle(Adata,Afs,fs,[minfilt maxfilt],.5,Atime+timedif/24,DN);
Jig = [JX JY JZ J];
catch
    warning('Error running TagJiggle, perhaps acc sample rate is lower than 180?  Can adjust high-pass filter in above lines to try again');
    Jig = nan(length(p),4);
    J = nan(size(p)); JX = J; JY = J; JZ = J;
end

% speedP = Paddles; speedP(speedP == 0) = nan;

% use this to examine the two metrics of turbulent flow.  They should
% align, else you may have an offset issue between the data and the 
% acoustics (and likely video)
JJ = J; JJ(isnan(JJ)) = 0; JJ = runmean(JJ,fs);
D = flownoise; D(isnan(D)|isinf(D)) = min(D(~isinf(D)));  D = runmean(D,fs);
figure; if sum(isnan(D))~=length(D); plotyy(DN,JJ,DN,D); else plot(DN,JJ,DN,D); end
legend('JiggleRMS','FlownoiseRMS')
JigRMS = Jig;
save([prhloc prhfile],'flownoise','At','DN','tagon','JigRMS','-append'); %'ODN','audstart','audon','camon','')
%


%%
% set threshold parameters
minDepth = 5;
minPitch = 45;
minSpeed = 1;

slips = [tag1 tag1; tag2 tag2];
speedper = [tag1 tag2];

% speedEnds([1 4 5 end-1:end]) = [];
% speedper = [1 430000; 430000 speedper(end)];

% if no flownoise, check if there is a speed sensor else only use one metric
if sum(isnan(flownoise)) == length(flownoise)
    RMS2 = []; lab = '';% could set RMS2 = Jig(:,4); lab = 'magJ'; if you want to compare the multiaxes model jig to the overall magnitude model
    try paddles = data.Speed; 
        RMS2 = decimateM(paddles,ofs,Hzs.SHz,df,length(JJ),'paddles',true); 
        RMS2 = runmean(RMS2,fs);
        lab = 'PW';
        Paddles = RMS2;
    catch
    end
else
    RMS2 = flownoise; lab = 'FN';
end
if ~exist([fileloc 'SpeedPlots//'],'dir'); mkdir([fileloc 'SpeedPlots//']); end
binSize = 1/fs*5;
if exist('Jig','var') && sum(isnan(Jig(:,1)))~=length(Jig(:,1))
    [~,speed,speedstats] = SpeedFromRMS3(Jig(:,1:3),'JJ',RMS2,lab,fs,p,pitch,roll,DN,speedper,slips,tagon,binSize,0.5,minDepth,minPitch,minSpeed,.2);
    X = Jig(:,1); Y = Jig(:,2); Z = Jig(:,3); Mag = Jig(:,4);
    for fig = [1 301:300+size(speedstats.r2used,1)]
        saveas(fig,[fileloc 'SpeedPlots//fig' num2str(fig) '.bmp']);
    end
else
    JJ = nan(size(p)); speed=table(JJ);
    X = nan(size(p)); Y = X; Z = X; Mag = X; 
end
    JigRMS = table(X, Y, Z, Mag);

%

if ~isempty(RMS2) %sum(isnan(flownoise)) ~= length(flownoise)
    s = input('Would you like to recalibrate speed from 2nd variable using its own sections (1 = yes, 2 = no- click no if current calibration is good)? ');
    if s == 1
        disp('Can quit out of this and start cell again later if the results don''t seem to be improving');
        [~,speedFN,speedstatsFN] = SpeedFromRMS3(RMS2,lab,[],'',fs,Depth,pitch,roll,DN,speedper,slips,tagon,.5,0.5,minDepth,minPitch,minSpeed,.2);
        if sum(isnan(JigRMS.X)) == length(JigRMS.X)
            speedstats = speedstatsFN;
        end
        oi = speedFN.Properties.VariableNames;
        oi(cellfun(@(x) strcmp('section',x), oi)) = {[lab 'FNsection']};
        oi(cellfun(@(x) strcmp('sectionUsed',x), oi)) = {[lab 'sectionUsed']};
        speedFN.Properties.VariableNames = oi;
        for i = 1:length(oi); speed.(oi{i}) = speedFN.(oi{i}); end
        speedstats.FN.Models = speedstatsFN.Models;
        speedstats.FN.ModelFits = speedstatsFN.ModelFits;
        speedstats.FN.Thresh = speedstatsFN.Thresh;
        speedstats.FN.r2used = speedstatsFN.r2used;
        speedstats.FN.sections_end_index = speedstatsFN.sections_end_index;
        for fig = [1 301:300+size(speedstats.r2used,1)]
            saveas(fig,[fileloc 'SpeedPlots//' lab 'fig' num2str(fig) '.bmp']);
        end
    end
else
    speed.FN = nan(size(speed.JJ));
    try speed.FNP68 = []; speed.FNP95 = []; speed.FN95 = []; speed.FNr2 = [];catch; end
end


if ~exist('speedstats','var'); speedstats = struct(); end
save([prhloc prhfile],'speed','speedstats','-append'); %'
