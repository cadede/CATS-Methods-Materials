function resampleprh(newfs,prhloc,INFOloc)
%based on old file "adjustcatssamplerate"
% newfs = new sample rate (must be divisible by original sample rate of
% data
%prhloc = file location of prh file
% INFOloc = file location of raw data (INFO file)

    try cd(prhloc); catch; end
        [prhfile, prhloc]  = uigetfile('*prh.mat','Choose prh file');

try cd(INFOloc); catch; end
        [INFOfile, INFOloc]  = uigetfile('*INFO.mat','Choose INFO file with prh making information');
load([INFOloc INFOfile]);%,'noc
VARS = load([INFOloc INFOfile]);%,
VARS = fieldnames(VARS);

if abs(Hzs.datafs/newfs - round(Hzs.datafs/newfs)) > .001; error(['new sample rate (' num2str(newfs) ') does not divide evenly into original sample rate (' num2str(Hzs.datafs) ')' ]); end

trunc = dir([prhloc '*truncate.mat']);
try rawfile = [trunc.name(1:end-12) 'truncate.mat'];
    load([prhloc rawfile] ,'data');
    disp(['Loaded Truncated raw data file: ' rawfile]);
catch
    cf = pwd; cd(prhloc);
    [rawfile,rawloc]=uigetfile('*.mat', 'select raw data from which prh was made (likely *truncate.mat in prh or raw folder)');
    %         rawfile = [Info.name(1:end-8) '.mat'];
    load([rawloc rawfile],'data');
    cd(cf);
end
ifs = fs;
df = Hzs.datafs/newfs;
iDN = DN;
DN = (DNorig(1):1/newfs/24/60/60:DNorig(end))';
% [Depth,CAL] = pressurecal(data,DN,CAL,false,ofs,df,tagon,Hzs.pHz);
slips = round(slips*newfs/ifs);
Wchange = round(Wchange*newfs/ifs);
for k = 1:length(calperiodI); calperiodI{k} = round(calperiodI{k}*newfs/ifs); end
Depth = decdc(data.Pressure,df);
[fs,Mt_bench,At_bench,Gt,DN,Temp,Light,LightIR,Temp1,tagondec,camondec,audondec,tagslipdec] = decimateandapplybenchcal(data,Depth,CAL,ofs,DN,df,Hzs,tagon,camon,audon,tagslip);
[p,At,Mt,Gt,T,TempI,Light] = applyCal2(data,DN,CAL,camondec,ofs,Hzs,df);
Mt = fir_nodelay(Mt,128,2/(fs/2),'low');
[Aw,Mw,Gw] = applyW(W,slips(1:end-1,2),slips(2:end,1),At,Mt,Gt);
[pitch,roll,head] = calcprh(Aw,Mw,dec);
%     try speed = applySpeed(JigRMS,'JJ',flownoise,'FN',tagondec,p,pitch,roll,fs,speedstats); catch; warning('Speed could not be calculated'); end
% CAL.info = 'Bench cals used for G, 3d in situ cal used for M, A and p. If A3d is empty, bench cal was used. If temp was used in Mag cal, there will be a "temp" variable in the structure; use appycalT to apply that structure to mag and temp data.  Axes must be corrected to NED before applying 3d cals, but not before applying original style bench cals since they take that into account';
tagon = tagondec; camon = camondec; tagslip = slips; 


load([prhloc prhfile],'speed','INFO','speedstats');
ispeed = speed;
if ifs<newfs
        disp('Resampling speed only.  For new sampling rates higher than the original data, this essentially remains at the same sampling rate as the original prh.')
%         disp('To recalculate speed from raw data, rerun script with "true" as the second input (will be slower)');
end
if length(speed.JJ)>length(p)
    speed(length(p)+1:end,:) = [];
else speed.JJ(length(p)) = nan;
end
for sp = {'FN','JJ','SP','JRMS','JJr2','FNr2','section','JJP68','JJP95','JJ95','FNP68','FNP95','FN95'}
    try ispeedFN = ispeed.(sp{1}); catch; continue; end
    speedFN = nan(length(p),size(ispeedFN,2));
    for k = 1:size(ispeedFN,2); oi = resample(timeseries(ispeedFN(:,k),iDN),DN);  speedFN(:,k) = oi.data; end

    if strcmp(sp{1},'section'); speed.(sp{1}) = round(speedFN);
    else speed.(sp{1}) = speedFN;
    end
end
oi = ispeed.section; oi(isnan(oi)) = [];
for i =  unique(oi)'
    Ii = find(ispeed.section==i,1);
    Ii(2) = find(ispeed.section == i,1,'last');
    In = find(speed.section==i,1);
    In(2) = find(speed.section == i,1,'last');
    speed.sectionUsed(In(1):In(2)) = ispeed.sectionUsed(Ii(1));
end
flownoise = nan(1);
load([prhloc prhfile],'flownoise');
flownoise = resample(timeseries(flownoise,iDN),DN); flownoise = flownoise.data;
iJigRMS = JigRMS;
if length(JigRMS.X)>length(p)
    JigRMS(length(p)+1:end,:) = [];
else JigRMS.X(length(p)) = nan;
end
hs = JigRMS.Properties.VariableNames;
for i = 1:size(hs)
    oi = resample(timeseries(iJigRMS.(hs{i}),iDN),DN);
    JigRMS.(hs{i}) = oi.data;
end
speedstats.FN = flownoise;
speedstats.sections_end_index = round(speedstats.sections_end_index*newfs/ifs);

INFO.calperiod = cellfun(@(x) DN(x),calperiodI,'uniformoutput',false);
INFO.calperiodI = calperiodI;
INFO.Wchange = Wchange; %tagslip indices used for tag rotation to animal frame (like Wcalperiods, but trying to find the actual slip)
INFO.tagslip = slips;
try INFO.TempInternal = TempI; catch; end
if exist('nopress','var') && nopress; INFO.NoPressure = true; end % if there's a tag with a messed up pressure sensor

load([prhloc prhfile],'vidDN','vidNam','vidDurs','viddeploy','GPS');
GPSI = find(~isnan(GPS(:,1)));
iGPS = GPS;
GPSIn = round(GPSI*newfs/ifs);
GPS = nan(length(DN),2);
GPS(GPSIn,:) = iGPS(GPSI,:);

nprhfile = [whaleName ' ' num2str(fs) 'Hzprh.mat'];
save([prhloc nprhfile],'Aw','At','Gw','Gt','fs','pitch','roll','head','p','T','Light','Mt','Mw','GPS','DN','speed','speedstats','JigRMS','tagon','camon','vidDN','vidNam','vidDurs','viddeploy','flownoise','INFO','audon');
prhfile = nprhfile;
if exist('Paddles','var'); save([fileloc prhfile],'Paddles','-append'); end
save([INFOloc INFOfile(1:end-8) 'Info' num2str(newfs) 'Hz.mat'],'prhfile','INFO');
VAR = '';
for k = 1:length(VARS); VAR = [VAR ', ''' VARS{k} '''']; end
eval(['save([INFOloc INFOfile(1:end-8) ''Info'' num2str(newfs) ''Hz.mat'']' VAR ',''-append'')']); 

disp('new prh file and INFO file saved successfully')


