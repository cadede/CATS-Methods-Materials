function [GPS,GPSerr,UTM] = GPS2prh (DN,p,fs,tagon,data,prhname,maxDist,GPSrate,timedif)
% DN,p,fs,tagon from prh file, data is the raw data with GPS info
% maxDist is optional and is the maximum distance in km the whale could
% have traveled (gets rid of points furhter than this
% GPSrate is the downsampled GPS rate (in seconds)

if nargin<8; GPSrate = 10; end % 10 means find the best GPS point in every 10 second window (once per surfacing).
if nargin<7 || isempty(maxDist); maxDist = 50000; end %~30 miles

%
dataHz = round(1/(nanmedian(diff(data.Time(1:101,1)))*24*60*60));
if dataHz>fs % resample down to 
    data = data(1:dataHz/fs:end,:);
end

sd = sqrt(sum(data.GPSsd.^2,2)); %mag of all 2d std dev values.  Not sure if this is best metric or not
csd = nanstd(sd);

t1 = DN(tagon); t2 = t1(end); t1 = t1(1);
ts5 = find(data.Date+data.Time>t1+5/24/60&data.Date+data.Time<t2-5/24/60); %use time during deployment to reduce risk of weird spieks
if isempty(ts5); ts5 = find(tagon); end 
% find time offset (from data file to prh file, should only be if tag was set incorrectly)
Tp = data.Pressure;
if sum(Tp) == 0 || sum(Tp==30) == length(Tp); %Tp = load([fileloc '\' F{j}{~cellfun(@isempty,cellfun(@(x) regexp(x,'prh'),F{j},'UniformOutput',false))}],'p','DN','INFO','fs'); Tp.timedif = Tp.INFO.timedif;
    TpDN = DN-timedif/24; [~,a1] = min(abs((data.Date+data.Time)-TpDN(1)));
    dfs = round(1./mean((data.Time(50:60)-data.Time(49:59))*24*60*60)); a2 = a1+dfs/fs*(length(p))-1;
    tTp = nan(size(data.Time)); tTp(a1:dfs/fs:a2) = p; tTp = fixgaps(tTp); tTp([1:a1-1 a2+1:end]) = 0;
    Tp = tTp; clear tTp;
end
data.Pressure = Tp;
try
    [pI,pH] = peakfinder(p,(max(p)-min(p))/4,(max(p)-min(p))/2+min(p));
    [PI,PH] = peakfinder(data.Pressure,(max(data.Pressure)-min(data.Pressure))/4,(max(data.Pressure)-min(data.Pressure))/2+min(data.Pressure));
    if length(pI)~=length(PI); numP = min([length(pI),length(PI),5]);
        [~,pSI] = sort(pH,'descend'); [~,PSI] = sort(PH,'descend');
        pI = pI(pSI(1:numP)); PI = PI(PSI(1:numP));
    end
    DNp = DN(pI);
    DNP = data.Date(PI)+data.Time(PI);
    offset = mean(DNP-DNp)*24;
    roffset = round(mean(DNP-DNp)*24*2)/2;
    if abs(offset-roffset)*60*60>20;  error('yup'); end
    offset = roffset;
catch
    try offset = -INFO.timedif;
    catch
        disp('Minimum dive depth is different in data and prh file, find another metric');
        offset = -input(['Time offset (from excel file ' prhname(1:end-12) '.xls in metadata folder) = ']);
    end
end

disp(['Time offset used for ' prhname ': ' num2str(-offset)]);

GPSDN = data.Date+data.Time-offset/24;

GPSdata = [data.Lat data.Long];
GPS = nan(length(DN),2);
GPSerr = nan(length(DN),3);
% kDN = GPSDN(find(GPSDN>t1-2/24/60,1)); %k1 = find(GPSDN == kDN,1);
% oi = find(~isnan(GPSdata(GPSDN<t1,1)),fs*60,'last');
% oi = oi(1);%ensure there are some GPS points pre deployment
% kDN = max(min(GPSDN(find(GPSDN>t1-2/24/60,1)),GPSDN(oi)),DN(1)); % find the first spot 2 minutes before tag on 
kDN = max(DN(1),GPSDN(1));
clat = nanmedian(data.Lat); clong = nanmedian(data.Long);
% pastt2 = 0;
while kDN<min(GPSDN(end),DN(end))
    k1 = find(GPSDN>=kDN,1);
    k1 = find(~isnan(GPSdata(k1:end,1)),1)+k1 -1;
    if isempty(k1); break; end
    k2 = find(GPSDN<GPSDN(k1)+GPSrate/24/60/60,1,'last'); % finds all the GPS points in 10 second groups
    tGPS = GPSdata(k1:k2,:);
    DIST = arrayfun(@(x,y) distance(x,y,clat,clong),tGPS(:,1),tGPS(:,2));
    DIST = arrayfun(@(x) distdim(x,'degrees','m'),DIST);
    %if any(DIST>1000); disp(max(DIST)); end
    tGPS(DIST>maxDist,:) = nan;
    tGPSerr = data.GPSsd(k1:k2,:);
    if sum(~isnan(tGPS(:,1))) > 5
        DIST = arrayfun(@(x,y) distance(x,y,nanmedian(tGPS(:,1)),nanmedian(tGPS(:,2))),tGPS(:,1),tGPS(:,2));
        DIST = arrayfun(@(x) distdim(x,'degrees','m'),DIST);
        tGPS(DIST>mean(DIST)+std(DIST),:) = nan; %gets rid of points further than 1 std from the mean distance from the median
    end
    tsd = sd(k1:k2); % now find the best sds
    tGPS(tsd>min(tsd)+csd,:) = nan;
    gps = nanmean(tGPS);
    kDN = GPSDN(min(k2+1,length(GPSDN)));
    kk2 = find(DN<=kDN,1,'last');
    [~,b] = min(p(kk2-GPSrate*fs+1:kk2)); %put it at the surfacing of this 10 second interval
    b = b(1)+kk2-GPSrate*fs+1 -1;
    GPS(b,:) = gps;
    tGPSerr(isnan(tGPS(:,1)),:) = nan;
    GPSerr(b,:) = nanmean(tGPSerr);
%     if kDN>t2; pastt2 = pastt2+GPSrate; end  
%     if pastt2>=120; break; end % if you have two minutes of GPS data post recovery
end

if nargout == 3
    UTM = nan (size(GPS));
    [x,y] = deg2utm(GPS(~isnan(GPS(:,1)),1),GPS(~isnan(GPS(:,1)),2));
    UTM(~isnan(GPS(:,1)),:) = [x y];
end

I = find(~isnan(GPSdata(:,1)),1,'last');
if GPSDN(I)>t2; GPS(end,:) = GPSdata(I,:); GPSerr(end,:) = data.GPSsd(I,:); end %put the last known GPS at the end of the prh file
                