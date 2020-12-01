% David Cade
% version 9.21.20
% Goldbogen Lab
% Stanford University

% 1) turn on tag, record times of calibrations and "tag type" in BenchTemplate.xlsx
% It's okay to use separate files if different calibration experiments were
% performed on different days
% 2) run importCATSdata to create a *.mat file for each set of calibration
% data (ideally if all calibrations were performed on one day, this is just
% one file)

%% 1. Test the static accelerometer position
% Inputs: xlsx sheet with calibration times, *.mat file with data
% Outputs: "A" variable with acclerometer values.  To check if axis conventions are correct and units are in m/s^2,
% A variable should be (approx):
% A = [9.8 0 0
% -9.8 0 0
% 0 9.8 0
% 0 -9.8 0
% 0 0 9.8
% 0 0 -9.8];  Eventually we will reduce these to be "1" instead of 9.8.
% Also check that the displayed tagnum, tagtype and fs (sampling rate) are correct
% then check the plots to see if the accelerometer magnitude is close to 1
% throughout the length of the cal file
clear all; close all;
cf = pwd;
[filename,fileloc]=uigetfile('*.mat', 'select the *.mat file with the Acc/Gyro positions');
cd(fileloc);
[benchfile,benchloc] = uigetfile('*.xls*','select xls file with calibration start times');
cd(cf);

data =  load([fileloc filename]); data = data.data;
[benchdata, txtdata] = xlsread([benchloc benchfile]);
%
tagtype = txtdata{3,2}
tagnum = benchdata(1,2)
starts = benchdata(13:3:30,2);
positions = benchdata(13:3:30,1);
stops = benchdata(13:3:30,3);
A = nan(size(positions,1),3); n = 1;
times = data.Time;
fs = round(10/(median(diff(times*24*60*60))))/10 % sampling rate to the nearest tenth of a second, then converted to Hz

acc = [data.Acc1 data.Acc2 data.Acc3];
% add in any missing pieces, sometimes there are gaps in the data
skippeddata = find(diff(times*24*60*60)>1.5*1/fs); % find spots with missing data and replace with nans.  NaNs are NOT accounted for later, so if you have nans, may want to redo data.
for i = length(skippeddata):-1:1
    numnew = length(times);
    times = [times(1:skippeddata(i)); (times(skippeddata(i))+1/fs/24/60/60:1/fs/24/60/60:times(skippeddata(i)+1)-1/fs/24/60/60)'; times(skippeddata(i)+1:end)];
    numnew = length(times)-numnew;
    disp(['WARNING: Missing data of length ' num2str(numnew) ' at time ' datestr(times(skippeddata(i)),'HH:MM:SS.fff')]);
    acc = [acc(1:skippeddata(i),:); nan(numnew,3); acc(skippeddata(i)+1:end,:)];
end

ssI = nan(length(starts),2); %start, stop index in the calibration file for each round
for i = 1:length(starts)
    [~,b] = min(abs(times-starts(i)));
    [~,b2] = min(abs(times-stops(i)));
    ssI(i,:) = [b b2];
end
ssI7 = [ssI(:,1)+7*fs ssI(:,2)-7*fs]; % an index 7 seconds before and after the listed start time (basically wiggle room)
[axA, axM, axG] = axisconventions(tagtype); %fix axes to be NED orientation
acc = acc*axA;
for i = 1:length(ssI7(:,1))
    rows7 = ssI7(i,1):ssI7(i,2);
    A(n,:) = mean(acc(rows7,:));
    if any(std(acc(rows7,:))/max(abs(mean(acc(rows7,:))))>.02);
        disp(['WARNING: Check row ' num2str(n) ', standard deviation of reading = ' num2str(round(max(std(acc(rows7,:))/max(abs(mean(acc(rows7,:)))))*1000)/10) '% of max']);
    end
    n = n+1;
end
A
%
tol = @(B) sum((sqrt(sum(B.^2,2))-1).^2);
calaccCATS = @(inputs,A) tol((A-repmat(inputs(4:6),length(A(:,1)),1))./repmat(inputs(1:3),length(A(:,1)),1));
options = optimset('MaxFunEvals',10^8,'MaxIter',10^8);
[f, feval] = fminsearch(@(inp) calaccCATS(inp,A),[max(A) mean([max(A); min(A)])],options); %if having trouble converting, perhaps need acc cals with axis pointed both towards and away from gravity
gcal = abs(f(1:3)); gconst = f(4:6);

%check gcal against the median mag and the mean mag.  They should be close
%to 1 if there is not a lot of movement.
Gs = repmat(gcal,length(acc),1);
Gc = repmat(gconst,length(acc),1);
M = sqrt(sum(((acc-Gc)./Gs).^2,2));
figure(9); clf; plot(M);
title('Overall magnitude of accelerometer data (should be very close to 1 when the tag is not moving)');
acal = axA*diag(1./gcal);
aconst = (axA*gconst')';
if abs(nanmedian(M)-1)>.01
    warning(['Calibration of accelerometers is inconsistent, median of M = ' num2str(nanmedian(M))]);
end
if abs(nanmean(M)-1)>.02
    warning (['Calibration of accelerometers is inconsistent, mean of M = ' num2str(nanmean(M))]);
end

if strcmpi(tagtype(1:5),'TDR10')
    save([fileloc 'TDR10cal' num2str(tagnum) '.mat'],'acal','aconst');
end

%% 2.  Now enter files where gyro calibration is found (only run this cell if you have gyro data)
% input: enter (or reenter) the mat file and xls file for the gyro cals
% output: gycal gyconst.  Check plots ensure axes line up and values are as
% expected
% Compass data looks like sine waves (with maxima circled to calculate
% rotation rate), Gyro data should be stable with one axis at a high
% positive value and two axes near 0.  If there are a lot of wobbles, redo
% the calibration

clearvars -except fileloc filename tagtype data tagnum benchfile benchloc benchdata acal aconst txtdata

if ~(strcmpi(tagtype,'acousonde') || strcmpi(tagtype(1:5),'TDR10'));
    cf = pwd; if ischar(fileloc); cd(fileloc); end
    [filename2,fileloc2]=uigetfile('*.mat', 'select mat file with gyro calibrations','multiselect','on');
    cd(fileloc2);
    [benchfile2,benchloc2] = uigetfile('*.xls*','select xls file with calibration start times','multiselect','on');
    cd(cf);
    %
    if ~strcmp(filename2,filename); filename = filename2; fileloc = fileloc2;
        data =  load([fileloc filename]); data = data.data;
    end
    if ~strcmp(benchfile2,benchfile); benchfile = benchfile2; benchloc = benchloc2;
        [benchdata, ~] = xlsread([benchloc benchfile]);
    end

        starts = benchdata(13:30,2); %13:30 are the row numbers in the excel file
        positions = benchdata(13:30,1);
        stops = benchdata(13:30,3);
        rpms =benchdata(13:30,4);
else gyconst = nan(1,3); gycal = nan(3,3);
end

times = data.Time;

% newaxis is just an index to look up what each axis is rotated to when the
% tag itself is rotated into a non-normal position
newaxis = [
    1 2 -3
    1 -2 3
    1 3 2
    1 -3 -2
    -1 2 3
    -1 -2 -3
    -1 3 -2
    -1 -3 2
    3 2 1
    3 -2 -1
    3 1 -2
    3 -1 2
    -3 1 2
    -3 -1 -2
    -3 2 -1
    -3 -2 1
    2 3 -1
    2 -3 1
    2 1 3
    2 -1 -3
    -2 3 1
    -2 -3 -1
    -2 1 -3
    -2 -1 3];
%if newz and newy, then newx.  
%
fs = round(10/(median(diff(times*24*60*60))))/10 % sampling rate to the nearest tenth of a second, then converted to Hz

if ~strcmp(tagtype,'TDR10'); comp = [data.Comp1 data.Comp2 data.Comp3];
    if ~strcmp(tagtype,'acousonde'); gyro = [data.Gyr1 data.Gyr2 data.Gyr3]; else gyro = nan(size(comp)); end
end
acc = [data.Acc1 data.Acc2 data.Acc3];

% add in any missing pieces, sometimes there are gaps in the data
skippeddata = find(diff(times*24*60*60)>1.5*1/fs); % find spots with missing data and replace with nans.  NaNs are NOT accounted for later, so if you have nans, may want to redo data.
for i = length(skippeddata):-1:1
    numnew = length(times);
    times = [times(1:skippeddata(i)); (times(skippeddata(i))+1/fs/24/60/60:1/fs/24/60/60:times(skippeddata(i)+1)-1/fs/24/60/60)'; times(skippeddata(i)+1:end)];
    numnew = length(times)-numnew;
    if ~strcmp(tagtype,'TDR10')
        comp = [comp(1:skippeddata(i),:); nan(numnew,3); comp(skippeddata(i)+1:end,:)];
        gyro = [gyro(1:skippeddata(i),:); nan(numnew,3); gyro(skippeddata(i)+1:end,:)];
    end
    acc = [acc(1:skippeddata(i),:); nan(numnew,3); acc(skippeddata(i)+1:end,:)];
    disp(['WARNING: Missing data of length ' num2str(numnew) ' at time ' datestr(times(skippeddata(i)),'HH:MM:SS.fff')]);
    
end
if ~strcmp(tagtype,'acousonde')&& ~strcmp(tagtype,'TDR10'); % take out
    ssI = nan(length(starts),2); %start, stop index in the calibration file for each round
    for i = 1:length(starts)
        [~,b] = min(abs(times-starts(i)));
        [~,b2] = min(abs(times-stops(i)));
        ssI(i,:) = [b b2];
    end
    ssI5 = [ssI(:,1)+5*fs ssI(:,2)-5*fs]; % an index 5 seconds before and after the listed start time (basically wiggle room)
end
%
[axA, axM, axG] = axisconventions(tagtype); %fix axes to be NED orientation
acc = acc*axA;
% sets up acc calibration to use to calibrate gyros based on subtracted
% centripetal acceleration when axes are off

gcal = (1./diag(axA^-1*acal))';
gconst = (axA^-1*aconst')';
Gs = repmat(gcal,length(acc),1);
Gc = repmat(gconst,length(acc),1);

if ~strcmp(tagtype,'TDR10')
    comp = comp*axM;
    gyro = gyro*axG;
end

%
% to calibrate gyroscopes, need the orientation information of the axes.
% If it is tilted, the rotation rate(calculated from the magnetometers) must be adjusted.
gycal = cell(length(unique(positions)),1); gyconst = gycal;
figure(10); clf;
k = ceil(length(rpms(rpms>0))/3); kk = 0;
for pos = unique(positions)'
    % figure out orientation of the position, rotation to be normal face
    % forward position (makes calculations easier)
    row = find(positions==pos & rpms == 0);
    rows5 = ssI5(row,1):ssI5(row,2); % the rows 5 seconds after start and before stop
    rows = ssI(row,1):ssI(row,2);
    %tells you which way the axis mostly facing gravity is oriented
    [~,newz] = max(abs(mean(acc(rows5,:))));
    direction = sign(mean(acc(rows5,newz))); % sign should be negative in R hand orientation, else flip it
    newz = -direction*newz;
    row2 = find(positions==pos & rpms >0,1,'first');
    rows52 = ssI5(row2,1):ssI5(row2,2);
    [~,newy] = max(abs(mean(acc(rows52,:))-mean(acc(rows5,:)))); % the biggest difference accel is the axis facing towards the rotation gyro axis on the record player
    newy = sign(mean(acc(rows52,newy))-mean(acc(rows5,newy)))*newy; % normal direction is facing towards
    if abs(newz) == abs(newy); error(['error in determining orientation of calibration position ' num2str(pos)]); end
    newx = newaxis(newaxis(:,1) == newz & newaxis(:,2) == newy,3);
    gy = gyro(:,abs([newx newy newz])).*repmat(sign([newx newy newz]),length(gyro(:,1)),1);
    co = comp(:,abs([newx newy newz])).*repmat(sign([newx newy newz]),length(gyro(:,1)),1);
    ac = (acc(:,abs([newx newy newz])).*repmat(sign([newx newy newz]),length(gyro(:,1)),1)-Gc(:,abs([newx newy newz])))./Gs(:,abs([newx newy newz])); %ac now is calibrated
    gy = filterCATS(gy,round(fs/8),round(fs/4),.5); % gets rid of random spikes (could use smaller threshold?  other times we use .05)
    co = filterCATS(co,round(fs/8),round(fs/4),.5);
    ac = filterCATS(ac,round(fs/8),round(fs/4),.5);
    
    p = asin(mean(ac(rows5,1))); % positive values indicate whale is facing up. (positive pitch). see rotateVr
    r = atan2(-mean(ac(rows5,2)),-mean(ac(rows5,3))); %negative gets carried through.  Important to get the right quadrant
    
    oigy = gy; oico = co;
    for ii = 1:3
        oigy(:,ii) = runmean(inpaint_nans(gy(:,ii)),floor(fs/12)); n = 0;
        oico(:,ii) = runmean(inpaint_nans(co(:,ii)),floor(fs/5)); % compass needs smoothing to ensure a similar spot is found each revolution
    end
    
    % filter magnetometer (old, no longer needed with filter cats, but need to figure out how to remove it)
    baddiff = nan(3,1); magmax = baddiff;
    for ax = 1:3
        [loc1, mag1] = peakfinder(oico(rows52,ax)); % can use .9*(max(oico(rows52,ax))-min(oico(rows52,ax)))+min(oico(rows52,ax)) as a threshold if axis was off kilter
        [loc2, mag2] = peakfinder(oico(rows52,ax),[],[],-1);
        magmax(ax) = median(mag1);
        baddiff(ax) = (median(mag1)-median(mag2))/2; %a threshold to look for filtering the data
        co(find(abs(diff(co(:,ax)))>baddiff(ax))+1,ax) = nan; % gets rid of weird spikes
        co(:,ax) = inpaint_nans(co(:,ax));
    end
    Omag = nanmean(co(rows5,:)); % the magnetometer readings during this orientation
    % get all the spots in the gyro rotations
    nrpm = 0; gycal{pos} = nan(length(rpms(rpms>0&positions==pos))*2,3);
    for rpm = rpms(rpms>0 & positions == pos)'
        nrpm = nrpm+1;
        rowm = find(positions==pos & rpms == rpm,1,'first');
        rows5m = ssI5(rowm,1):ssI5(rowm,2);
        
        samps = nan(3,1); % the number of samples between peaks (basically how many samples between revolutions)
        kk = kk+1; subplot(3,k,kk); hold on;
        [~,ax] = min(std(co(rows5m,:))); % should be 3, finds the axis that is rotated the least and removes it as a gauge of how quickly the tag is spinning.
        pax = 1:3; pax = pax(pax~=ax);
        for ax = pax
            diffs = diff(peakfinderplot(oico(rows5m,ax),(max(oico(rows5m,ax))-min(oico(rows5m,ax)))/2,magmax(ax)-.4*baddiff(ax)));
            if isempty(diffs); diffs = diff(peakfinderplot(oico(rows5m,ax),(max(oico(rows5m,ax))-min(oico(rows5m,ax)))/2)); end
            if isempty(diffs); error('Peakfinder failed to detect peaks in mag data'); end
            diffs(diffs>1.5*mean(diffs)) = []; % ensures no peaks (rotations of magnetometer) are skipped in calibration
            samps(ax) = mean(diffs); % this used to not be the smoothed compass, whiched seemed silly.
        end
        %         disp([samps std(co(rows5m,:))']);
        samp = round(nanmean(samps)); % one will be a nan
        if any(abs(samps-samp)>0.05*samp); error('Check rotation rate calculations, compass peaks do not line up.'); end
        
        %         plot(oico(rows5m,pax),'g');
        [ax,h1,h2] = plotyy(1:length(rows5m),oico(rows5m,pax),1:length(rows5m),oigy(rows5m,:)); set(h2,'linewidth',3); set(h1,'color','g');
        set(ax(1),'yticklabel',[]);
        ys = get(gca,'ylim'); xs = get(gca,'xlim');
        text(xs(2),ys(1),[num2str(pos) '-' num2str(rpm)],'verticalalignment','bottom','horizontalalignment','right');
        %         clear oig;
        gyval = nan(ceil(length(rows5m)/samp),3);
        % find the gyro values every time the tag is back in the position
        % from which stillness was measured
        for s = rows5m(1):samp:rows5m(end)
            [~,b] = min(sum((oico(s:s+samp-1,:)-repmat(Omag,samp,1)).^2,2));
            n = n+1;
            gyval(n,:) = mean(oigy(s+b-1-fs/10:s+b+fs/10-1,:)); %1/10 of a second to each side of the smoothed gyro val.
        end
        rotrate = 1/(samp/fs)*2*pi; % radians per second
        %         H = [cos(h) -sin(h) 0; sin(h) cos(h) 0; 0 0 1];
        P = [cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)];
        R = [1 0 0; 0 cos(r) -sin(r); 0 sin(r) cos(r)];
        gyconst{pos} = mean(gy(rows5,:)); %when rpm = 0 in this position.
        oi = (mean(gyval)-gyconst{pos})*R'*P'/rotrate; %either rotate the values to the rotation axis and divide by rotrate
        gycal{pos}(nrpm*2-1,abs(newz)) = oi(3);
        oi = (mean(gyval)-gyconst{pos})./([0 0 rotrate]*P*R); % or rotate the rotrate to the tag axes
        gycal{pos}(nrpm*2,abs(newz)) = oi(3); % three was the rotated axis in this frame
        gyconst{pos}(abs([newx newy newz])) = gyconst{pos}.*sign([newx newy newz]);
    end
end
title({'Left axis is magnetometer data (peaks are used to calculate rotation rate),'; 'right axis is gyros (axis of interest should be flat, other two should be ~ 0, 45 rpm should be higher than 33 rpm'});

gycal = 1./nanmean(vertcat(gycal{:}));
I = [1 0 0; 0 1 0; 0 0 1];
gyconst = mean(vertcat(gyconst{:}));

gyconst = (axG*gyconst')'
gycal = axG*diag(gycal)

if sum(sum(isnan(gyconst))) == length(gyconst(:)); disp('No Gyros found'); gycal = nan(3); gyconst = nan(1,3); end

%
%% 3. magnetometers (load new files)
% Ensure calibration period is properly identified with no spikes in data
% that could mess up calibraton
% input: enter (or reenter) the mat file and xls file for the gyro cals
% output: magcalon magconston, magcaloff magconstoff. These are initial
% calibrations and represent an alternate method to the spherical
% calibrations that end up being employed.
% In figure 8 black line should be pretty close to upper red line.  All other axes
% should be well constrained by red lines.

cf = pwd; if ischar(fileloc); cd(fileloc); end
[filename2,fileloc2]=uigetfile('*.mat', 'select magnetometer data file');
cd(fileloc2);
[benchfile2,benchloc2] = uigetfile('*.xls*','select xls file with calibration start times');
cd(cf);
%
    if ~strcmp(filename2,filename); filename = filename2; fileloc = fileloc2;
        data =  load([fileloc filename]); data = data.data;
    end
    if ~strcmp(benchfile2,benchfile); benchfile = benchfile2; benchloc = benchloc2;
        [benchdata, txtdata] = xlsread([benchloc benchfile]);
    end
%

    starts = benchdata(17:18,8);
    %     positions = benchdata(~isnan(benchdata(:,13)),12);
    stops = benchdata(17:18,9);
    cameras = txtdata(18:19,7);
  times = data.Time;

%
timescal = data.Date + data.Time;

fs = round(10/(median(diff(times*24*60*60))))/10  % sampling rate to the nearest hundredth of a second, then converted to Hz

comp = [data.Comp1 data.Comp2 data.Comp3];
if ~strcmp(tagtype,'acousonde'); gyro = [data.Gyr1 data.Gyr2 data.Gyr3]; else gyro = nan(size(comp)); end
acc = [data.Acc1 data.Acc2 data.Acc3];

% add in any missing pieces, sometimes there are gaps in the data
skippeddata = find(diff(times*24*60*60)>1.5*1/fs); % find spots with missing data and replace with nans.  NaNs are NOT accounted for later, so if you have nans, may want to redo data.
for i = length(skippeddata):-1:1
    numnew = length(times);
    times = [times(1:skippeddata(i)); (times(skippeddata(i))+1/fs/24/60/60:1/fs/24/60/60:times(skippeddata(i)+1)-1/fs/24/60/60)'; times(skippeddata(i)+1:end)];
    timescal = [timescal(1:skippeddata(i)); (timescal(skippeddata(i))+1/fs/24/60/60:1/fs/24/60/60:timescal(skippeddata(i)+1)-1/fs/24/60/60)'; timescal(skippeddata(i)+1:end)];
    
    numnew = length(times)-numnew;
    comp = [comp(1:skippeddata(i),:); nan(numnew,3); comp(skippeddata(i)+1:end,:)];
    gyro = [gyro(1:skippeddata(i),:); nan(numnew,3); gyro(skippeddata(i)+1:end,:)];
    acc = [acc(1:skippeddata(i),:); nan(numnew,3); acc(skippeddata(i)+1:end,:)];
     disp(['WARNING: Missing data of length ' num2str(numnew) ' at time ' datestr(times(skippeddata(i)),'HH:MM:SS.fff')]);

end
comp0 = comp;
ssI = nan(length(starts),2); %start, stop index in the calibration file for each round
for i = 1:length(starts)
    [~,b] = min(abs(times-starts(i)));
    [~,b2] = min(abs(times-stops(i)));
    ssI(i,:) = [b b2];
end
ssI7 = [ssI(:,1)+7*fs ssI(:,2)-7*fs]; % an index 7 seconds before and after the listed start time (basically wiggle room)

comp = comp*axM;

gyro = (gyro-ones(size(gyro))*diag(gyconst))*gycal; % use calibration from earlier. Gives radians/sec
acc = (acc-ones(size(acc))*diag(aconst))*acal; % gives value in gs

% # of offs must equal # of ons, usually just one of each
ons = find(strcmp(cameras,'on') == 1);
offs = find(strcmp(cameras,'off') == 1);
maxs = cell(1,2);
maxs{1} = nan(length(ons),3); maxs{2} = maxs{1}; mins = maxs;
%
j = 1;
    fig = figure(1); clf;
    XYZ = 'XYZ';
    if j ~=1; set(fig, 'windowStyle','docked'); end
    ind = min(ssI([ons(j) offs(j)],1)):max(ssI([ons(j) offs(j)],2));
    indoffon = {ssI7(offs(j),1):ssI7(offs(j),2); ssI7(ons(j),1):ssI7(ons(j),2)};
    inds = cell(1,2); inds{1} = nan(2,3); inds{2} = inds{1}; indsI = inds;
    hold on;
    ys = [min(min(comp)) max(max(comp))];
          clear ts; clear cs;
    cord = get(gca,'colororder');
    offon = {'OFF' 'ON'};
    %
    for k = 1:2
        if k == 2 && starts(k) == starts(1); continue;  end
        xs = [max(1,indoffon{k}(1)-30*fs) min(length(timescal),indoffon{k}(end)+30*fs)];
        xlim(timescal(xs));
        oi = datestr(get(gca,'xtick'),'HH:MM:SS');
        set(gca,'xticklabel',oi);
        for i = 1:3
            ys = [min(comp(xs(1):xs(2),i)) max(comp(xs(1):xs(2),i))]; ys(1) = ys(1)-diff(ys)/5; ys(2) = ys(2)+diff(ys)/5;
            %             ys = [-1000 1000];
            ylim(ys);
            p1 = cell(0);
              for ii = 1:length(p1); try delete (p1{ii});catch; end; end;  clear p1;
             p1{1} = patch(timescal([ssI7(offs(j),:) ssI7(offs(j),[2 1])])',[ys(1) ys(1) ys(2) ys(2)],[255 255 100]/255);
            if k == 2;  for ii = 1:length(p2); try delete (p2{ii}); catch; end; end; clear p2; p2{1} = patch(timescal([ssI7(ons(j),:) ssI7(ons(j),[2 1])])',[ys(1) ys(1) ys(2) ys(2)],[175 238 238]/255); end;
            title({'Press enter if no spikes and boundaries capture rotations,' 'else click start and end of rotation periods, excluding any abnormal spikes' '(must be an even number of clicks)'});
            ts = text(timescal(indoffon{k}(1)-30*fs),min(get(gca,'ylim')),[XYZ(i) ' cam ' offon{k} ' boundaries'],'verticalalignment','bottom');
            cs = plot(timescal(ind),comp(ind,i),'color',cord(i,:));
            goodbs = false;
            inds{k}(:,i) = timescal([indoffon{k}(1); indoffon{k}(end)]); indsI{k}(:,i) = [indoffon{k}(1); indoffon{k}(end)];
            while ~goodbs
                [x,~] = ginput();
                if floor(length(x)/2) + .1 < length(x)/2; tt = text(timescal(indoffon{k}(1)-30*fs),max(get(gca,'ylim')),'start and end boundaries must be paired'); continue;
                else try clear tt; catch; end; end;
                if isempty(x); goodbs = true; continue;
                else inds{k}(1:length(x),i) = x; 
                    for ii = 1:length(x); 
                        [~,indsI{k}(ii,i)] = min(abs(timescal-inds{k}(ii,i))); 
                    end
                end
                if k == 1;
                    for ii = 1:length(p1); delete (p1{ii}); end; clear p1
                    for ii = 1:2:length(indsI{k}(:,i))
                        if isnan(indsI{k}(ii,i)); continue; end
                        p1{ii} = patch(timescal([indsI{k}(ii,i) indsI{k}(ii+1,i) indsI{k}(ii+1,i) indsI{k}(ii,i)])',[ys(1) ys(1) ys(2) ys(2)],[255 255 100]/255);
                    end
                elseif k == 2;
                    for ii = 1:length(p2); delete (p2{ii}); end; clear p2
                    for ii = 1:2:length(indsI{k}(:,i))
                        if isnan(indsI{k}(ii,i)); continue; end
                        p2{ii} = patch(timescal([indsI{k}(ii,i) indsI{k}(ii+1,i) indsI{k}(ii+1,i) indsI{k}(ii,i)])',[ys(1) ys(1) ys(2) ys(2)],[255 255 100]/255);
                    end
                end
                delete([cs;ts]);
                  ts = text(timescal(indoffon{k}(1)-30*fs),min(get(gca,'ylim')),[XYZ(i) ' cam ' offon{k} ' boundaries'],'verticalalignment','bottom');
                  cs = plot(timescal(ind),comp(ind,i),'color',cord(i,:));
             end
             delete([cs;ts]);
        end
        indsI{k}(indsI{k}==0) = nan;
    end
    
    %     delete([p1 p2 ts]);
    for k = 1:2
        if k == 2 && starts(k) == starts(1); continue;  end
        figure(k);
        clf; s = nan(3,1);
        for i = 1:3
            s(i) = subplot(3,1,i); hold on;
            title(['Cam ' offon{k}])
            ys = [min(min(comp)) max(max(comp))];
            clear p1
            
            for ii = 1:2:length(indsI{k}(:,i))
                if isnan(indsI{k}(ii,i)); continue; end
                p1{ii} = patch(timescal([indsI{k}(ii,i) indsI{k}(ii+1,i) indsI{k}(ii+1,i) indsI{k}(ii,i)])',[ys(1) ys(1) ys(2) ys(2)],[255 255 100]/255);
            end
            cs = plot(timescal(ind),comp(ind,i),'color',cord(i,:));
            xlim(timescal([ind(1) ind(end)]));
            ylim([min(comp(indsI{k}(1,i):indsI{k}(find(~isnan(indsI{k}(:,1)),1,'last'),i),i)) max(comp(indsI{k}(1,i):indsI{k}(find(~isnan(indsI{k}(:,1)), 1, 'last' ),i),i))]);
            text(timescal(ind(1)),max(get(gca,'ylim'))-.1*diff(get(gca,'ylim')),XYZ(i),'fontweight','bold','fontsize',18);
            oi = datestr(get(gca,'xtick'),'HH:MM:SS');
            set(gca,'xticklabel',oi);
        end
        
    end
%     legend([p1; p2],'Cam off', 'Cam on','orientation','horizontal','location','south');
    
    % ylim([min(min(comp(ind,:)))*1.1 max(max(comp(ind,:)))*1.1]);
    
    for k = 1:2
        for i = 1:3
            for ii = 1:2:length(indsI{k}(:,1))
                j = (ii+1)/2;
                if isnan(indsI{k}(ii,i)); maxs{k}(j,i) = nan; else
            maxs{k}(j,i) = max(comp(indsI{k}(ii,i):indsI{k}(ii+1,i),i)); 
            mins{k}(j,i) = min(comp(indsI{k}(ii,i):indsI{k}(ii+1,i),i));
                end
            end
        end
        maxs{k} = max(maxs{k});
        mins{k} = min(mins{k});
    end
    
clear GPS elev
GPS = benchdata(3:4,2);
elev = benchdata(5,2);

try 
    if DV(1,1)<2015; str = '2010'; elseif DV(1,1)<2020; str = '2015'; else str = '2020'; end
    [~,~,~,~,b] = wrldmagm(elev,GPS(1),GPS(2),decyear(data.Date(1)),str); % Best guess for 2015 dates until they update the script
catch
 disp('Go to https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm and enter the magnetic field strength (in nT) for the place where calibration was performed')
    b = input('? ');
end
b = b*10^-3;
Mx = maxs{1}; Mn = mins{1};
foff(1:3) = 2*b./(Mx-Mn);
foff(4:6) = (Mx+Mn)/2;
Mx = maxs{2}; Mn = mins{2};
fon(1:3) = 2*b./(Mx-Mn);
fon(4:6) = (Mx+Mn)/2;

% Test the original calibrations

magcaloff = axM*diag(foff(1:3));
magconstoff = (axM*foff(4:6)')';
if starts(2) == starts(1); 
    magcalon = magcaloff; magconston = magconstoff;
else magcalon = axM*diag(fon(1:3)); magconston = (axM*fon(4:6)')';
end

compToff = (filterCATS(comp0,ceil(fs/8),round(fs),.05)-repmat(magconstoff,size(comp0,1),1))*magcaloff;
compTon = (filterCATS(comp0,ceil(fs/8),round(fs),.05)-repmat(magconston,size(comp0,1),1))*magcalon;

figure(8); clf;
s(1) = subplot(2,1,1);
ind = min(ssI(:,1)):max(ssI(:,2));
plot(sqrt(sum(compToff.^2,2)),'k','linewidth',2);hold on;
plot(compToff);
plot([0 length(compToff)],[b b],'r',[0 length(compToff)],[-b -b],'r');
legend('Magoff','Xoff','Yoff','Zoff','orientation','horizontal','location','southeast');
xlim([ind(1) ind(end)]);
ylim([-60 60]);
for j = 1:length(offs)
    p1 = patch([ssI7(offs(j),:) ssI7(offs(j),[2 1])]',[-900 -900 900 900],[255 255 100]/255,'facealpha',.6);
end
s(2) = subplot(2,1,2);
ind = min(ssI(:,1)):max(ssI(:,2));
plot(sqrt(sum(compTon.^2,2)),'k','linewidth',2);hold on;
plot(compTon);
plot([0 length(compTon)],[b b],'r',[0 length(compTon)],[-b -b],'r');
legend('Magon','Xon','Yon','Zon','orientation','horizontal','location','southeast');
xlim([ind(1) ind(end)]);
ylim([-60 60]);
for j = 1:length(ons)
    p2 = patch([ssI7(ons(j),:) ssI7(ons(j),[2 1])]',[-900 -900 900 900],[175 238 238]/255,'facealpha',.6);
end
%% 4. Add in spherical calibrations (Mark Johnson method) for magnetometer and accelerometer

At = acc; % these are already calibrated 
Mt = comp;
try camon = data.Camera>0; 
catch; camon = false(size(data.Acc1));
warning('cam on periods not automatically detected, for display, assuming cam is off')
end



for k = 1:2
    if k == 2 && starts(2) == starts(1);  Acal3d0cam =  Acal3d0;   Mcal3d0cam = Mcal3d0; continue; end
    figure (4+k); clf
    I = false(size(At(:,1))); % make an index that uses all the previously selected sections
    for ii = 1:2:length(indsI{k}(:,1));
        for i =1:3;
            if isnan(indsI{k}(ii,i)); continue; end
            I(indsI{k}(ii,i):indsI{k}(ii+1,i)) = true;
        end
    end
    s1 =  subplot(211);
    camstarts = find(diff(camon) == 1); if camon(1); camstarts = [1; camstarts]; end
    camends = find(diff(camon) == -1); if camon(end); camends = [camends; length(camon)]; end
    for ii = 1:length(camstarts)
        camp = patch([camstarts(ii) camends(ii) camends(ii) camstarts(ii)],[-1000 -1000 1000 1000],[10 10 10]/255,'facealpha',.1);
    end
    hold on;
    calstarts = find(diff(I) == 1); if I(1); calstarts = [1; camstarts]; end
    calends = find(diff(I) == -1); if I(end); calends = [camends; length(camon)]; end
    if k == 1; col = [255 255 100]/255;
    else col = [175 238 238]/255;
    end
    for ii = 1:length(calstarts)
        calp = patch([calstarts(ii) calends(ii) calends(ii) calstarts(ii)],[-1000 -1000 1000 1000],col);
    end
      disp('ideal is low residual (< 5 %) and high balance (> 20%)');
    [~,Acalnew,~] = spherical_calwk(At(I,:),1);
    
    Atnew = (At*diag(Acalnew.poly(:,1))+repmat(Acalnew.poly(:,2)',size(At,1),1))*Acalnew.cross;

    ap = plot(Atnew);  hold on; am = plot(sqrt(sum(Atnew.^2,2)),'k');
    ylim([-2 2])
    if k == 1
        Acal3d0 = Acalnew;
    elseif k == 2
        Acal3d0cam = Acalnew;
    end
    try legend([camp calp ap' am],'Cam on','Calibration period','X','Y','Z','Magnitude','position','Eastoutside')
    catch; legend([calp ap' am],'Calibration period','X','Y','Z','Magnitude','position','Eastoutside')
    end
        title (['Cam ' offon{k} ' cal'],'fontsize',14);
    
       [~, Mcalnew,~] = spherical_calwk(Mt(I,:),b,'cross');
    Mtnew = (Mt*diag(Mcalnew.poly(:,1))+repmat(Mcalnew.poly(:,2)',size(Mt,1),1))*Mcalnew.cross;
    s2 = subplot(212); hold on
     for ii = 1:length(camstarts)
        camp = patch([camstarts(ii) camends(ii) camends(ii) camstarts(ii)],[-1000 -1000 1000 1000],[10 10 10]/255,'facealpha',.1);
    end
    for ii = 1:length(calstarts)
        calp = patch([calstarts(ii) calends(ii) calends(ii) calstarts(ii)],[-1000 -1000 1000 1000],col);
    end
    plot(Mtnew);  hold on; plot(sqrt(sum(Mtnew.^2,2)),'k');
    ylim([-100 100])
    if k == 1
        Mcal3d0 = Mcalnew;
    elseif k == 2;
        Mcal3d0cam = Mcalnew;
    end
     try legend([camp calp ap' am],'Cam on','Calibration period','X','Y','Z','Magnitude','position','Eastoutside')
     catch; legend([calp ap' am],'Calibration period','X','Y','Z','Magnitude','position','Eastoutside')
     end
         linkaxes([s1 s2],'x')
end

%% 5 save cal file
% Can enter factory or other calibrations here for pressure, Temperature
% etc.
% for older plug in tags:
% pcal = 100/32768/.0980665; pconst = 16384;
% Tconst = - 21*333.87; Tcal = 1/333.87;
% acal =   1 / 8192  ; % * 9.80665; % last term for gravity
% magcal = .15;
% gycal = ( 1 / 32.8 ) * ( 3.141592653 / 180 );

str = 'CATS'; % can set a prefix (e.g. AcousondeCal32;
if ~exist('Tcal','var'); Tcal = 1; Tconst = 0; end
if ~ exist('pcal','var'); pcal = 1; end
pconst = nanmean(data.Pressure(1:180*fs)); %give pressure above water in first few minutes (before tag gets hot)

% if calibration files are in a "bench test" folder, this puts the cal
% files one level up;
dirup = max(strfind(fileloc(1:strfind(fileloc,'Bench')),'\'));
if isempty(dirup); dirup = length(fileloc); end

save([fileloc(1:dirup) str 'cal' num2str(tagnum) '.mat'],'gycal','gyconst','acal','aconst','magcaloff','magconstoff','magcalon','magconston','pcal','pconst','Tcal','Tconst');
if exist('Acal3d0','var') && ~isempty(Acal3d0)
    save([fileloc(1:dirup) str 'cal' num2str(tagnum) '.mat'],'Acal3d0','Acal3d0cam','-append');
end
if exist('Mcal3d0','var') && ~isempty(Mcal3d0)
    save([fileloc(1:dirup) str 'cal' num2str(tagnum) '.mat'],'Mcal3d0','Mcal3d0cam','-append');
end
