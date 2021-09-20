function data = importOPENTAGdata
%%
% David Cade
% Nov 2018 (may need adjustments for more recent data)
%
%
% when there are gaps in the data, which happens at the end of every file,
% this code assumes that it is a gap and not a sampling error.
% it assumes that sensors are sampled evenly and discretely, and there is error
% checking included to ensure this.
    pcal = 1/111.377; %manufacturer calibration
    pconst = 1010;

    fileloc = uigetdir('select folder with INER and PTMP csv files');
    fileloc = [fileloc '\'];
    D = dir(fileloc);
    D = {D.name};
    D = D(~cellfun(@isempty,cellfun(@(x) strfind(x,'.csv'),D,'uniformoutput',false)))';
    dates = cellfun(@(x) datenum(x([10:17 23:31]),'HH_MM_SS_dd_mm_yy'),D);
    [dates,I] = sort(dates);
    D = D(I);
    warning('off','all');
    dataI = readtable([fileloc D{1}],'readvariablenames',true);
    dataP = readtable([fileloc D{2}],'readvariablenames',true);
    times = datenum(dataI.FileTime,'HH:MM:SS'); times = times-floor(times);
    DN = dates(1) + dataI.TimeFromStart_s_/24/60/60; UTC = (dates(1) - datenum(dataI.FileDate(1),'mm/dd/yy') - times(1))*24;
%     DN = datenum(dataI.FileDate,'mm/dd/yy')+times; UTC = (dates(1)-DN(1))*24;
    if abs(UTC-round(UTC))>.01; error('check UTC offset'); else UTC = round(UTC); disp(['UTC = ' num2str(UTC)]); end
    fs = 1/mean(diff(dataI.TimeFromStart_s_)); disp(['fs = ' num2str(fs)]);
    if abs(fs-round(fs))>.01; error('skipped data'); end
%     nexthour = ceil(DN(10)*24)/24-1/fs/24/60/60;
    nextfile = dates(3);
    DN = (DN(1):1/fs/24/60/60:nextfile-1/fs/24/60/60)';
    if length(DN)<size(dataI,1);  %dataI = dataI(1:length(DN),:);
        error('data extends beyond next file start'); 
    end
    UTCDN = DN-UTC/24;
    data = table (UTCDN,floor(DN),DN-floor(DN),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),...
        'VariableNames',{'UTCDN','Date','Time','Acc1','Acc2','Acc3','Comp1','Comp2','Comp3','Gyr1','Gyr2','Gyr3','Pressure','Temp'});
    b = size(dataI,1);
%     data.Date = floor(DN + UTC/24); data.Time = DN+UTC/24-floor(DN+UTC/24); %DN+UTC/24 - data.Date + dataI.dataI.TimeFromStart_s_ - floor(dataI.TimeFromStart_s_);
    data.Acc1(1:b,1) = dataI.AccelIntX; data.Acc2(1:b,1) = dataI.AccelIntY; data.Acc3(1:b,1) = dataI.AccelIntZ;
    data.Comp1(1:b,1) = dataI.MagX; data.Comp2(1:b,1) = dataI.MagY; data.Comp3(1:b,1) = dataI.MagZ;
    data.Gyr1(1:b,1) = dataI.GyroX; data.Gyr2(1:b,1) = dataI.GyroY;  data.Gyr3(1:b,1) = dataI.GyroZ;
%     times = datenum(dataP.FileTime,'HH:MM:SS'); times = times-floor(times);
%     DNp = datenum(dataP.FileDate,'mm/dd/yy') + times + UTC/24;
    DNp = dates(2) + dataP.TimeFromStart_s_/24/60/60;
%     DN = UTCDN+UTC/24;
    pfs = 1/mean(diff(dataP.TimeFromStart_s_)); disp(['pfs = ' num2str(pfs)]);
    if fs/pfs ~= round(fs/pfs); error('sample rates do not divide evenly'); end
    % For some files with time steps that are inconsistent, we're assuming samples are evenly taken and the timestamp itself
    % does not matter.
    [a,I] = min(abs(data.Date+data.Time-DNp(1))); if a>5/24/60/60; error('some data is missing'); end
    if (I-1+fs/pfs*size(dataP,1))>size(data,1); 
%         while (I-1+fs/pfs*size(dataP,1))>size(data,1); dataP(end,:) = []; end
        error('Pressure extends past the end of the next file');
    end
    for i = 1:fs/pfs; data.Pressure(I+i-1:fs/pfs:I-1+fs/pfs*size(dataP,1)) = dataP.Pressure; data.Temp(i:fs/pfs:fs/pfs*size(dataP,1)) = dataP.Temperature; end
%     I = cell2mat(arrayfun(@(x) find(abs(DN-x)<.1/pfs/fs/24/60/60,1),DNp,'uniformoutput',false));
%     data.Pressure(I) = dataP.Pressure; data.Temp(I) = dataP.Temperature;
%     for i = 1:fs-1; data.Pressure(I+i) = dataP.Pressure; data.Temp(I+i) = dataP.Temperature; end
% %     for i = 1:length(I); try I2 = I(i+1)-1; catch; I2 = length(data.Date); end; data.Pressure(I(i):I2) = dataP.Pressure(i); data.Temp(I(i):I2) = dataP.Temperature(i); end
    for i = 3:2:length(D)
        %%
        dataI = readtable([fileloc D{i}],'readvariablenames',true);
        dataP = readtable([fileloc D{i+1}],'readvariablenames',true);
        %         times = datenum(dataI.FileTime,'HH:MM:SS'); times = times-floor(times);
        %         DN = datenum(dataI.FileDate,'mm/dd/yy')+times; UTC2 = (dates(i)-DN(1))*24; if abs(UTC2-round(UTC))>.1; error('check UTC offset'); end
        DN = dates(i) + dataI.TimeFromStart_s_/24/60/60-dataI.TimeFromStart_s_(1)/24/60/60;
        if i<length(D)-1; nextfile = dates(i+2); else nextfile = max(dates(i)+size(dataI,1)/fs/24/60/60,dates(i+1)+size(dataP,1)/pfs/24/60/60); end
        DN = (DN(1):1/fs/24/60/60:nextfile-1/fs/24/60/60)';
        if length(DN)<size(dataI,1);  %dataI = dataI(1:length(DN),:);
                    error('data extends beyond next file start');
        end
        fs2 = 1/mean(diff(dataI.TimeFromStart_s_)); if abs(fs-fs2)>.01 || abs(fs2-round(fs2))>.01; error('Missing data'); end
        %         times = datenum(dataP.FileTime,'HH:MM:SS'); times = times-floor(times);
        %         DNp = datenum(dataP.FileDate,'mm/dd/yy') + times + UTC/24;
        DNp = dates(i+1) + dataP.TimeFromStart_s_/24/60/60-dataP.TimeFromStart_s_(1)/24/60/60;
        pfs2 =  1/mean(diff(dataP.TimeFromStart_s_)); if abs(pfs-pfs2)>.01 || abs(pfs2-round(pfs2))>.01; error('Missing data'); end
%         prevhour = floor(DN(10)*24)/24;
        
%         lastend = data.Date(end)+data.Time(end);%-UTC/24+1/fs/24/60/60;
%         if abs(lastend-prevhour)>1/24/60/60; error('start/end not matching'); else [prevhour,lastend] = max([prevhour;lastend]); if lastend == 2; addpressure = true; else addpressure = false; end; end
%         nexthour = max(ceil(DN(end-10)*24)/24-1/fs/24/60/60,DNp(end)-UTC/24);
%         if nexthour <= DN(end); ex = sum(DN >= nexthour); if ex > fs; error ('checktiming'); end; nexthour = nexthour+(ex-1)/fs/24/60/60; end
%         UTCDN = (prevhour:1/fs/24/60/60:nexthour)';
%         [~,a] = min(abs(UTCDN-DN(1))); 
a = 1;
b = size(dataI,1); b = b+a-1;
%         if abs(UTCDN-DN(a))>5/24/60/60; error('First time point is more than 5 seconds off'); end
        UTCDN = DN - UTC/24;
        data1 = table (UTCDN,floor(DN),DN-floor(DN),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),...
            'VariableNames',{'UTCDN','Date','Time','Acc1','Acc2','Acc3','Comp1','Comp2','Comp3','Gyr1','Gyr2','Gyr3','Pressure','Temp'});
%         data1.Date = floor(DN); data1.Time = DN-floor(DN); %DN+UTC/24 - data.Date + dataI.dataI.TimeFromStart_s_ - floor(dataI.TimeFromStart_s_);
        data1.Acc1(a:b,1) = dataI.AccelIntX; data1.Acc2(a:b,1) = dataI.AccelIntY; data1.Acc3(a:b,1) = dataI.AccelIntZ;
        data1.Comp1(a:b,1) = dataI.MagX; data1.Comp2(a:b,1) = dataI.MagY; data1.Comp3(a:b,1) = dataI.MagZ;
        data1.Gyr1(a:b,1) = dataI.GyroX; data1.Gyr2(a:b,1) = dataI.GyroY;  data1.Gyr3(a:b,1) = dataI.GyroZ;
        %         I = cell2mat(arrayfun(@(x) find(abs(DN-x)<.1/fs/24/60/60,1),DNp,'uniformoutput',false));
        %         data1.Pressure(I) = dataP.Pressure; data1.Temp(I) = dataP.Temperature;
        [a,I] = min(abs(data1.Date+data1.Time-DNp(1))); if a>5/24/60/60; error('some data is missing'); end
        if (I-1+fs/pfs*size(dataP,1))>size(data1,1); 
%             while (I-1+fs/pfs*size(dataP,1))>size(data1,1); dataP(end,:) = []; end
            error('Pressure extends past the end of the next file'); 
        end
        for ii = 1:fs/pfs; data1.Pressure(I+ii-1:fs/pfs:I-1+fs/pfs*size(dataP,1)) = dataP.Pressure; data1.Temp(ii:fs/pfs:fs/pfs*size(dataP,1)) = dataP.Temperature; end
%         I2 = I(end); %I = I(1:end-1);
%         for ii = 1:fs-1; I(end) = I(end)-1; data1.Pressure(I+ii) = dataP.Pressure; data1.Temp(I+ii) = dataP.Temperature; end
%         for ii = I2:size(data1.Pressure,1); data1.Pressure(ii) = data.Pressure(end); end
%         if addpressure; data1.Pressure(1:find(DNp(1)>=DN,1)-1) = data.Pressure(end); end
        data = [data; data1];
    end
        
    oi = regexp(fileloc,'\');
    newfileloc = fileloc(1:oi(end-1));
    Adata = [data.Acc1 data.Acc2 data.Acc3]; Atime = data.Date+data.Time;
    warning('on','all');
    lastwarn('');
    try
        save([newfileloc D{1}(1:end-3) 'mat'],'data','Adata','Atime');
%         if ~isempty(notes); save([newfileloc file(1:end-3) 'mat'],'notes','-append'); end
        if ~isempty(lastwarn)
            error(lastwarn);
        end
    catch %v7.3 allows for bigger files, but makes a freaking huge file if used when you don't need it
       save([newfileloc D{1}(1:end-3) 'mat'],'data','Adata','Atime','-v7.3');
%         else save([newfileloc file(1:end-3) 'mat'],'data','Adata','Atime','notes','-v7.3'); end
        disp('Made a version 7.3 file in order to include all');
    end
if exist('pconst','var'); p = (data.Pressure-pconst)*pcal; else p = data.Pressure; end
if any(data.Pressure>10)&&min(data.Pressure)<.5*max(data.Pressure) % plot something if it went underwater
    if ~exist('fs','var'); fs = round(1/((data.Time(50)-data.Time(49))*60*60*24)); end
    oi = max(1, find(p>3,1,'first')-60*fs):min(find(p>4,1,'last')+60*fs,length(p));
    if ~isempty(oi)
        oi1 = oi(1:floor(length(oi)/3)); oi2 = oi(floor(length(oi)/3):2*floor(length(oi)/3)); oi3 = oi(2*floor(length(oi)/3):end);
        %         data.Date = floor(DN); data.Time = DN-floor(DN);
        time = data.Date + data.Time;
        figure(1); clf; set(1, 'units','normalized','outerposition',[0 0 1 1]);
        subplot(3,1,1); plot(data.Date(oi1)+data.Time(oi1),p(oi1)); ylim([-5 max(p(oi1))]); xlim(time([oi1(1) oi1(end)])); set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM'),'ydir','rev');
        subplot(3,1,2); plot(data.Date(oi2)+data.Time(oi2),p(oi2)); ylim([-5 max(p(oi2))]); xlim(time([oi2(1) oi2(end)])); set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM'),'ydir','rev');
        subplot(3,1,3); plot(data.Date(oi3)+data.Time(oi3),p(oi3)); ylim([-5 max(p(oi3))]); xlim(time([oi3(1) oi3(end)])); set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM'),'ydir','rev');
        %         print(1,[fileloc 'TDR.jpg'],'-djpeg');
        saveas(1,[newfileloc 'TDR3.bmp']);
        saveas(1,[newfileloc 'TDR3.fig']);
    end
end
disp(newfileloc);
disp(['data on: ' datestr(data.Date(1)+data.Time(1),'mm/dd/yy HH:MM:SS') ' data off: ' datestr(data.Date(end)+data.Time(end),'mm/dd/yy HH:MM:SS')]);
figure(2); clf;
s1=subplot(311);
plot([data.Acc1 data.Acc2 data.Acc3]); title('Accelerometer');
s2=subplot(312);
plot([data.Comp1 data.Comp2 data.Comp3]); title('Magnetometer');
s3 = subplot(313);
plot(data.Pressure); title('Pressure');
linkaxes([s1 s2 s3],'x');