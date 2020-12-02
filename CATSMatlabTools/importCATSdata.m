function [data, Adata, Atime] = importCATSdata(fileloc, filename,FS,importAll)
%
% David Cade
% version 11.23.2020
% Goldbogen Lab
% Stanford University


% Matlab packages required: Signal Processing Toolbox
%
% Reads tag data and outputs a matlab formatted data table "data", as well as
% "Adata" (the accelerometer data at the original sample rate) and "Atime" 
% a time stamp for each Adata point.  "data" is sampled at the highest 
% non-accelerometer sample rate recorded, and Acc data in "data" is downsampled
% to that rate.  

% Prerequisites: 
% 1) Data downloaded from tag in split csv format (recommended ~100 MB, but
% can handle smaller or larger if there is sufficient RAM).  IF AND ONLY IF
% DATA IS IN A SINGLE CSV > 200 MB, then:
    % use a csv splitter (e.g.
    % https://download.cnet.com/CSV-Splitter/3000-2074_4-75910188.html) to
    % split it.  Include the header row in each new csv file, and put the
    % output files in a "csvs" folder within the same folder as the large csv
    % file.
% 2) Ensure that the ".txt" file with sampling rate information is in the
% same folder as the downloaded csv file (or files).
% 3) folder structure.  There is some flexibility here, but ease of use
% will be facilitated if all data directly from the tag (csvs, ubx, txt,
% bin, and videos) are in a single folder labeled "raw" within whatever
% organizing folder (typically the deployment ID).

% script format:
% importCATSdata(); % prompts to select file
% importCATSdata(fileloc, filename,FS);
% importCATSdata(fileloc, filename);
% importCATSdata(FS);
% importCATSdata([],[],[],true); % sets the importAll parameter to
% true.  setting this to true ensures that all csvs are read.  If false,
% data stops once pressure remains at the surface for an entire csv read.
% [data, Adata, Atime] = importCATSdata(...);
%
%set FS if you want to set the frequency of data table (you will have to
%ensure that the values divide evenly).  Default is to use the maximum non
%acceleration value.  FS = new frequency
%
%filename is the original csv (if large) or the first csv if files are
%split from the tag.  Can also select a later file if you know that the
%first xx files are from before the deployment.


%% Section 1, set up files to run through and import
% nargin = 0; %uncomment this line if you are running from the code (not as a function)
deletecsvs = true; %if you want to delete csv file after it is created from the split file and then read
if nargin <2 || isempty(fileloc) || isempty(filename)
    if nargin == 1; FS = fileloc; end
    [filename,fileloc] = uigetfile('*.csv','select original csv file, if it''s big, split into parts in a ''csvs'' folder');
end
if nargin<4; importAll = false;
end
% 
if iscell(filename); fname = filename{1}; else fname = filename; end
%
% looks for a 'csv' folder with the split csv files.
fileloccsv = [fileloc 'csvs\']; simple = false; %simple indicates whether the file was huge and you broke it up into pieces
if ~exist(fileloccsv,'dir'); simple = true; fileloccsv = fileloc;
else
    FILES = dir(fileloccsv); FILES ={FILES.name};
    [~,b] = max(cellfun(@(x) sum(x(1:min(length(x),length(fname))) == fname(1:min(length(x),length(fname)))),FILES));
    fname = char(FILES(b)); fname = [fname(1:end-3) 'csv'];%ensures that if you chose a text file it has the same format as the csv file (sometimes the time is slightly different
end
DIR = dir(fileloccsv);

disp('Can watch the csvs folder with partial files to gauge progress (the currently uploading file is a copy).');
disp('When processing is complete, if not all files were read a copy of the last csv read will be left in the csvs folder.');

% allows you to start the import at a file number above 0
if regexp(filename(end-7:end-4),'_\d\d\d')
    i = strfind(filename,'_');
    i = str2num(filename(i+1:i+3));
else
    i = 0;
end
mini = i;
notes = [];
Adata = nan(0,3); rownum = 1; Atime = nan(0,1);
noPress = false;
extraData = [];

%% Section 2 runs through the csv files and reads them into the data frame
while any(strcmp({DIR.name},[fname(1:end-3) num2str(i,'%03u')])) || any(strcmp({DIR.name},[fname(1:end-7) num2str(i,'%03u') '.csv'])) || any(strcmp({DIR.name},[fname(1:end-4) '_' num2str(i,'%03u') '.csv'])) || simple
    % this reads in a variety of file names.  To add new file name types, add in a new try catch format
    try file = DIR(strcmp({DIR.name},[fname(1:end-3) num2str(i,'%03u')])).name; froot = fname(1:end-4);
    catch; try file = DIR(strcmp({DIR.name},[fname(1:end-7) num2str(i,'%03u') '.csv'])).name; froot = fname(1:end-8);
        catch; try file = DIR(strcmp({DIR.name},[fname(1:end-4) '_' num2str(i,'%03u') '.csv'])).name; froot = fname(1:end-4);
            catch; file = fname; froot = fname(1:end-4); end;
        end;
    end
    try copyfile([fileloccsv file],[fileloccsv file '.csv']); % makes a copy of the file
        % this checks to see if the user has deleted any files since the
        % process started (e.g. if it is noted that all files after a
        % certain point are not relevant)
    catch er
        if ~exist([fileloccsv file],'file');
            warning(['csvs greater than ' num2str(i-1) ' not found.  Saving as is.']);
            break
        end
        throw (er)
    end
    % This is a check.  If the data is corrupted , this throws us out of
    % the loop to save the data we already have.
    try if baddataEnd; break; end
    catch
    end
    %This stops the import if at some point there are at least 100 samples
    %deeper than 10 m (i.e. the file was on an animal) and the current file
    %has no change in pressure greater than 1 m (e.g. the tag is floating
    %at the surface)
    if ~importAll
        try if max(dataT.Pressure)<mean(dataT.Pressure)+1 && sum(data.Pressure>10)>100
                disp(['Tag appears off, check profile but saving as is (through csv ' num2str(i-1) '). If tag is not yet off, rerun importCATSdata with 4th argument set to true.']);
                break
            end
        catch
        end
    end
    % if this is the first file
    % most of the "try"s are to allow older file versions to work if they
    % do not have the newer data columns.  If your file has additional data
    % columns that should be imported, may need to add them in here.
    if i == mini
        % check the version of matlab and how readtable might read data.
        vers = version('-release'); if strcmp(vers(end),'a'); vers = str2num(vers(1:4)); else vers = str2num(vers(1:4))+0.1; end
       if vers<2020
        data = readtable([fileloccsv file '.csv'],'headerlines',0,'readvariablenames',false);
       else 
           data = readtable([fileloccsv file '.csv'],'headerlines',0,'readvariablenames',false,'Format','auto');
       end
            headers = data{1,:};
        data(1,:) = [];
         
        try headers(~cellfun(@isempty,strfind(headers,'Light'))) = {'Light1' 'Light2'}; catch; headers(~cellfun(@isempty,strfind(headers,'Light'))) = {'Light'}; end
        % you may wish to be more descriptive with your "Temp" heading.
        try headers(~cellfun(@isempty,strfind(headers,'Pressure'))) = {'Pressure' 'Temp'}; catch;...
                try headers(~cellfun(@isempty,strfind(headers,'Depth'))) = {'Pressure' 'Temp'}; catch; ...
                    try headers(~cellfun(@isempty,strfind(headers,'Depth'))) = {'Pressure'};  headers(~cellfun(@isempty,strfind(headers,'Temperature (depth'))) = {'Temp'};  catch; ...
                        headers(~cellfun(@isempty,strfind(headers,'Pressure'))) = {'Pressure'}; end;
                end;
        end
        if sum(cellfun(@(x) strcmp(x,'Pressure'),headers)) == 0;
           warning('No Pressure data detected, continue? 1 = yes, 2 = no');
           x = input('?');
           if x ~=1;
               error ('Fix pressure data');
           else noPress = true;
           end            
        end
        %This adds a numeral to every subsequently found version of
        %Temperature that is in the csv
        if sum(strcmp(headers,'Temp')) == 0; headers(~cellfun(@isempty,strfind(headers,'Temp'))) = {'Temp'}; else
            ii = 1; while sum(~cellfun(@isempty,strfind(headers,'Temperature'))) > 0 || sum(~cellfun(@isempty,strfind(headers,'Temp. (mag'))) > 0; headers(max(find(~cellfun(@isempty,strfind(headers,'Temp'))&~cellfun(@(x) strcmp(x,'Temp'),headers),ii,'first'))) = {['Temp' num2str(ii)]}; ii = ii+1; end
        end
        % Magnetometers are described as "Comp" for compass
        try headers(~cellfun(@isempty,strfind(headers,'Comp'))) = {'Comp1' 'Comp2' 'Comp3'};  catch; headers(~cellfun(@isempty,strfind(headers,'Magnet'))) = {'Comp1' 'Comp2' 'Comp3'}; end
        headers(~cellfun(@isempty,strfind(headers,'Gyr'))) = {'Gyr1' 'Gyr2' 'Gyr3'};
        headers(~cellfun(@isempty,strfind(headers,'Acc'))) = {'Acc1' 'Acc2' 'Acc3'};
        try headers(~cellfun(@isempty,strfind(headers,'Speed'))) = {'Speed'}; catch; end
        delcol = [];
        if sum(~cellfun(@isempty,strfind(headers,'Date'))) == 2; % if there is both local and UTC time.
            delcol(1) = find(~cellfun(@isempty,strfind(headers,'Date')),1);
            delcol(2) = find(~cellfun(@isempty,strfind(headers,'Time')),1);
        end
        % delete columns you don't want to import
        try headers(~cellfun(@isempty,strfind(headers,'System error'))) = {'SystemError'};
            delcol = [delcol find(~cellfun(@isempty,strfind(headers,'SystemError')),1)];
        catch; end
        try delcol = [delcol find(~cellfun(@isempty,strfind(headers,'Flags')),1)];
        catch; end
        % GPS columns with more extensive headers are imported below. There
        % are some versions with just "GPS" that did not have useful data
        try delcol = [delcol find(strcmp(headers,'GPS'))];
        catch; end
        try delcol = [delcol find(~cellfun(@isempty,strfind(headers,'CC vid. size')),1)];
        catch; end
        data(:,delcol) = []; headers(delcol) = [];
        if ~isempty(strfind(headers{~cellfun(@isempty,strfind(headers,'Date'))},'UTC')); UTCflag = true; else UTCflag = false; end
        headers(~cellfun(@isempty,strfind(headers,'Date'))) = {'Date'};
        headers(~cellfun(@isempty,strfind(headers,'Time'))) = {'Time'};
        try headers(~cellfun(@isempty,strfind(headers,'GPS'))) = {'GPSDate' 'GPSTime' 'GPSsat1' 'GPSsat2'}; catch; try headers(~cellfun(@isempty,strfind(headers,'GPS'))) = {'GPSDate' 'GPSTime' 'GPSsat'}; catch; end; end
        try headers(~cellfun(@isempty,strfind(headers,'BATT'))) = {'BATTv' 'BATTmA' 'BATTmAh'}; catch;  try headers(~cellfun(@isempty,strfind(headers,'BATT'))) = {'BATTv' 'BATTmA'}; catch; try headers(~cellfun(@isempty,strfind(headers,'BATT'))) = {'BATTv'}; catch; end;  end; end
        try headers(~cellfun(@isempty,strfind(headers,'Camera time'))) = {'CamTime'}; catch; end
        try headers(~cellfun(@isempty,strfind(headers,'CC status'))) = {'CamOn'}; catch; end
        try headers(~cellfun(@isempty,strfind(headers,'CC vid. '))) = {'VidSize'}; catch; end
        data.Properties.VariableNames = headers;
        
        [accHz,gyrHz,magHz,pHz,lHz,GPSHz,UTC,THz,T1Hz] = sampledRates(fileloc,file);
        if ~UTCflag; UTC = 0; end
        if exist('FS','var') && ~isempty(FS); fs = FS; %if you preset the maxFS, else use the max of the others
        else fs = max([gyrHz magHz pHz lHz GPSHz]);
        end
        data = data(find(~strcmp(data.Acc1,'0'),1,'first'):end,:);
        dataT = data;
    else % if not on the first csv imported
        if vers<2020
        dataT = readtable([fileloccsv file '.csv'],'headerlines',0,'readvariablenames',false);
        else 
           dataT = readtable([fileloccsv file '.csv'],'headerlines',0,'readvariablenames',false,'Format','auto');
       end 
      
        dataT(:,delcol) = [];
        if any(cellfun(@any,cellfun(@(x) strfind(x,'Acc'),table2cell(dataT(1,:)),'uniformoutput',false))) % if the csv has headers
            dataT(1,:) = []; % this line gets rid of any copied headers
        end
        if noPress; dataT.Pressure = zeros(size(dataT,1),1); end
        dataT.Properties.VariableNames = headers;
        DT = datenum(dataT.Date,'dd.mm.yyyy');
        if any(abs(diff(DT))>5);
            stupidi = find(abs(diff(DT))>1,1);
            dataT = dataT(1:stupidi,:);
            disp(['There is bad data of some kind starting at line ' num2str(stupidi + 1) ' in csv ' num2str(i-1) '. Stopping import at this point.']);
            baddataEnd = true;
        else baddataEnd = false;
        end
    end
    if isempty(extraData); extraData = dataT(1,:); extraData(1,:) = []; end
    dataT = [extraData; dataT]; % add in any leftover rows from last reading.
    
    % this checks for errors in a specific row.  If you get this error
    % check your csv.  If it is just a couple of bad rows, this smooths
    % them out with interpolated values.  If possible, try to redownload
    % the csvs.
    oiT = char(dataT.Time);
    goodms = true;
    badrows = [];
    for iii = 1:size(dataT,2)
        if ~any(strcmp(headers{iii},{'Time','GPSTime','GPSDate'})) && iscell(dataT{1,iii}(1))
            if any(cellfun(@isempty, dataT{:,iii}) | ~cellfun(@isempty, regexp(dataT{:,iii},'[x,:]')))
                badrow = find(cellfun(@isempty, dataT{:,iii}) | ~cellfun(@isempty, regexp(dataT{:,iii},'[x,:]')));
                badrows = [badrows; badrow];
                goodms = false;
            end
        end
    end
    badrows = unique(badrows); badrows = sort(badrows);
    for iii = 1:length(badrows)
        disp(['Bad row in csv ' num2str(i) ', row: ' num2str(badrows(iii)) ', row replaced with mean accelerometer of immediately surrounding rows, and other values replaced with subsequent row.']);
    end
    
    % this is for an older version of the csv files in which the ms part of
    % the time stamps were messed up (they had an extra digit).  should not
    % be relevant for most versions, though having the check does not
    % inhibit the import of newer data
    if ~goodms
        for iii = 1:length(badrows)
            ms = str2num(oiT(badrows(iii)-1,10:end));
            ms = ms/10^ceil(log10(max(ms))); %thi
            d = datestr( datenum(oiT(badrows(iii)-1,1:8))+(ms+1/accHz)/24/60/60,'HH:MM:SS');
            %             dms = diff(str2num(oiT(badrows(iii)-2:badrows(iii)-1,10:end)))
            oiT(badrows(iii),:) = [d '.' oiT(badrows(iii)-1,10:end)];
            %             timescal(badrows(iii)) = datenum(oiT(badrows(iii),1:8),'HH:MM:SS');
            %             if ms == 0;
            dataT.Date(badrows(iii)) = dataT.Date(badrows(iii)-1);
            %             else dataT.Date(badrows(iii)) = dataT.Date(badrows(iii)+1);
            %             end
            for iv = 3:size(dataT,2);
                dataT{badrows(iii),iv} = {'NaN'};
            end
            dec = num2str((length(char(dataT.Acc1(badrows(iii)-1)))-strfind(char(dataT.Acc1(badrows(iii)-1)),'.')));
            dataT.Acc1(badrows(iii)) = {sprintf(['%0.' dec 'f'], nanmean(str2num(char(dataT.Acc1(badrows(iii)-1:badrows(iii)+1)))))};
            dataT.Acc3(badrows(iii)) = {sprintf(['%0.' dec 'f'], nanmean(str2num(char(dataT.Acc3(badrows(iii)-1:badrows(iii)+1)))))};
            dataT.Acc2(badrows(iii)) = {sprintf(['%0.' dec 'f'], nanmean(str2num(char(dataT.Acc2(badrows(iii)-1:badrows(iii)+1)))))};
        end
    end
    
    timescal = datenum(oiT(:,1:8));
    
    DN = timescal-floor(timescal)+datenum(dataT.Date,'dd.mm.yyyy');
    isb = arrayfun(@(x) strcmp(x,' '),oiT(:,end)); % some bad imports skipped a few seconds of data by going from 10.9 to 10.10 seconds (e.g.)
    if any(isb)
        bad10 = find(~isb);
        DN(bad10) = DN(bad10)+1/24/60/60; % add a second back
        oiT(bad10,10:end-1) = oiT(bad10,11:end); %move the values forward
        oiT(:,end) = []; % get rid of the extra column
        warning(['SOME TIMESTAMPS STARTING ROW ' num2str(bad10(1)) ' AND ENDING ROW ' num2str(bad10(end)) ' IN CSV ' num2str(i) ' DON''T COUNT ACCURATELY.  SHOULD BE FIXED BUT HIGHLY RECOMMEND REDOWNLOADING DATA TO ENSURE NOTHING WAS MISSED OR FIXED INACCURATELY.']);
    end
    
    % similarly, this fixes old versions with poorly formatted timestamps.
    ms = str2num(oiT(:,10:end));
    ms = ms/10^ceil(log10(max(ms))); %this accounts for some of the funny cats formatting (extra zeroes etc.)
    for iii = 1:length(badrows)
        ms(badrows(iii)) = ms(badrows(iii)-1)+1/accHz;
    end
    DN = DN + ms/24/60/60+UTC/24; % brings back to local time
    for iii = 1:length(badrows)
        if isempty(notes); notes = 'Bad data at t = '; end
        notes = [notes datestr(DN(badrows(iii)),'mm/dd HH:MM:SS.fff') ', '];
    end
    time = DN-floor(DN);
    DN(badrows(time(badrows)==0)) = DN(badrows(time(badrows)==0))+1; % accounts for any points when time switched to the next day that happened to be exactly at a bad row
    
    if iscell(dataT.Acc1);  oi = [str2num(char(dataT.Acc1)) str2num(char(dataT.Acc2)) str2num(char(dataT.Acc3))];
    else oi = [dataT.Acc1 dataT.Acc2 dataT.Acc3]; end
    oi2 = fixgaps(oi);
    oi(badrows,:) = oi2(badrows,:); clear oi2;
    d = diff(DN*24*60*60);
    skippeddata = find(d>1.5*1/accHz); %
    numgaps = 0;
    while ~isempty(skippeddata)
        b = skippeddata(1);
        oi = [oi(1:b,:); nan(1,3); oi(b+1:end,:)];
        oi2 = fixgaps(oi);
        oi(b+1,:) = oi2(b+1,:); clear oi2;
        DN = [DN(1:b); DN(b)+1/accHz/24/60/60; DN(b+1:end)];
        dataT = [dataT(1:b,:); dataT(b+1,:); dataT(b+1:end,:)];
        badrows(badrows>b) = badrows(badrows>b)+1;
        numgaps=numgaps+1;
        if numgaps>length(badrows)*5;
            warning([num2str(numgaps) ' gaps in data fixed.  Stopping since more than 5x number of bad rows fixed.  Check data for inaccuracies']);
            break;
        end
        d = diff(DN*24*60*60);
        skippeddata = find(d>1.5*1/accHz); clear d;
    end
    if numgaps>0; disp([num2str(numgaps) ' additional gaps in data filled']); end
    
    % this line ensures all data decimates to the right length, if there
    % are extra points, carry them to the next csv;
    if mod(size(dataT,1),accHz/fs) ~=0
        lastrows = size(dataT,1)-mod(size(dataT,1),accHz/fs)+1:size(dataT,1);
        extraData = dataT(lastrows,:);
        dataT(lastrows,:) = [];
        DN(lastrows,:) = [];
        oi(lastrows,:) = [];
    else
        extraData = [];
    end
    dataT.Time = DN-floor(DN);
    dataT.Date = floor(DN);
    
    Atime = [Atime; DN];
    Adata = [Adata; oi];
    
    % downsample data
    acc = decdc(oi,accHz/fs); %used to round these values, but doesn't work if they are pre-calibrated.
    % for other, non acc data, sample it (data are repeated values);
    % takes into account different kinds of imported data (numbers or
    % strings)
    % FS is set at the beginning or is set to be the highest sampled data
    % that's not the accelerometer (often 50 Hz).
    if exist('FS','var') && ~isempty(FS)
        if ~iscell(dataT.Gyr1); gyr = [dataT.Gyr1 dataT.Gyr2 dataT.Gyr3]; else gyr = [str2num(char(dataT.Gyr1)) str2num(char(dataT.Gyr2)) str2num(char(dataT.Gyr3))]; end
        gyr = gyr(1:accHz/fs:end,:);
        if ~iscell(dataT.Comp1); comp = [dataT.Comp1 dataT.Comp2 dataT.Comp3]; else comp = [str2num(char(dataT.Comp1)) str2num(char(dataT.Comp2)) str2num(char(dataT.Comp3))]; end
        comp = comp(1:accHz/fs:end,:);
    end
    if noPress && ~any(cellfun(@(x) strcmp(x,'Pressure'),headers)); dataT.Pressure = zeros(size(dataT.Comp1)); headers{end+1} = 'Pressure'; end
    if ~iscell(dataT.Pressure); p = dataT.Pressure; else p = str2num(char(dataT.Pressure)); end
    p = p(1:accHz/fs:end);
    samples = abs(round(dataT.Time*24*60*60*fs)-dataT.Time*24*60*60*fs)<1/(accHz*1.5); %find the values closest to the rounded frequencies and sample those.
    for ii = length(badrows):-1:1
        if samples(badrows(ii))
            nextgoodrow = badrows(ii)+1;
            while ismember(nextgoodrow,badrows)
                nextgoodrow = nextgoodrow+1;
            end
            dataT(badrows(ii),3:end) = dataT(nextgoodrow,3:end);
        end
    end
    
    dataT = dataT(samples,:);
    
    try
        dataT.CamOn = ~cellfun(@isempty, cellfun(@(x) strfind(x,'R'),dataT.CamOn,'uniformoutput',false));
    catch
    end
    
    for j = 1:size(dataT,2) %assumes the first two columns are the dates
        if ~any(strcmp(headers{j},{'Date','Time','GPSDate','GPSTime','GPSsat','CamOn'})) && iscell(dataT{1,j}(1))
            dataT.(headers{j}) = str2num(char(dataT.(headers{j})));
        end
    end
%     %now accounted for by carrying over data to the next csv
%     % sometimes decdc is one short from the sampled points, this tries to
%     % fit it in as best as possible and then fill in the gap, also only
%     % resamples p, comp, gyr if you set FS to be less than the maximum of
%     % those
%     I = size(dataT,1)-size(acc,1)+1;
%     I2 = sum(abs(dataT.Pressure(I:end)-p));
%     I1 = sum(abs(dataT.Pressure(1:end-I+1)-p));
%     adjcols = any([ ~cellfun(@isempty,strfind(headers,'Comp'));  ~cellfun(@isempty,strfind(headers,'Acc'));  ~cellfun(@isempty,strfind(headers,'Gyr'));  ~cellfun(@isempty,strfind(headers,'Pressure'))]);
%     if size(dataT,1) == size(acc,1) || I2<=I1
%         dataT.Acc1(I:end) = acc(:,1); dataT.Acc2(I:end) = acc(:,2); dataT.Acc3(I:end) = acc(:,3);
%         if exist('FS','var') && ~isempty(FS)
%             dataT.Comp1(I:end) = comp(:,1); dataT.Comp2(I:end) = comp(:,2); dataT.Comp3(I:end) = comp(:,3);
%             dataT.Gyr1(I:end) = gyr(:,1); dataT.Gyr2(I:end) = gyr(:,2); dataT.Gyr3(I:end) = gyr(:,3);
%             dataT.Pressure(I:end) = p;
%             % if there was a leftover point
%             if I2<I1 ; for jj = 1:length(adjcols); if adjcols(jj)&&i~=mini; dataT{1,jj} = nan; end; end; end
%         end
%     else
%         dataT.Acc1(I:end) = acc(:,1); dataT.Acc2(I:end) = acc(:,2); dataT.Acc3(I:end) = acc(:,3);
%         if exist('FS','var') && ~isempty(FS)
%             dataT.Comp1(1:end-I+1) = comp(:,1); dataT.Comp2(1:end-I+1) = comp(:,2); dataT.Comp3(1:end-I+1) = comp(:,3);
%             dataT.Gyr1(1:end-I+1) = gyr(:,1); dataT.Gyr2(1:end-I+1) = gyr(:,2); dataT.Gyr3(1:end-I+1) = gyr(:,3);
%             dataT.Pressure(1:end-I+1) = p;
%         end
%         for jj = 1:length(adjcols); if adjcols(jj)&&i~=mini; dataT{end,jj} = nan; end; end;
%     end
%     testcol = find(adjcols,1,'first');
%     if i~=mini && (isnan(dataT{1,testcol}) || isnan(data{end,testcol}))
%         for jj = 1:length(adjcols); if adjcols(jj)
%                 col = inpaint_nans([data{end-10:end, jj};dataT{1:11,jj}]);
%                 data{end-10:end,jj} = col(1:11);
%                 dataT{1:11,jj} = col(12:end);
%             end; end
%     end
    
    if ~isnan(GPSHz) && iscell(dataT.GPSDate(1))
        I = find(cellfun(@isempty,dataT.GPSDate),1,'last'); if isempty(I); I = 0; end
        if I == size(dataT,1); dataT.GPSDate = nan(size(dataT.GPSDate)); dataT.GPSTime = nan(size(dataT.GPSTime));
        else
            try oiD = [nan(I,1); datenum(dataT.GPSDate(I+1:end),'dd.mm.yyyy')];
                oi = char(dataT.GPSTime(I+1:end));
                timescal = datenum(oi(:,1:11)); %changed to try to capture fractional seconds
                dataT.GPSTime = [nan(I,1); timescal - floor(timescal)];
            catch; oiD = [nan(I,1); str2num(char(dataT.GPSDate(I+1:end)))];
                dataT.GPSTime = [nan(I,1); str2num(char(dataT.GPSTime(I+1:end)))];
            end
            dataT.GPSDate = oiD;
        end
    end
    
    DN = dataT.Date+dataT.Time;
    d = diff(DN*24*60*60);
    skippeddata = find(d>1.5*1/fs);
    if ~isempty(badrows)
        if ~isempty(skippeddata)
            disp(['Results in ' num2str(length(skippeddata)) ' gaps in downsampled data.  Should be fixed during prh process.']);
        else
            disp('No gaps found in downsampled data');
        end
    end
    
    if i == mini; data = dataT; else
        data = [data;dataT];
    end
    i = i+1;
    if deletecsvs; delete([fileloccsv file '.csv']); end
    simple = false;
    
end
try if strcmp(fileloc(end-3:end-1),'raw'); newfileloc = fileloc(1:end-4); else newfileloc = fileloc; end; catch; newfileloc = fileloc; end
lastwarn('');
if regexp(file(end-7:end-4),'_\d\d\d')
    file = [file(1:end-8) file(end-3:end)];
end
if ~exist('THz','var'); THz = nan; end; if ~exist('T1Hz','var'); T1Hz = nan; end
Hzs = struct('accHz',accHz,'gyrHz',gyrHz,'magHz',magHz,'pHz',pHz,'lHz',lHz,'GPSHz',GPSHz,'UTC',UTC,'THz',THz,'T1Hz',T1Hz);
try
    save([newfileloc file(1:end-3) 'mat'],'data','Adata','Atime','Hzs');
    if ~isempty(notes); save([newfileloc file(1:end-3) 'mat'],'notes','-append'); end
    if ~isempty(lastwarn)
        error(lastwarn);
    end
catch %v7.3 allows for bigger files, but makes a freaking huge file if used when you don't need it
    if isempty (notes); save([newfileloc file(1:end-3) 'mat'],'data','Adata','Atime','Hzs','-v7.3');
    else save([newfileloc file(1:end-3) 'mat'],'data','Adata','Atime','notes','Hzs','-v7.3'); end
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
plot(data.Date+data.Time,[data.Acc1 data.Acc2 data.Acc3]); title('Accelerometer');
s2=subplot(312);
plot(data.Date+data.Time,[data.Comp1 data.Comp2 data.Comp3]); title('Magnetometer');
s3 = subplot(313);
plot(data.Date+data.Time,data.Pressure); title('Pressure'); set(s3,'ydir','rev');
linkaxes([s1 s2 s3],'x');
set([s1 s2 s3],'xticklabel',datestr(get(s1,'xtick'),'HH:MM:SS'))
disp('For figure 2, use "set([s1 s2 s3],''xticklabel'',datestr(get(s1,''xtick''),''HH:MM:SS''))" to reset axis labels if you zoom in');
if any(sum(diff([data.Acc1 data.Acc2 data.Acc3])==0)>.25*size(data,1)); error('CHECK ACCELEROMETER GRAPH, may have dropped an axis'); end
Mdiff = diff([data.Comp1 data.Comp2 data.Comp3]);
Mdiff2 = [[data.Comp1(1:end-2) data.Comp2(1:end-2) data.Comp3(1:end-2)] - [data.Comp1(3:end) data.Comp2(3:end) data.Comp3(3:end)]; [0 0 0]];
Mdiff3 = [[data.Comp1(1:end-3) data.Comp2(1:end-3) data.Comp3(1:end-3)] - [data.Comp1(4:end) data.Comp2(4:end) data.Comp3(4:end)]; [0 0 0; 0 0 0]];
Mdiff4 = [[data.Comp1(1:end-4) data.Comp2(1:end-4) data.Comp3(1:end-4)] - [data.Comp1(5:end) data.Comp2(5:end) data.Comp3(5:end)]; [0 0 0; 0 0 0; 0 0 0]];
Mdiff5 = [[data.Comp1(1:end-5) data.Comp2(1:end-5) data.Comp3(1:end-5)] - [data.Comp1(6:end) data.Comp2(6:end) data.Comp3(6:end)]; [0 0 0; 0 0 0; 0 0 0; 0 0 0]];
if any(sum(Mdiff == 0 & Mdiff2 == 0 & Mdiff3 == 0 & Mdiff4 == 0 & Mdiff5 == 0)>.25*size(data,1)); error('CHECK MAGNETOMETER GRAPH, may have dropped an axis'); end

try
    if sum(data.Pressure>5)<fs*120; error('less than two minutes of dives below 5 m'); end
    tagon = gettagon(data.Pressure,fs,data.Date(1)+data.Time(1),[data.Acc1 data.Acc2 data.Acc3]);
    save([newfileloc file(1:end-3) 'mat'],'tagon','-append');
catch
    disp('No pressure data (or less than two minutes of dives below 5 m) or error in gettagon function, could not determine tag on time')
end
