function [data, Adata, Atime] = importTDR10toCATS(fileloc, filename)
%
% David Cade and James Fahlbusch
% version 11.27.2020
% Goldbogen Lab
% Stanford University
%
% This script imports data from the Wildlife Computers TDR10F to the new 
% CATS data processing standards (11/2020), which includes the addition of 
% Adata, Atime and Hzs to the output of the script
%
% format:
% importTDR10toCATS(fileloc, filename);
% importTDR10toCATS();
% [data, Adata, Atime] = importTDR10toCATS(...);
%
% filename is the Archive CSV file outputted by DAP Processor
%
% NOTE: Raw data from TDR10 is in GMT; No conversion to Local time takes place in this script

% nargin = 0; %uncomment this line if you are running from the code (not as a function)
if nargin <2;
    [filename,fileloc] = uigetfile('*.csv','Select the Wildlife Computers Archive CSV file (exported from DAP Processor).');
end

%% 1. Import Section
        opts = detectImportOptions([fileloc filename]);
        opts.ExtraColumnsRule = 'ignore';
        data = readtable([fileloc filename],opts);
        headers = opts.VariableNames;
        %Change the uncorrected depth field name to pressure
        headers{strcmp(headers,'Depth')} = 'Pressure';
        %Change CorrectedDepth field name to Depth
        % headers{~cellfun(@isempty,strfind(headers,'CorrectedDepth'))} = 'Depth';
        %Change LightLevel field name to Light
        headers{strcmp(headers,'LightLevel')} = 'Light';
        try headers{~cellfun(@isempty,strfind(headers,'ExternalTemperature'))} = 'Temp'; catch; end
        try headers{~cellfun(@isempty,strfind(headers,'InternalTemperature'))} = 'TempInternal'; catch; end
        try headers{~cellfun(@isempty,strfind(headers,'DepthSensorTemperature'))} = 'TempDepthInternal'; catch; end
        headers(~cellfun(@isempty,strfind(headers,'intA'))) = {'Acc1' 'Acc2' 'Acc3'};
        %Uncomment to manually set the headers
%       headers = {'Date','Pressure','TempInternal','Temp','TempDepthSensor','LightLevel','Acc1','Acc2','Acc3','BatteryVoltage','Wet_Dry','Depth','Wet_DryThreshold','Events'};
        data.Properties.VariableNames = headers;
        disp('Section 1 (Import) finished');
%% 2. Data Processing
% This section creates a Matlab Datenumber variable and checks for gaps in
% the ACC data.
    % Events (e.g. 'Fast GPS attempt') are logged as separate lines in the data 
    % and must be removed. These are saved in a variable called 'events'
    events = [data(find(data.Events ~= ""),find(strcmp(data.Properties.VariableNames,'Time'))),data(find(data.Events ~= ""),find(strcmp(data.Properties.VariableNames,'Events')))];
    data = data(find(data.Events == ""),:);
    notes = 'Variables sampled at lower rates (e.g. Temp, Light, TempDepthInternal) are filled with repeated values to match the sampling rate of Pressure.';
    % Matlab has issues importing time that sometimes has milliseconds and
    % sometimes doesn't (i.e. at the top of every second, Wildlife Computers
    % omits the decimal seconds instead of .00000)
    % Need to process with milliseconds and without seperetely and combine
    dateMill = datetime(data.Time, 'inputformat', 'MM/dd/yyyy HH:mm:ss.SSSSS'); 
    dateFull = datetime(data.Time, 'inputformat', 'MM/dd/yyyy HH:mm:ss');
    dateNaT = isnat(dateMill);
    dateMill(dateNaT) = dateFull(dateNaT);
    %plot(diff(dateMill)) % not working
    datetimeUTC = dateMill;
    clearvars dateFull dateMill dateNaT 
    % This generates enormous files when saving data. 
    % data.datetimeUTC = dateMill;
    figure(1); clf;
    plot(datetimeUTC)
    disp('First 5 time stamps:');
    disp(datestr(datetimeUTC(1:5),'dd-mmm-yyyy HH:MM:SS.fff'));
    disp('Last 5 time stamps:');
    disp(datestr(datetimeUTC(end-5:end),'dd-mmm-yyyy HH:MM:SS.fff'));
    disp('Data duration (hours): ');
    disp(num2str(diff(datenum([datetimeUTC(1);datetimeUTC(end)]))*24, '%.2f'))
    data.Time = [];
    %Creates a matlab DateNumber field that can be displayed using:
    %datestr(data.DN,'dd-mmm-yyyy HH:MM:SS.fff')
    data.DN = datenum(datetimeUTC);
    % Create data.Date and data.Time for CATS PRH processing    
    data.Date = floor(data.DN);
    data.Time = data.DN-floor(data.DN);
    %Checks for incomplete lines from data
    skip = isnan(data.Acc1)|isnan(data.Acc2)|isnan(data.Acc3);
    if sum(skip)~=0; error('Missing ACC Data, need to address'); end 
    disp(['Number of missing data points (Pressure): ' num2str(sum(isnan(data.Pressure)))]);
    oiD = diff(data.DN);
    BAD = find(diff(data.DN)*24*60*60>median(oiD*24*60*60)*1.5);
    disp(['Number of time gaps: ' num2str(size(BAD,1))]);
    if length(BAD) >=1
        for i = 1:size(BAD,1)
            data.Date(BAD(i)+1) = data.Date(BAD(i));
            data.DN(BAD(i)+1) = data.Date(BAD(i)+1) + data.Time(BAD(i)+1);
            disp(datestr(data.DN(BAD(i):BAD(i)+2),'dd-mmm-yyyy HH:MM:SS.fff'));
            if isempty(notes); notes = 'Bad data at t = '; end
            notes = [notes datestr(DN(BAD(i)),'mm/dd HH:MM:SS.fff') ', '];
        end
    end
    
    %determine sampling rate
    fs = 1/((data.DN(end)-data.DN(1))*24*60*60/length(data.DN));
    disp(['Sampling Rate is: ', num2str(fs)]);
    %check if sampling rate is noninteger and if so round
    if abs(round(fs)-fs)>.01;  disp('Sampling rate is non-integer, but using rounded value'); disp('fs was '); disp(num2str(fs)); fs = round(fs); disp('fs is now '); disp(num2str(fs)); 
    else
    %if the sampling rate is <.01 difference from integer value, we still need
    %to round to integer and display
    disp('Unrounded fs was '); disp( num2str(fs)); fs = round(fs); disp('fs is now '); disp( num2str(fs));
    end   
   
    if abs(fs-(1/((data.DN(end) - data.DN(1))*24*60*60/size(data,1))))>.1; 
        disp('Sampling rate is non-integer, but using rounded value'); 
    end     

    %Plot the difference between subsequent lines of the datenumber to
    %check for time errors. The difference should be consistent, meaning
    %the plot should look like a solid block.
    figure(2); clf; 
    annotation('textbox', [0, .05, 0, 0], 'string', 'This figure should look like a solid block. If it does not, look for time errors','FitBoxToText','on')
    plot(diff(data.DN));
    
    clearvars BAD oiD skip datetimeUTC
    % Remove unneeded variables
    try data(:,'OneMinuteLightLevel') = []; catch end
    try data(:,'SmoothedLightLevel') = []; catch end
    try data(:,'Events') = []; catch end
    try data(:,'BatteryVoltage') = []; catch end
    try data(:,'Wet_Dry') = []; catch end
%   Create repeated values for variables sampled at lower rates (e.g.Temp, Light)
% Temp
    % Create Temp at fs Hz by repeating 1Hz values
    % Temp1Hz is the actual sampling rate of Temp
    Temp1Hz = data.Temp(~isnan(data.Temp));%1hz vector of Temp
    TempDN = data.DN(~isnan(data.Temp));
    % Determine tempFS
    THz = 1/((TempDN(end)-TempDN(1))*24*60*60/length(TempDN));
    disp(['Temperature sampling rate is: ', num2str(THz) 'Hz']);
    %check if sampling rate is noninteger and if so round
    if abs(round(THz)-THz)>.01;  disp('Temperature sampling rate is non-integer, but using rounded value'); disp('THz was '); disp(num2str(THz)); 
        THz = round(THz); disp('THz is now '); disp(num2str(THz)); 
    else
    %if the sampling rate is <.01 difference from integer value, we still need
    %to round to integer and display
        THz = round(THz);
    end   
    % TempRep will have the 1Hz data repeated fs times per second
    TempRep = nan(size(data.Temp));  
    for i = 1:round(fs)  
        TempRep(i:round(fs):length(Temp1Hz)*round(fs)) = Temp1Hz; 
    end
    % find NaN's at end
    nanStart = find(isnan(TempRep),1,'first'); % location of first NaN
    lastGood = find(~isnan(TempRep),1,'last'); % location of last good temp value
    % throw an error if there are too many NaN's being replaced at end.
    if length(TempRep)-nanStart > fs 
        error('There is an issue with the up/downsampling process of Temp') ;
    else      
        TempRep(nanStart:end) = TempRep(lastGood);
    end
    data.Temp = TempRep(1:length(data.Temp));
    clearvars TempRep TempDN lastGood nanStart;
% Light    
    % Light1Hz is the actual sampling rate of Light
    Light1Hz = data.Light(~isnan(data.Light));%1hz vector of Light
    LightDN = data.DN(~isnan(data.Light));
    % Determine lightFS
    lHz = 1/((LightDN(end)-LightDN(1))*24*60*60/length(LightDN));
    disp(['LightLevel sampling rate is: ', num2str(lHz) 'Hz']);
    %check if sampling rate is noninteger and if so round
    if abs(round(lHz)-lHz)>.01;  disp('Temperature sampling rate is non-integer, but using rounded value');
        disp('lHz was '); disp(num2str(lHz)); 
        lHz = round(lHz); disp('lHz is now '); disp(num2str(lHz)); 
    else
    %if the sampling rate is <.01 difference from integer value, we still need
    %to round to integer and display
        lHz = round(lHz);
    end  
    % LightRep will have the 1Hz data repeated 32 times per second
    LightRep = nan(size(data.Light));  
    for i = 1:round(fs)  
        LightRep(i:round(fs):length(Light1Hz)*round(fs)) = Light1Hz; 
    end
    % find NaN's at end
    nanStart = find(isnan(LightRep),1,'first'); % location of first NaN
    lastGood = find(~isnan(LightRep),1,'last'); % location of last good Light value
    % throw an error if there are too many NaN's being replaced at end.
    if length(LightRep)-nanStart > fs 
        error('There is an issue with the up/downsampling process of Light') ;
    else      
        LightRep(nanStart:end) = LightRep(lastGood);
    end
    data.Light = LightRep(1:length(data.Light));
    clearvars LightRep LightDN lastGood nanStart;
% TempInternal
    % Create TempInternal at fs Hz by repeating 1Hz values
    % Temp1Hz is the actual sampling rate of TempInternal
    Temp1Hz = data.TempInternal(~isnan(data.TempInternal));%1hz vector of Temp
    TempDN = data.DN(~isnan(data.TempInternal));
    % Determine tempFS
    T1Hz = 1/((TempDN(end)-TempDN(1))*24*60*60/length(TempDN));
    disp(['Internal Temperature sampling rate is: ', num2str(T1Hz) 'Hz']);
    %check if sampling rate is noninteger and if so round
    if abs(round(T1Hz)-T1Hz)>.01;  disp('Internal Temperature sampling rate is non-integer, but using rounded value'); disp('T1Hz was '); disp(num2str(T1Hz)); 
        T1Hz = round(T1Hz); disp('T1Hz is now '); disp(num2str(T1Hz)); 
    else
    %if the sampling rate is <.01 difference from integer value, we still need
    %to round to integer and display
        T1Hz = round(T1Hz);
    end   
    % TempRep will have the 1Hz data repeated fs times per second
    TempRep = nan(size(data.TempInternal));  
    for i = 1:round(fs)  
        TempRep(i:round(fs):length(Temp1Hz)*round(fs)) = Temp1Hz; 
    end
    % find NaN's at end
    nanStart = find(isnan(TempRep),1,'first'); % location of first NaN
    lastGood = find(~isnan(TempRep),1,'last'); % location of last good temp value
    % throw an error if there are too many NaN's being replaced at end.
    if length(TempRep)-nanStart > fs 
        error('There is an issue with the up/downsampling process of TempInternal') ;
    else      
        TempRep(nanStart:end) = TempRep(lastGood);
    end
    data.TempInternal = TempRep(1:length(data.TempInternal));
    clearvars TempRep TempDN lastGood nanStart;    
% TempDepthInternal
    % Create TempDepthInternal at fs Hz by repeating 1Hz values
    % Temp1Hz is the actual sampling rate of TempDepthInternal
    Temp1Hz = data.TempDepthInternal(~isnan(data.TempDepthInternal));%1hz vector of TempDepthInternal
    TempDN = data.DN(~isnan(data.TempDepthInternal));
    % Determine tempFS
    TDIHz = 1/((TempDN(end)-TempDN(1))*24*60*60/length(TempDN));
    disp(['Depth Sensor Temperature sampling rate is: ', num2str(TDIHz) 'Hz']);
    %check if sampling rate is noninteger and if so round
    if abs(round(TDIHz)-TDIHz)>.01;  disp('Depth Sensor Temperature sampling rate is non-integer, but using rounded value'); disp('TDIHz was '); disp(num2str(TDIHz)); 
        TDIHz = round(TDIHz); disp('TDIHz is now '); disp(num2str(TDIHz)); 
    else
    %if the sampling rate is <.01 difference from integer value, we still need
    %to round to integer and display
        TDIHz = round(TDIHz);
    end   
    % TempRep will have the 1Hz data repeated fs times per second
    TempRep = nan(size(data.TempDepthInternal));  
    for i = 1:round(fs)  
        TempRep(i:round(fs):length(Temp1Hz)*round(fs)) = Temp1Hz; 
    end
    % find NaN's at end
    nanStart = find(isnan(TempRep),1,'first'); % location of first NaN
    lastGood = find(~isnan(TempRep),1,'last'); % location of last good temp value
    % throw an error if there are too many NaN's being replaced at end.
    if length(TempRep)-nanStart > fs 
        error('There is an issue with the up/downsampling process of TempDepthInternal') ;
    else      
        TempRep(nanStart:end) = TempRep(lastGood);
    end
    data.TempDepthInternal = TempRep(1:length(data.TempDepthInternal));
    clearvars TempRep TempDN lastGood nanStart;    
    
    disp(head(data,5));
    % Save Hzs (determined from data)
    Hzs = struct('accHz',fs,'gyrHz',fs,'magHz',fs,'pHz',fs,'lHz',lHz,'GPSHz',fs,'UTC',0,'THz',THz,'TDIHz',TDIHz);
    
    disp('Section 2 (Data Processing) finished');
        
%% 2a. Truncate Data (optional)
% To reduce file sizes, select the period of data to save by zooming in and selecting the boundaries of
% time on the whale.
% output: tagon and tagoff variables in workspace
    Depth = decdc(data.Pressure,round(fs/4));% if resolution is too fine, the slope of p is hard to determine
    tagonI = find(Depth>.25); % CATS tags start before they are on the whale, TDR starts in the water. 
    [s, e] = consec(tagonI);
    tagonI = s(find(e-s>=fs,1,'first')); % first time you have 1 second deeper than 1 m
    tagonI = find(Depth(1:tagonI)<min(Depth(1:tagonI))+.1,1,'last'); % changed a value
    tagoff = find(Depth<3,1,'last');
    tagoff = round(tagoff*fs/4); 
    if tagonI ~= 1
        tagonI = tagonI*round(fs/4); %adjust back to undecimated values
    end
    Depth = data.Pressure;   

    figure(41); clf;
    %maximize the figure
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
    plot(Depth); hold on; set(gca,'ydir','rev'); ylim([0 max(data.Pressure)]);
    title('Truncate Data');
    pat = nan(size(tagonI));
    for i = 1:length(tagonI)
        pat(i) = rectangle('position',[tagonI(i) -10 tagoff(i)-tagonI(i) 600],'facecolor',[1 1 .4]); %,'parent');%,AX{jj}); patch([tagonI(i) tagoff(i) tagoff(i) tagonI(i)],[-10 -10 600 600],[255 255 100]/255);
    end
    oi = get(gca,'children'); oi=[oi(i+1:end); oi(1:i)];
    set(gca,'children',oi);
    annotation('textbox', [0.05, .05, 0, 0], 'string', 'Yellow box highlights the portion of data that will be saved','FitBoxToText','on')

    TEX = text(.5,max(Depth),'LEFT click to CHANGE the boundary of the data, RIGHT click to ZOOM in. Press Enter when done.','verticalalignment','bottom','fontweight','bold','horizontalalignment','left');
    [x,~,button] = ginput(1);
    zoom reset; fact = 20;
    while ~isempty(button)
        % regraph
        delete(pat); pat = nan(size(tagonI));
        for i = 1:length(tagonI)
            pat(i) = rectangle('position',[tagonI(i) -10 tagoff(i)-tagonI(i) 600],'facecolor',[1 1 .4]); 
        end;  oi = get(gca,'children'); oi=[oi(i+1:end); oi(1:i)]; set(gca,'children',oi);
        if button == 3
            fact = fact/2;
            xlim([x-fact*60*fs x+fact*60*fs]); %surrounds the point by fact minutes
            delete (TEX); TEX = text(min(get(gca,'xlim')),max(get(gca,'ylim')),'LEFT click to CHANGE a boundary for tag on/tag off, RIGHT click to ZOOM in, press "a" to ADD a set of boundaries, press "d" to delete a boundary','verticalalignment','bottom','fontweight','bold','horizontalalignment','left');
            [x,~,button] = ginput(1);
            if button == 3; continue; end
        end
        if button == 1
            [~,e] = min(abs([tagonI tagoff]-x));
            if e>length(tagonI); tagoff(e-length(tagonI)) = x;
            else tagonI(e) = x; end
            zoom out; delete(pat); pat = nan(size(tagonI));
            for i = 1:length(tagonI);
                pat(i) = rectangle('position',[tagonI(i) -10 tagoff(i)-tagonI(i) 600],'facecolor',[1 1 .4]); 
            end;  oi = get(gca,'children'); oi=[oi(i+1:end); oi(1:i)]; set(gca,'children',oi);
            delete (TEX); TEX = text(min(get(gca,'xlim')),max(get(gca,'ylim')),'LEFT click to CHANGE a boundary for tag on/tag off, RIGHT click to ZOOM in, press "a" to ADD a set of boundaries, press "d" to delete a boundary','verticalalignment','bottom','fontweight','bold','horizontalalignment','left');
            [x,~,button] = ginput(1); fact = 20; continue;
        end
    end
    %
    if tagonI < 100
        tagonI = 1;
    end
    tagoff = round(tagoff, 0);
    
    oi = zeros(size(Depth));
    for i = 1:length(tagonI)
        oi(tagonI(i):tagoff(i)) = 1;
    end
    tagon = logical(oi);
    % remove the data the was off-whale
    data(~tagon,:)=[]; 
    clearvars pat fact oi useold x TEX s e button
    disp('Section 2a (Truncate Data) finished');
    
%% 3. Save Imported File
% Create Adata, Atime and Hzs
Atime = data.DN;
Adata = [data.Acc1 data.Acc2 data.Acc3];

 try
    save([fileloc filename(1:end-4) '.mat'],'data','events','Adata','Atime','Hzs');
    if ~isempty(notes); save([fileloc filename(1:end-3) 'mat'],'notes','-append'); end
     if ~isempty(lastwarn)
         error(lastwarn);
     end
 catch %v7.3 allows for bigger files, but makes a freaking huge file if used when you don't need it
     save([fileloc filename(1:end-4) '.mat'],'data','events','Adata','Atime','Hzs','-v7.3');
     if ~isempty(notes); save([fileloc filename(1:end-3) 'mat'],'notes','-append'); end
     disp('Made a version 7.3 file in order to include all');
 end
 disp('Section 3 (Save) finished');
%% 4. Plot Depth Profile and Save as BMP
if exist('pconst','var'); p = (data.Pressure-pconst)*pcal; else p = data.CorrectedDepth; p2 = data.Pressure; end
if any(data.Pressure>10)&&min(data.Pressure)<.5*max(data.Pressure) % plot something if it went underwater
    if ~exist('fs','var'); fs = round(1/((data.DN(50)-data.DN(49))*60*60*24)); end
    oi = max(1, find(p>3,1,'first')-60*fs):min(find(p>4,1,'last')+60*fs,length(p)); 
    if ~isempty(oi)
        oi1 = oi(1:floor(length(oi)/3)); 
        oi2 = oi(floor(length(oi)/3):2*floor(length(oi)/3)); 
        oi3 = oi(2*floor(length(oi)/3):end);
        
        time = data.DN;
        figure(1); clf; 
        set(1, 'units','normalized','outerposition',[0 0 1 1]);
        subplot(3,1,1); plot(time(oi1),p(oi1)); ylim([-5 max(p(oi1))]); xlim(time([oi1(1) oi1(end)])); set(gca,'xticklabel',datestr(get(gca,'xtick'),'mm/dd/yy HH:MM'),'ydir','rev');
        subplot(3,1,2); plot(time(oi2),p(oi2)); ylim([-5 max(p(oi2))]); xlim(time([oi2(1) oi2(end)])); set(gca,'xticklabel',datestr(get(gca,'xtick'),'mm/dd/yy HH:MM'),'ydir','rev');
        subplot(3,1,3); plot(time(oi3),p(oi3)); ylim([-5 max(p(oi3))]); xlim(time([oi3(1) oi3(end)])); set(gca,'xticklabel',datestr(get(gca,'xtick'),'mm/dd/yy HH:MM'),'ydir','rev');
        legend('Depth','Location','best');
        saveas(1,[fileloc filename(1:end-4) '_TDR3.bmp']);
    end
end
disp('Section 4 (Plot) finished');
