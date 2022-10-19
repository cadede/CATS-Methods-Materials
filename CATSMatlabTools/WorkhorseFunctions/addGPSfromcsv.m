%% 13. Combine fastGPS data with PRH file
% based on Add_SDA_GPS by James Fahlbusch (c) 2019

function addGPSfromcsv(fileloc,INFO)
try D = dir([fileloc '/raw/']);
    D = {D.name}; D = D(3:end);
    iscsv = cellfun(@(x) strcmp(x(end-6:end),'gps.csv'),D,'uniformoutput',false);
    if sum(cellfun(@isempty,iscsv)) == 0
       D = D(cellfun(@(x) x,iscsv));
    else
        D = D(~isempty(iscsv));
    end
    if length(D) == 1; filenameGPS = D{1}; filelocGPS = [fileloc '/raw/']; else error('can''t find file'); end
catch
    cdir = pwd; cd(fileloc);
    [filenameGPS,filelocGPS] = uigetfile('*.csv','select GPS.csv file for this deployment');
    cd (cdir);
end
% tag 70 is fastGPS 1231, 71 is 1232;
% if str2num(INFO.tagnum) == 70; gpsID = 1231; elseif str2num(INFO.tagnum) == 71; gpsID = 1232;
% else gpsID = input('Input fastGPS ID # for this deployment: ');
% end
% disp(['gpsID = ' num2str(gpsID)]);

% optsG = detectImportOptions([filelocGPS filenameGPS]);
% optsG.ExtraColumnsRule = 'ignore';
dataGPS = readtable([filelocGPS filenameGPS],'ReadVariableNames',false);
% headersG = {'GPSdate' 'GPSTime' 'GPSDate' 'GPSTime' 'Lat' 'Long' 'Alt' 'Err1' 'Err2' 'Err3' 'numsats' 'numsats2'}; %optsG.VariableNames;
headersG = {'GPSDate' 'GPSTime' 'UTCoffset_s' 'Latitude' 'Longitude' 'Alt1','Alt2','DistTrav'};
% csv has is exported from R script so all values have a common format
% dataGPS = readtable([filelocGPS filenameGPS],optsG);
% need to process with milliseconds and without seperetely and combine
%dateMill = datetime(dataGPS.Date, 'InputFormat', 'MM/dd/yyyy HH:mm:ss.SS');
try dataGPS.Properties.VariableNames([1 2 6 15 16 17 18 31]) = headersG;
catch
    dataGPSN = readtable([filelocGPS filenameGPS],'ReadVariableNames',true,'Range','1:2');
    headersN = dataGPSN.Properties.VariableNames;
    headersN{find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Date'),headersN,'uniformoutput',false)),1)} = 'GPSDate';
    headersN{find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Time'),headersN,'uniformoutput',false)),1)} = 'GPSTime';
    dataGPS.Properties.VariableNames = headersN;
    
end
% dataGPS(dataGPS.ID~=gpsID,:) = [];
% dateFull = datenum(dataGPS.GPSDate,'dd.mm.yyyy')+datenum(dataGPS.GPSTime,'HH:MM:SS.fff')-floor(datenum(dataGPS.GPSTime)); %datetime(dataGPS.DateTimeUTC, 'inputformat', 'M/d/yyyy HH:mm:ss');
GPSDN = datenum(dataGPS.GPSDate)+datenum(dataGPS.GPSTime);%-floor(datenum(dataGPS.GPSTime)); %datetime(dataGPS.DateTimeUTC, 'inputformat', 'M/d/yyyy HH:mm:ss');
% dataGPS.datetimeUTC = dateFull;
% dataGPS(:,[2 7:10])=[]; % remove unused columns
dataGPS.DNUTC = GPSDN;% not sure what UTC offset does, but would also have to adjust data to match-dataGPS.UTCoffset_s/24/60/60;%datenum(dataGPS.datetimeUTC); % Convert to number format if needed
disp('First 10 Timestamps of GPS Data (UTC)');
try disp(datestr(dataGPS.DNUTC(1:10),'dd-mmm-yyyy HH:MM:SS.fff')); catch
    disp(datestr(dataGPS.DNUTC(1:end),'dd-mmm-yyyy HH:MM:SS.fff'));
end
% clearvars dateFull headersG optsG
% cd(filelocGPS)
%Reload PRH file
D = dir(fileloc); D = {D.name};
% D = [{fileloc} D];
% oi = dir(fileloc); oi2 = {oi.name};
% F = [{oi2(~vertcat(oi.isdir))} F];
prhfile = D{~cellfun(@isempty,cellfun(@(x) strfind(x,'prh.mat'),D,'uniformoutput',false))};
% prhfile = F{j}(~cellfun(@isempty,prhfile));

% [prhfile,prhfileloc]=uigetfile('*.mat', 'select PRH file');
load([fileloc prhfile]);

% disp(['TagTurnedOn: ' datestr(DN(1),'mm/dd/yy HH:MM:SS.fff')]);
disp(['TagOnAnimal (local): ' datestr(DN(find(tagon,1)),'mm/dd/yy HH:MM:SS.fff')]);
disp(['TagOffAnimal (local): ' datestr(DN(find(tagon,1,'last')),'mm/dd/yy HH:MM:SS.fff')]);
disp(['EndData (local): ' datestr(DN(end),'mm/dd/yy HH:MM:SS.fff')]);

%add GPS to PRH file
% try
    GPSoffset = INFO.UTC;%+INFO.timedif;
% catch
%     GPSoffset = INFO.timedif;
% end
dataGPS.GPSDN = dataGPS.DNUTC+GPSoffset/24; %This converts GPS to same offest as PRH (local Time)

Lat = runmean(dataGPS.Latitude,10); Long = runmean(dataGPS.Longitude,10);
Lat = Lat(1:10:end); Long = Long(1:10:end);
dataGPS = dataGPS(1:10:end,:); dataGPS.Lat = Lat; dataGPS.Long = Long;
try repI = find(diff(dataGPS.DistTrav)==0)+1;
catch
    repI = find(abs(diff(dataGPS.Latitude))>.0001|abs(diff(dataGPS.Longitude))>.0001)+1;
end
% repI = find(abs(diff(dataGPS.Lat)) < .0001 & abs(diff(dataGPS.Long)) <.0001)+1;
dataGPS(repI,:) = [];
% Remove points from recovery
dataGPS(dataGPS.GPSDN>max(DN)|dataGPS.GPSDN<min(DN),:) = []; % remove GPS locations outside of PRH data (sometimes hits from recovery)

figure(641); clf; plot(DN,-p,dataGPS.GPSDN,0,'rs')
set(gca,'xticklabel',datestr(get(gca,'xtick'),'mm/dd/yy HH:MM:SS'));
saveas(641,[fileloc INFO.whaleName '_DiveGPS.bmp']);

GPSdata = [dataGPS.Lat dataGPS.Long];
GPSall = nan(length(DN),2); % Create an empty matrix to store gpsdata
% GPSall(find(tagon,1,'first'),1:2) = GPS(1,1:2); %this adds the deployment location from the XLSX file
GPSDN = dataGPS.GPSDN;
GPSerr = nan(max(length(GPSDN),length(DN)),3);
% add GPS locations to the correct location of the PRH
for b =1:length(dataGPS.GPSDN)
    [~,a] = min(abs(DN-GPSDN(b)));
    GPSall(a,1:2) = GPSdata(b,1:2);
%     GPSerr(a,1:3) = [dataGPS.Err1(b) dataGPS.Err2(b) dataGPS.Err3(b)];
end
GPSall(1,1:2) = GPS(1,1:2);
%save a UTM version of the location data
% UTM = nan (size(GPSall));
% [x,y] = deg2utm(GPSall(~isnan(GPSall(:,1)),1),GPSall(~isnan(GPSall(:,1)),2));
% UTM(~isnan(GPSall(:,1)),:) = [x y];
GPS = GPSall;
% Save data to PRH file
save([fileloc prhfile],'GPS','GPSerr','-append');
% Display a map of positions (now redundant to main text)
% whalename = INFO.whaleName;% prhfile(1:strfind(prhfile,' ')-1);
% [fig,ax] = plotMap(GPS,DN,tagon,10);
% surfs = findsurfacings(p,fs,tagon,60);
% xs = get(ax,'xlim'); ys = get(ax,'ylim');
% text(xs(1)+diff(xs)/20,ys(2)-diff(ys)/20,{['Number of surfacings: ' num2str(length(surfs))]; ['Number of GPS hits on surfacings: ' num2str(size(GPS(tagon&~isnan(GPS(:,1))),1))]},'parent',ax);
% set(fig,'units','normalized','outerposition',[0 0 1 1])
% % Save Map
% if ~exist([fileloc '\QL\'],'dir'); mkdir([D{j} '\QL\']); end
% % savefig(fig,[D{j} '\QL\' whalename ' Map.fig']);
% savefig(fig,[fileloc '\QL\' whalename ' Map.fig']);
% saveas(fig,[fileloc whalename ' Map.bmp']);
% timedif = GPSoffset;


% % Save GPS Data
% % Add Deployment Location
% Tnew = dataGPS(1,:);
% Tnew.Lat = GPS(find(tagon,1,'first'),1);  Tnew.Long = GPS(find(tagon,1,'first'),2);
% Tnew.datetimeUTC = datestr(DN(find(~isnan(GPS(:,1)),1,'first'))); Tnew.DNUTC = DN(find(~isnan(GPS(:,1)),1,'first'))-(timedif/24);
% Tnew.GPSDN = DN(find(~isnan(GPS(:,1)),1,'first'));
% dataGPS = [Tnew;dataGPS];
% dataGPS.dateTime = datestr(dataGPS.GPSDN,'mm/dd/yy HH:MM:SS.fff'); % add local time column
% writetable(dataGPS,[prhfileloc whalename ' GPS.xlsx'],'Sheet',1);
% 
% %create KML
% lat= GPS(~isnan(GPS(:,1)),1); lat(1) = [];
% long= GPS(~isnan(GPS(:,1)),2); long(1) = [];
% tstart = DN(~isnan(GPS(:,1))); tstart(1) = [];
% tstop = DN(~isnan(GPS(:,1))); tstop(1) = [];
% %convoluted system, but it works
% kmlstr = '';
% for ii = 1:length(lat)
%     kmlstr = [kmlstr ge_point(long(ii),lat(ii),0,'timeSpanStart',datestr(tstart(ii)-INFO.timedif/24,'yyyy-mm-ddTHH:MM:SSZ'),'timeSpanStop',datestr(tstop(ii)-INFO.timedif/24,'yyyy-mm-ddTHH:MM:SSZ'))];
% end
% KML = ge_folder('Points',kmlstr);
% ge_output([prhfileloc whalename '.kml'],KML);