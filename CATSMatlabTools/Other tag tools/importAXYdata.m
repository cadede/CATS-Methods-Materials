function [data, Adata, Atime, Hzs] = importAXYdata(FS,fillnans)
%
% import AXY data from technosmart
% David Cade
% version 7.14.2025
% Goldbogen Lab
% Stanford University

% Matlab packages required: Signal Processing Toolbox
%
% script format:
% importCATSdata(); % prompts to select file
% importCATSdata(FS);% if FS is entered, final data will have FS sample
% rate, else data will be at the same sample rate as the accelerometer
% 
% 
% fillnans if true does a linear interpolation between points for data sampled at
% less than FS, else it repeats those values.

if nargin<2
    fillnans = false;
end

dbstop if error
  [file,fileloc] = uigetfile('*.csv','select csv file with AGM data');
  UTC = input('UTC offset? ');
  
%   if exist([depthloc 'Acc.txt'],'file'); accloc = depthloc; accfile = 'Acc.txt';
%   else
%       c = pwd; cd(depthloc)
%       [accfile,accloc] = uigetfile('*.txt','select txt file with accelerometer data');
%       cd(c);
%   end
  
  data = readtable([fileloc file],'DatetimeType','text');
  data.Date = datenum(data.Date,'dd/mm/yyyy'); % may change this in future downloads
  headers = data.Properties.VariableNames;
  headers(contains(headers,'Depth','ignorecase',true)) = {'Pressure'};
  headers(contains(headers,'Press','ignorecase',true)) = {'Pressure'};
   headers(contains(headers,'Temp','ignorecase',true)) = {'Temp'};
   try headers(contains(headers,'Comp','ignorecase',true)) = {'Comp1' 'Comp2' 'Comp3'};
   catch headers(contains(headers,'Mag','ignorecase',true)) = {'Comp1' 'Comp2' 'Comp3'};
   end
        try headers(contains(headers,'Gyr','ignorecase',true)) = {'Gyr1' 'Gyr2' 'Gyr3'}; catch; end
        try headers(find(contains(headers,'Acc','ignorecase',true),3)) = {'Acc1' 'Acc2' 'Acc3'};
        catch
            headers(cellfun(@(x) strcmp(x,'X'),headers)) = {'Acc1'};
            headers(cellfun(@(x) strcmp(x,'Y'),headers)) = {'Acc2'};
            headers(cellfun(@(x) strcmp(x,'Z'),headers)) = {'Acc3'};
        end
    data.Properties.VariableNames = headers;
    try DN = datenum(data.Timestamp)+UTC/24;%+18/24/60/60;
    catch; DN = datenum(data.Date)+datenum(data.Time)+UTC/24;
    end
    data.Date = floor(DN);
    data.Time = DN-floor(DN);
    Afs = round(1/mean(diff(DN(10:60))*24*60*60));
Adata = [data.Acc1 data.Acc2 data.Acc3]; Atime = DN; 

Hzs = struct();%'accHz',accHz,'gyrHz',gyrHz,'magHz',magHz,'pHz',pHz,'lHz',lHz,'GPSHz',GPSHz,'UTC',UTC,'THz',THz,'T1Hz',T1Hz,'datafs',datafs);
Hzs.accHz = Afs; 
oAfs = Afs;
if exist('FS','var') 
    resamp = true;
    Afs = FS;
    toresamp = {'Acc1' 'Acc2' 'Acc3'};
    tempfs = floor(Hzs.accHz/FS)*FS; if tempfs == 0; tempfs = FS; end
    tempDN = (DN(1):1/tempfs/24/60/60:DN(end))';
    %     data2 = data(1,:); data2(1,:) = [];
    for ii = 1:3
        tempdf = tempfs/FS;
        VAR = timeseries(data.(toresamp{ii}),data.Date+data.Time);
        VAR = resample(VAR,tempDN);
        if any(isnan(VAR.data)); VAR.data = edgenans(inpaint_nans(VAR.data)); end
        %         if isempty(data2); data2.Acc1 = nan(size(decdc(VAR.data,tempdf))); end
        if ii == 1
            newDN = DN(1)+(0:length(decdc(VAR.data,tempdf))-1)/24/60/60/FS;
            data2 = table(floor(newDN)',newDN'-floor(newDN'),'VariableNames',{'Date','Time'});
        end
        data2.(toresamp{ii}) = decdc(VAR.data,tempdf);
    end
    
   else
    resamp = false;
end

Hzs.datafs = Afs; 
strs = {'p','T','mag','gyr'};
datastrs = {'Pressure','Temp','Comp','Gyr'};
datalength = size(data,1);

for h =1:length(strs)
    if h>2; d = [datastrs{h} '1']; numvar = 3; else; d = datastrs{h}; numvar = 1; end
    s = [strs{h} 'Hz'];
    try data.(d); catch; continue; end
    Hz = round(sum(~isnan(data.(d)))/length(data.(d))*oAfs); 
    Hzs.(s) = Hz;
    for k = 1:numvar
        if k >1; d(end) = num2str(k); end
        if ~fillnans
            I = find(~isnan(data.(d)));
            for n = 1:oAfs/Hz-1
                data.(d)(I+n) = data.(d)(I);
            end
        else
            data.(d) = edgenans(fixgaps(data.(d)));
        end
        if resamp
            I = 1:size(data2,1);
            I = round(I*oAfs/FS);
            % if I(end)>size(data,1); I(end) = size(data,1); end
            data2.(d)= data.(d)(I);
        end
    end
end
data(datalength+1:end,:) = [];
Hzs.UTC = UTC
ODN = data.Date(1)+data.Time(1);
if resamp; data = data2; end
try data.TagID = []; catch; end
try data.Metadata = []; catch; end
lastwarn('');
    try
        save([fileloc file(1:end-4) '.mat'],'data','Adata','Atime','Hzs','ODN');
        if ~isempty(lastwarn)
            error(lastwarn);
        end
    catch %v7.3 allows for bigger files, but makes a freaking huge file if used when you don't need it
          save([fileloc file(1:end-4) '.mat'],'data','Adata','Atime','Hzs','ODN','-v7.3');
        disp('Made a version 7.3 file (large data format) file');
    end



