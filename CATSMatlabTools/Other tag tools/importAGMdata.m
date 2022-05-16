function [data, Adata, Atime, Hzs] = importAGMdata(FS)
%
% import Little Leonardo data
% David Cade
% version 5.31.2021
% Goldbogen Lab
% Stanford University

% Matlab packages required: Signal Processing Toolbox
%
% script format:
% importCATSdata(); % prompts to select file
% importCATSdata(FS);% if FS is entered, final data will have FS sample
% rate, else data will be at the same sample rate as the accelerometer
dbstop if error
  [file,fileloc] = uigetfile('*.csv','select csv file with AGM data');
  
%   if exist([depthloc 'Acc.txt'],'file'); accloc = depthloc; accfile = 'Acc.txt';
%   else
%       c = pwd; cd(depthloc)
%       [accfile,accloc] = uigetfile('*.txt','select txt file with accelerometer data');
%       cd(c);
%   end
  
  data = readtable([fileloc file]);
  headers = data.Properties.VariableNames;
  headers(contains(headers,'Depth','ignorecase',true)) = {'Pressure'};
   headers(contains(headers,'Temp','ignorecase',true)) = {'Temp'};
   headers(contains(headers,'Comp','ignorecase',true)) = {'Comp1' 'Comp2' 'Comp3'};  
        headers(contains(headers,'Gyr','ignorecase',true)) = {'Gyr1' 'Gyr2' 'Gyr3'};
        headers(find(contains(headers,'Acc','ignorecase',true),3)) = {'Acc1' 'Acc2' 'Acc3'};
    data.Properties.VariableNames = headers;
    DN = datenum(data.Timestamp);%+18/24/60/60;
    data.Date = floor(DN);
    data.Time = DN-floor(DN);
    Afs = 100;
Adata = [data.Acc1 data.Acc2 data.Acc3]; Atime = DN; 

Hzs = struct();%'accHz',accHz,'gyrHz',gyrHz,'magHz',magHz,'pHz',pHz,'lHz',lHz,'GPSHz',GPSHz,'UTC',UTC,'THz',THz,'T1Hz',T1Hz,'datafs',datafs);
Hzs.accHz = Afs; Hzs.pHz = 1; Hzs.THz = 1; Hzs.datafs = Afs; Hzs.magHz = Afs; Hzs.gyrHz = Afs;
ODN = data.Date(1)+data.Time(1);
   save([fileloc file(1:end-4) '.mat'],'data','Adata','Atime','Hzs','ODN');


