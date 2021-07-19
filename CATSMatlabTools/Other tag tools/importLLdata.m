function [data, Adata, Atime, Hzs] = importLLdata(FS)
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

  [depthfile,depthloc] = uigetfile('*.txt','select txt file with depth data');
  
  if exist([depthloc 'Acc.txt'],'file'); accloc = depthloc; accfile = 'Acc.txt';
  else
      [accfile,accloc] = uigetfile('*.txt','select txt file with accelerometer data');
  end
  
  data = readtable([accloc accfile],'headerlines',6);
  data.Properties.VariableNames = {'Acc1' 'Acc2' 'Acc3'};
  ddata = readtable([depthloc depthfile],'headerlines',6);
  depth = interp2length(ddata.Depth,1,100,size(data,1));
  Temp = interp2length(ddata.Temp,1,100,size(data,1));
  Camon = logical(interp2length(ddata.Video,1,100,size(data,1)));
DN = datenum([2021 05 22 12 57 00]);
I = (0:size(data,1)-1)';
DN = DN+I/100/24/60/60;
data.Date = floor(DN);
data.Time = DN-floor(DN);
data.Pressure = depth;
data.Temp = Temp;
data.CamOn = Camon;
Adata = [data.Acc1 data.Acc2 data.Acc3];
Atime = DN;
Hzs = struct();%'accHz',accHz,'gyrHz',gyrHz,'magHz',magHz,'pHz',pHz,'lHz',lHz,'GPSHz',GPSHz,'UTC',UTC,'THz',THz,'T1Hz',T1Hz,'datafs',datafs);
Hzs.accHz = 100; Hzs.pHz = 1; Hzs.THz = 1; Hzs.datafs = 100;
data.Comp1 = nan(size(data.Acc1));
data.Comp2 = data.Comp1; data.Comp3 = data.Comp2;
data.Gyr1 = data.Comp1; data.Gyr2 =data.Comp1; data.Gyr3 = data.Comp1;
   save([depthloc 'data.mat'],'data','Adata','Atime','Hzs');


