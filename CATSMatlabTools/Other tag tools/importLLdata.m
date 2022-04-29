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
      c = pwd; cd(depthloc)
      [accfile,accloc] = uigetfile('*.txt','select txt file with accelerometer data');
      cd(c);
  end
  
  data = readtable([accloc accfile],'headerlines',6);
  fid = fopen([accloc accfile]);
  heads = textscan(fid,'%s',10,'headerlines',4);
  fclose(fid);
  try data.Properties.VariableNames = {'Acc1' 'Acc2' 'Acc3'};
  catch; data.Properties.VariableNames = {'Acc1' 'Acc2' 'Acc3' 'Comp1' 'Comp2' 'Comp3'};
  end
  ddata = readtable([depthloc depthfile]);%,'headerlines',6);
  depth = interp2length(ddata.Depth,1,100,size(data,1));
  Temp = interp2length(ddata.Temp,1,100,size(data,1));
  try Camon = logical(interp2length(ddata.Video,1,100,size(data,1))); catch; end
  try Speed = interp2length(ddata.Speed,1,100,size(data,1)); catch; end
% DN = datenum([2021 05 22 12 57 00]);
try DN = datenum([str2num(heads{1}{3}(1:4)) str2num(heads{1}{4}(1:regexp(heads{1}{4},'/')-1)) str2num(heads{1}{5}) 0 0 0]) ...
        + datenum(heads{1}{8},'HH:MM:SS')-floor(datenum(heads{1}{8},'HH:MM:SS'));
catch; disp('Could not read start time, input starttime as datevec:');
    DV = input('?'); DN = datenum(DV);
end
        
I = (0:size(data,1)-1)';
DN = DN+I/100/24/60/60;
data.Date = floor(DN);
data.Time = DN-floor(DN);
data.Pressure = depth;
data.Temp = Temp;
try data.CamOn = Camon; catch; data.CamOn = false(size(data.Acc1)); end
Adata = [data.Acc1 data.Acc2 data.Acc3];
Atime = DN;
Hzs = struct();%'accHz',accHz,'gyrHz',gyrHz,'magHz',magHz,'pHz',pHz,'lHz',lHz,'GPSHz',GPSHz,'UTC',UTC,'THz',THz,'T1Hz',T1Hz,'datafs',datafs);
Hzs.accHz = 100; Hzs.pHz = 1; Hzs.THz = 1; Hzs.datafs = 100;
try t = data.Comp1(1)+1; catch; data.Comp1 = nan(size(data.Acc1));
data.Comp2 = data.Comp1; data.Comp3 = data.Comp2;
end
data.Gyr1 = nan(size(data.Acc1)); data.Gyr2 =data.Gyr1; data.Gyr3 = data.Gyr1;
   save([depthloc 'data.mat'],'data','Adata','Atime','Hzs');


