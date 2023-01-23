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
dbstop if error
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
  try Afs = 1000/str2num(heads{1}{find(cellfun(@(x) strcmp(x,'PERIOD'),heads{1}))+1});
  catch fid = fopen([accloc accfile]);
  heads = textscan(fid,'%s',30,'headerlines',4);
  fclose(fid);
  Afs = 1000/str2num(heads{1}{find(cellfun(@(x) strcmp(x,'PERIOD'),heads{1}))+1});
  end
     
  try data.Properties.VariableNames = {'Acc1' 'Acc2' 'Acc3'};
  catch; data.Properties.VariableNames = {'Acc1' 'Acc2' 'Acc3' 'Comp1' 'Comp2' 'Comp3'};
  end
  ddata = readtable([depthloc depthfile]);%,'headerlines',6);
  fid = fopen([depthloc depthfile]);
  dheads = textscan(fid,'%s',10,'headerlines',4);
  fclose(fid);
  try dfs = str2num(dheads{1}{find(cellfun(@(x) strcmp(x,'PERIOD'),heads{1}))+1});
  catch
      fid = fopen([depthloc depthfile]);
  dheads = textscan(fid,'%s',30,'headerlines',4);
  fclose(fid);
  dfs = str2num(dheads{1}{find(cellfun(@(x) strcmp(x,'PERIOD'),heads{1}))+1});
  end
  depth = interp2length(ddata.Depth,dfs,Afs,size(data,1));
  Temp = interp2length(ddata.Temp,dfs,Afs,size(data,1));
  try Camon = logical(interp2length(ddata.Video,1,Afs,size(data,1))); catch; end
  try Speed = interp2length(ddata.Speed,1,Afs,size(data,1)); catch; end
% DN = datenum([2021 05 22 12 57 00]);
try DN = datenum([str2num(heads{1}{3}(1:4)) str2num(heads{1}{4}(1:regexp(heads{1}{4},'/')-1)) str2num(heads{1}{5}) 0 0 0]) ...
        + datenum(heads{1}{8},'HH:MM:SS')-floor(datenum(heads{1}{8},'HH:MM:SS'));
catch; disp('Could not read start time, input starttime as datevec:');
    DV = input('?'); DN = datenum(DV);
end
        
I = (0:size(data,1)-1)';
DN = DN+I/Afs/24/60/60;
data.Date = floor(DN);
data.Time = DN-floor(DN);
data.Pressure = depth;
data.Temp = Temp;
try data.CamOn = Camon; catch; data.CamOn = false(size(data.Acc1)); end
try data.Speed = Speed; catch; data.Speed = nan(size(data.Acc1)); end
Adata = [data.Acc1 data.Acc2 data.Acc3];
Atime = DN;
Hzs = struct();%'accHz',accHz,'gyrHz',gyrHz,'magHz',magHz,'pHz',pHz,'lHz',lHz,'GPSHz',GPSHz,'UTC',UTC,'THz',THz,'T1Hz',T1Hz,'datafs',datafs);
Hzs.accHz = Afs; Hzs.magHz = Afs; Hzs.pHz = dfs; Hzs.THz = dfs; Hzs.datafs = Afs; Hzs.SHz = dfs;
try t = data.Comp1(1)+1; catch; data.Comp1 = nan(size(data.Acc1));
data.Comp2 = data.Comp1; data.Comp3 = data.Comp2;
end
data.Gyr1 = nan(size(data.Acc1)); data.Gyr2 =data.Gyr1; data.Gyr3 = data.Gyr1;
ODN = data.Date(1)+data.Time(1);
   save([depthloc 'data.mat'],'data','Adata','Atime','Hzs','ODN');


