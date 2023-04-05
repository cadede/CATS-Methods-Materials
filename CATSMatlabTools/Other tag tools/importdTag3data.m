% make dtag prh 2 cats data (requires d3 tag tools)
clearvars -except CAL
% first load in appropriate cal file 
% load data from swv files
% set these valu
folder = 'D:\Tag Data\EE\From Fleur\ONR Energetics of Exposure\ONR proposal\DATA\Beaked Whales\zc18_185a\';
fileprefix = 'zc185a';
% ODN = datenum([2015 07 31 8 4s6 52]);
df = 1;
name = 'zc18_185a';

 d3 = readd3xml([folder fileprefix '003.xml']);
 ODN = datenum(d3.CFG{1}.TIME,'yyyy,mm,dd,HH,MM,SS');
 disp(['Start Time: ' datestr(ODN)]);
 X=d3readswv(folder,fileprefix,df);
% 
 [ch_names,descr] = d3channames(X.cn);
 [p,CAL2,pfs,T] = d3calpressure(X,CAL,'full');
 figure (2); clf; ax=  plotyy(1:length(p),-p,1:length(p),T); legend('Pressure','Temperature'); set(ax(2),'ylim',[0 45]);
 pause;
 CAL = CAL2;
 [M,CAL2,mfs,Mz] = d3calmag(X,CAL,'full');
 mfs = mfs(1);
 figure(2); clf; plot(1:size(M,1),[M sqrt(sum(M.^2,2))]);
 pause;
 CAL = CAL2;
 [A,CAL2,afs] = d3calacc(X,CAL,'full');
 figure(2); clf; plot(1:size(A,1),[A sqrt(sum(A.^2,2))]);
 CAL = CAL2;

 minfs = min(X.fs); disp(['Will decimate data to ' num2str(minfs) ' Hz, Acc data will be saved separately at ' num2str(afs) ' Hz'])

%%
DN = (ODN:1/afs/24/60/60:ODN+(size(A,1)-1)/afs/24/60/60)';
Adata = A; Atime = DN;
% magI = ~cellfun(@isempty,cellfun(@(x) strfind(x,'MAG'),ch_names,'uniformoutput',false)) & cellfun(@isempty,cellfun(@(x) strfind(x,'TST'),ch_names,'uniformoutput',false));


Acc1 = decdc(A(:,1),afs/minfs);
Acc2 = decdc(A(:,2),afs/minfs);
Acc3 = decdc(A(:,3),afs/minfs);
Pressure = decdc(p,pfs/minfs);
Comp1 = decdc(M(:,1),mfs/minfs);
Comp2 = decdc(M(:,2),mfs/minfs);
Comp3 = decdc(M(:,3),mfs/minfs);
Temp = decdc(T,pfs/minfs);
datafs = minfs;

DN = (ODN:1/minfs/24/60/60:ODN+(size(Acc1,1)-1)/minfs/24/60/60)';
Date = floor(DN); Time = DN-floor(DN);
Hzs = struct('accHz',min(datafs,afs),'gyrHz',nan,'magHz',min(datafs,mfs),'pHz',min(datafs,pfs),'lHz',nan,'GPSHz',nan,'UTC',0,'THz',min(datafs,pfs),'T1Hz',nan,'datafs',datafs);
data = table(Date,Time,Pressure,Acc1,Acc2,Acc3,Comp1,Comp2,Comp3,Temp);

d3initialCAL = CAL;
save([folder name '.mat'],'data','Hzs','ODN','Adata','Atime','d3initialCAL')

% figure; plot(DN,A)
% acal = [1 0 0; 0 1 0; 0 0 -1];
% aconst = [0 0 0];
% magcalon = [1 0 0; 0 1 0; 0 0 -1]; magcaloff = magcalon; 
% magconston = [0 0 0]; magconstoff = magconston;
% save([folder 'DTAGcal1.mat'],'acal','aconst','magcalon','magcaloff','magconston','magconstoff')