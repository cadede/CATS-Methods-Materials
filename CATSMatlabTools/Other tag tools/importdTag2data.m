% make dtag prh 2 cats data (requires d3 tag tools)
clearvars -except CAL DEPLOY TAGON OTAB TAGLOC
% first load in appropriate cal file and using settagpath and d3loadcal to load "CAL" and "DEPLOY" (if a prh has been created
% already in a dtag format)
% loads data from swv files
% set these values

% data need to be in a folder structure such that swv files are in (e.g.)
% [folder \zc10\zc10_272a\ with "name" as the deployment ID
folder = 'G:\Shared drives\Goldbogen Transfer\Energetics of Exposure\All prh files\SoCal beaked\';
fileprefix = 'zc272a';
% ODN = datenum([2015 07 31 8 4s6 52]);
df = 1;
name = 'zc10_272a';
useexistingcalibration = true;

if useexistingcalibration; str = 'none';
else; str = 'full';
end

%  d3 = readd3xml([folder fileprefix '003.xml']);
%  ODN = datenum(d3.CFG{1}.TIME,'yyyy,mm,dd,HH,MM,SS');
try if size(TAGON,1)>size(TAGON,2); TAGON = TAGON'; end; 
    ODN = datenum(TAGON); catch; try ODN = datenum(DEPLOY.TAGON); 
catch; ODN = input('Tag start time as date vector? '); ODN = datenum(ODN);
end; end
 disp(['Start Time: ' datestr(ODN)]);
 settagpath('cal',folder);
settagpath('prh',folder);
settagpath('audio',folder);
 settagpath('raw',folder) % add raw directory to paths; rawdir is the directory to put the raw sensor files e.g., 'c:/tag/data/raw'
% if length(fileprefix)~=9;
%     if length(name) == 9; warning('DTAG2 can only handle 9 character deployment IDs. Renaming swv files with "name" ID');
%     D = dir(folder); D = {D(:).name};
%     I = ~cellfun(@isempty,cellfun(@(x) strfind(x,'.swv'),D,'uniformoutput',false));
%     I = find(I);
%     for ii = I;
%         oldfile = D{ii}
%         newfile = [name oldfile(length(fileprefix)+1:end)]
%         movefile([folder oldfile],[folder newfile]);
%     end
%     else; error('DTAG 2 can only handle 9 character deployment names');
%     end
% end

[X,fs,ch_names] = swvread(name,[],df) ; % read the .swv files
% saveraw(X,fs) % save a raw file%
%  [ch_names,descr] = d3channames(X.cn);
 [p,T,CAL2] = calpressure(X,CAL,str);
 pfs = fs;
 figure (2); clf; ax=  plotyy(1:length(p),-p,1:length(p),T); legend('Pressure','Temperature'); set(ax(2),'ylim',[0 45]);
 pause;
 CAL = CAL2;
 [M,CAL2] = calmag(X,CAL,str);
 mfs = fs;%mfs(1);
 figure(2); clf; plot(1:size(M,1),[M sqrt(sum(M.^2,2))]);
 pause;
 CAL = CAL2;
 [A,CAL2] = calacc(X,p,T,CAL,str);
 afs = fs;
 figure(2); clf; plot(1:size(A,1),[A sqrt(sum(A.^2,2))]);
 CAL = CAL2;

 minfs = fs;
%  minfs = min(X.fs); disp(['Will decimate data to ' num2str(minfs) ' Hz, Acc data will be saved separately at ' num2str(afs) ' Hz'])
disp(['Data will be saved at original sample rate (' num2str(fs) ' Hz)'])
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
Hzs = struct('Afs',afs,'accHz',min(datafs,afs),'gyrHz',nan,'magHz',min(datafs,mfs),'pHz',min(datafs,pfs),'lHz',nan,'GPSHz',nan,'UTC',0,'THz',min(datafs,pfs),'T1Hz',nan,'datafs',datafs);
data = table(Date,Time,Pressure,Acc1,Acc2,Acc3,Comp1,Comp2,Comp3,Temp);

d3initialCAL = struct('CAL',CAL,'OTAB',OTAB,'TAGON',TAGON,'TAGLOC',TAGLOC);
save([folder name(1:4) '\' name '\' name '.mat'],'data','Hzs','ODN','Adata','Atime','d3initialCAL')

% figure; plot(DN,A)
% acal = [1 0 0; 0 1 0; 0 0 -1];
% aconst = [0 0 0];
% magcalon = [1 0 0; 0 1 0; 0 0 -1]; magcaloff = magcalon; 
% magconston = [0 0 0]; magconstoff = magconston;
% save([folder 'DTAGcal1.mat'],'acal','aconst','magcalon','magcaloff','magconston','magconstoff')