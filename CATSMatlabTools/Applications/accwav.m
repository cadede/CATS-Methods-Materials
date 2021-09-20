% ACCWAV - filter accelerometer data and save as .wav for detecting
% vocalizations
% William Oestreich and David Cade
% Goldbogen Lab
% Stanford University
%
% INPUTS: 
% Calling this script will the user to select the truncate.mat and prh.mat
% files for a deployment of interest.
%
% OUTPUTS:
% The script will save a .wav file with number of channels equal to the
% number of accelerometer axes (i.e. 3 channels (x, y, and z) for triaxial
% accelerometer data.

% load tag data 
% PRH file
disp('Select PRH file of the tag deployment')
[prhfilename,prhfileloc]=uigetfile('*prh.mat', 'select PRH file of the tag deployment');
cd(prhfileloc);
% truncate file
disp('Select truncate file of the tag deployment')
[filename,fileloc]=uigetfile('*.mat', 'select truncate file of the tag deployment');
cd(fileloc);
disp('Loading PRH file'); 
load([prhfileloc prhfilename(1:end-3) 'mat']); 
whaleName = INFO.whaleName;
disp('Loading truncate file'); 
load([fileloc filename(1:end-3) 'mat']); 

%calculate sampling freq and round to make sure it is an integer value
Afs = round(1./mean((Atime(50:60)-Atime(49:59))*24*60*60));
    
%apply calibration
names =fieldnames(INFO.CAL);
for ii = 1:length(names)
    eval([names{ii} ' = INFO.CAL.' names{ii} ';']);
end
if ~exist('Acal','var') || isempty(Acal)
    A = (Adata-repmat(aconst,size(Adata,1),1))*acal;
else
    axA = (acal./abs(acal)); axA(isnan(axA)) = 0;
    A = Adata*axA;
    A = (A*diag(Acal.poly(:,1))+repmat(Acal.poly(:,2)',[size(A,1),1]))*Acal.cross;
end

%apply running mean filter
A = A - runmean(A,2*Afs);

% add start time to acc.wav file
t = Atime(1);


%save acc.wav
disp('Saving acc.wav'); 
audiowrite([whaleName,'_acc-',datestr(t,'yyyymmdd-HHMMSS-fff'),'.wav'],A,Afs);