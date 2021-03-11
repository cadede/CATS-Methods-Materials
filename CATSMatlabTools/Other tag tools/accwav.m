% ACCWAV - filter accelerometer data and save as .wav for detecting
% vocalizations
% William Oestreich and David Cade
% Goldbogen Lab
% Stanford University
%
% INPUTS: 
% Calling this script will the user to select the truncate.mat, prh.mat, 
% and calibration.mat files for a deployment of interest.
%
% OUTPUTS:
% The script will save a .wav file with number of channels equal to the
% number of accelerometer axes (i.e. 3 channels (x, y, and z) for triaxial
% accelerometer data.


% load tag data  
% truncate file
disp('select truncate file of the tag deployment')
[filename,fileloc]=uigetfile('*.*', 'select truncate file of the tag deployment');
cd(fileloc);
disp('Loading truncate file'); 
load([fileloc filename(1:end-3) 'mat'], 'Adata','Atime'); 
% PRH file
disp('select PRH file of the tag deployment')
[filename,fileloc]=uigetfile('*.*', 'select PRH file of the tag deployment');
cd(fileloc);
disp('Loading PRH file'); 
load([fileloc filename(1:end-3) 'mat'],'DN','INFO','tagon','fs','p'); 
whaleName = INFO.whaleName;
% calibration file
disp('select calibration file for tag');
[filename,fileloc]=uigetfile('*.*', 'select calibration file for tag'); 
disp('Loading calibration file')
load([fileloc filename(1:end-3) 'mat'])

%calculate sampling freq and round to make sure it is an integer value
Afs = round(1./mean((Atime(50:60)-Atime(49:59))*24*60*60));
    
%apply calibration
if exist('Acal','var') && ~isempty(Acal)
    axA = (acal./abs(acal)); axA(isnan(axA)) = 0;
    A = Adata*axA;
    A = (A*diag(Acal.poly(:,1))+repmat(Acal.poly(:,2)',size(A,1),1))*Acal.cross;
else
    A = (Adata-repmat(aconst,size(Adata,1),1))*acal;
end

%apply running mean filter
A = A - runmean(A,2*fs);

%save acc.wav
audiowrite([whaleName,'_acc','.wav'],A,Afs);