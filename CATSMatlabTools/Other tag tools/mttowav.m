%code to convert .mt files from Bprobes and Acousondes to .wav files for
%analysis in Triton
%mfm 2011-06-20
%
% put time stamp in filename for Triton (ie wav file time and LTSAs)
% converted p units from mPa (after MTRead) back to 16-bit A/D counts because
% wavwrite was incorrectly adjusting amplitude values for 32-bit files.
% smw 2012-12-12
%--------------------------------------------------------------------------
clear all;close all; clc
%TAG TYPE -----------------------------------------------------------------
prompt1={'Enter tag type (1=Acousonde, 2=Bprobe)','Enter Speed: (1=Slow, 2=High)'};
inl = inputdlg(prompt1);
flag = str2num(inl{1});
sptype = str2num(inl{2});
%DIRECTORY OF FILES TO PROCESS---------------------------------------------
start_path = 'F:\TRANSDEC2012\';
if flag==1 %acousdonde
 if sptype == 1
 inpath = uigetdir(start_path,'Select Directory for MT files');
 cd(inpath);D=dir('*S*.MT');
 elseif sptype == 2
 inpath = uigetdir(start_path,'Select Directory for MT files');
 cd(inpath);D=dir('*H*.MT');
 end
 outpath = uigetdir(inpath,'Select Directory for WAV files');
elseif flag==2
 inpath = uigetdir(start_path);cd(inpath);D=dir('*_Sound_*.MT');
end
%PROCESS MT FILES (loop)---------------------------------------------------
disp('please wait ...')
for ii = 1:length(D); %ii=1;
 [p,header,info] = MTRead([inpath '/' D(ii).name]);
 % p2 = p; % convert units to uPa from mPa
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % 32-bit wavwrite below doesn't preserve the mPa values, it
 % scales/normalizes in some unknown way, so convert mPa back to A/D counts
 calmax = str2num(header.calmax);
 calmin = str2num(header.calmin);
 % signed integer:
 bitmin = -(2^(str2num(header.samplebits)-1));
 bitmax = (2^(str2num(header.samplebits)-1)) - 1;
 multiplier = (calmax-calmin)/(bitmax-bitmin);
 %p = (p - bitmin).*multiplier + calmin; % this converts A/D counts to mPa
 %in MTRead.m
 p2 = (p - calmin)./multiplier + bitmin; % convert units from mPa back to A/D counts

 %----------------------------------------------------------------------
 % HEADER INFORMATION
 yy =str2num(header.year); mm=str2num(header.month); dd=str2num(header.day);
 hh=str2num(header.hours); m=str2num(header.minutes); ss=str2num(header.seconds);
 strt = datenum(yy, mm, dd, hh, m, ss);
 n=info.srate;
 msamp = (length(p2)/n)/60;
 ftstr = datestr(strt,'yymmdd-HHMMSS'); % file time string
 %write out wavfiles
 % pout = int32(p2);
 pout = int16(p2);
 max(p2);min(p2);
 if flag==1
 % outfileA =(D(ii).name(1:8));
 outfileA =[D(ii).name(1:8),'_',ftstr]; % put date/time string in file name for Triton 25
 % wavwrite(pout,n,32,outfileA)
 f = fullfile(outpath,outfileA);
 f = [f '.wav'];
 wavwrite(pout,n,16,f)
 elseif flag==2
 outfileB =(D(ii).name(1:22));
 outfileB = [outfileB '.wav'];
 wavwrite(pout,n,32,outfileB)
 end

 outdat = char(datestr(strt,'mm/DD/YYYY HH:MM:ss'));
 if flag == 1
 fprintf('%s %s \n',outdat, outfileA);
 elseif flag==2
 fprintf('%s %s \n',outdat, outfileB);
 end
 clear p pout p2 %header info outfile outdate
end
disp('Done')