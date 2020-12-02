function CATS2TrackPlot(head,pitch,roll,tagon,DN,fs,track,isgeo,whaleName,FS,filedest)

% David Cade
% version 3.19.18
% Goldbogen Lab
% Stanford University
% head- animal heading (in radians from true north)
% pitch- animal pitch in radians
% roll- positive is to the right, flip sign of roll for dtags
% tagon- an index variable indicating when the tag is on the whale
% DN- a set of date numbers for each sample point
% fs- sampling rate
% track- pseudotrack (3 columns of [N E depth])
% isgeo- is this a georeferenced pseudotrack?  (usually set to false since georeferenced tracks can have sharp curves
% whaleName- name for final file
% FS- final sampling rate (1.25 is good)
% filedest- folder location for the final file


folder = '';
if nargin<1 || isempty(head);
    [filename,fileloc] = uigetfile('*.*','Select prhfile file ');
    load([fileloc filename]);
    try [track,ptrack] = gtrack(pitch,head,p,fs,speed.JJ,tagon,DN,GPS,GPSerr);
        if nargin>=8 && ~isempty(isgeo) && isgeo
            track = track;
        else track = ptrack; isgeo = false;
        end
    catch
        track = nan(size(p));
        head = fixgaps(head); roll = fixgaps(roll); pitch = fixgaps(pitch);
        head(isnan(head)) = 0; roll(isnan(roll)) = 0; pitch(isnan(pitch)) = 0;
        sp = speed.JJ;
        sp2 = sp; sp2(isnan(sp2)) = min(sp2);
        sp2 = runmean(sp2,2*fs);
        sp2(p<=pthresh) = nan;
        sp2 = fixgaps(sp2);
        sp(p<=pthresh) = sp2(p<=pthresh); sp = fixgaps(sp); sp(isnan(sp)) = min(sp);
        t = ptrack(pitch(tagon),head(tagon),p(tagon),fs,[],sp(tagon));
        tt = t(:,1);
        t(:,1) = t(:,2); t(:,2) = tt;
        isgeo = false;
    end
  
end
    



if nargin<11
     [~,filedest] = uigetfile('*.*','Select file in the folder to place output file (where trackplot is), press cancel to use prh folder');
end
if nargin<10 || isempty(FS)
    FS = 1.25;
end


% fl = fileloc;
% fn = prhfile;
% fl2 = destloc;
% 
% % takes CATS prh files and removes the nans for use in TrackPlot
% 
% raw_smooth = .8; % Colin recommends a ratio of .8 of a second, I don't know for 10 Hz files .1 might be more appropriate (no smoothing)
% folder = '';%'files\'; % if you want to keep your track plot files in a subfolder of the trackplot files, put that here.
% 
% if exist('fl','var')&&exist('fn','var'); fileloc = fl; filename = fn; if exist('fl2','var'); outloc = fl2; else outloc = fl; end; 
% else
%     [filename,fileloc]=uigetfile('*mat', 'select prh files','multiselect','on');
%     cf = pwd;
%     cd(fileloc);
%     [~,outloc] = uigetfile('*.*','Select file in the folder to place output file (where trackplot is), press cancel to use prh folder');
%     cd(cf);
%     if ~any(outloc); outloc = fileloc; end
% end
% % load files
% 
% 
% if ischar(filename)
%     filename = {filename};
% end

% trackplot requires n

head = fixgaps(head); roll = fixgaps(roll); pitch = fixgaps(pitch);
head(isnan(head)) = 0; roll(isnan(roll)) = 0; pitch(isnan(pitch)) = 0;

dcm = angle2dcm(head(tagon),pitch(tagon),roll(tagon));
s = struct('v',squeeze(num2cell(dcm,[1 2])));
N = arrayfun(@(y) mtimes([1 0 0],y.v),s,'uniformoutput',false);
N = vertcat(N{:});
N(:,3) = -N(:,3);
oi = N(:,1); N(:,1) = N(:,2); N(:,2) = oi;
D = arrayfun(@(y) mtimes([0 0 -1],y.v),s,'uniformoutput',false);
D = vertcat(D{:});
D(:,3) = -D(:,3);
% D(:,2) = -D(:,2);
oi = D(:,1); D(:,1) = D(:,2); D(:,2) = oi;
M = [track(tagon,:) N D ones(size(N,1),2)];
M = runmean(M,round(1/FS*fs));
if any (any(isnan(M)))
    error('Nans present');
end
 M = decdc(M,fs/FS);

if isgeo
    geo = '';
else geo = '';
end

M(:,3) = -M(:,3); %pressure is negative in trackplot;

outfile = [whaleName 'tp' geo '.vec'];
txtfile = [whaleName 'tp' geo '.txt'];
nrec = size(M,1);
fidOut = fopen([filedest outfile],'w');
% Write header line
fprintf(fidOut,'Vers2 %0.2f\r\n',1/FS);
fprintf(fidOut,'NRecords %d\r\n',nrec);
fprintf(fidOut,'X \tY \tZ \tFx \tFy \tFz \tDx \tDy \tDz \tDacc \tNacc \r\n');
for i = 1:nrec
    fprintf(fidOut, '%5.3f  \t%5.3f  \t%5.3f  \t%5.3f  \t%5.3f  \t%5.3f  \t%5.3f \t%5.3f \t%5.3f \t%5.3f \t%5.3f \r\n', M(i,1),M(i,2),M(i,3),M(i,4),M(i,5),M(i,6),M(i,7),M(i,8),M(i,9),M(i,10),M(i,11));
end
fclose(fidOut);
fidOut = fopen([filedest txtfile],'w');
fprintf(fidOut,'smooth %0.0f\r\n',1);
fprintf(fidOut,'vec_file %s\r\n',[folder outfile]);
if strcmpi(whaleName(1:2),'mn'); AN = 'HUMPBACK'; elseif strcmpi(whaleName(1:2),'bw')||strcmpi(whaleName(1:2),'bp')||strcmpi(whaleName(1:2),'bb')||strcmpi(whaleName(1:2),'be'); AN = 'BLUE'; elseif strcmp(whaleName(1:2),'oo')||strcmp(whaleName(1:2),'rt'); AN = 'DOLPHIN'; else AN = 'HUMPBACK'; end
fprintf(fidOut,'animal %s\r\n',AN);
DV = datevec(DN(find(tagon,1)));
fprintf(fidOut,'tagon_hms %02.0f %02.0f %02.0f\r\n',DV(1,4:6));% we set trackplot above to only be for when the tag is on
% fprintf(fidOut,'raw_sample %.0f\r\n',fs);
% fprintf(fidOut,'speed %.1f\r\n',1);
% fprintf(fidOut,'declination %.0f\r\n',dec);

% fprintf(fidOut,'tagon_hms %02.0f %02.0f %02.0f\r\n',DV(1,4:6));% we set trackplot above to only be for when the tag is on
% A = filename{fn}(1:2);
% if strcmp(A,'mn'); AN = 'HUMPBACK'; elseif strcmp(A,'bw')||strcmp(A,'bp')||strcmp(A,'Bp')||strcmp(A,'Bw'); AN = 'BLUE'; elseif strcmp(A,'oo'); AN = 'DOLPHIN'; else AN = 'HUMPBACK'; end
% fprintf(fidOut,'animal %s\r\n',AN);
fprintf(fidOut,'name %s\r\n',whaleName);
fprintf(fidOut,'end');
fclose(fidOut);

end


function y=fixgaps(x)
% FIXGAPS Linearly interpolates gaps in a time series
% YOUT=FIXGAPS(YIN) linearly interpolates over NaN
% in the input time series (may be complex), but ignores
% trailing and leading NaN.
%

% R. Pawlowicz 6/Nov/99

y=x;

bd=isnan(x);
gd=find(~bd);

bd([1:(min(gd)-1) (max(gd)+1):end])=0;


y(bd)=interp1(gd,x(gd),find(bd)); 
end