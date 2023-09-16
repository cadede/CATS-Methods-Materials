% rename video files: cell 1 to rename stitched video clips (at the conclusion of writing prh data onto video files)
% , cell 2 can be used to rename non-cats video files so that they have the video start
% times embedded in the file name (this is usually done before prh process even starts).

% If video files have been made with data and video, and if they have been
% stitched together with Adobe Premiere, they will have an appendix usually
% like: ... (stitched clip).mp4.  This removes that appendix and replaces
% it with the timestamp of the start of the video.  This can also be used
% for files without the (stitched clip) appendix as long as they have the
% video number embedded.

[filename,fileloc] = uigetfile('*.mp4','Grab files to rename','multiselect','on');
if ~iscell(filename); filename = {filename}; end
d = dir(fileloc); D = {d.name};
prh = ~cellfun(@isempty,cellfun(@(x) regexp(x,'prh.mat$'),D,'uniformoutput',false));
for i = 1:length(filename)
    w = ~cellfun(@isempty,cellfun(@(x) regexp(x,['^' getWhaleID(filename{i})]),D,'uniformoutput',false));
    if any(w&prh); vidT = true; oi = D(prh&w); load([fileloc oi{1}],'vidDN','vidDurs','INFO'); else vidT = false; end
    I = strfind(filename{i},' (Stitched Clip)');
    vn = min(strfind(filename{i}, '('));
    dI = [min(strfind(filename{i},')')) max(strfind(filename{i},'-'))]; vn2 = max([dI(2) vn]); dI = dI(dI>vn);
    dI2 = dI(1); 
    dI = min(dI);
    vn = str2num(filename{i}(vn+1:dI-1));
    vn2 = str2num(filename{i}(vn2+1:dI2-1));
    if ~isempty(I); f = [filename{i}(1:I-1) '.mp4']; else f = filename{i}; end
    if vidT;
        if  i==1 ; v = mmread([fileloc filename{i}],[1 10]); dur = v.totalDuration;
            if dur < vidDurs(vn); vidDN(vn) = vidDN(vn) + (vidDurs(vn)-dur)/24/60/60; end
        end
        f = [f(1:end-4) ' (' datestr(vidDN(vn),'mm.dd HHMMSS') ').mp4'];
    end
    movefile([fileloc filename{i}],[fileloc f]);
end

disp(['All files renamed for ' INFO.whaleName]);
    
% Rename in a LOOP
% for n = 1:numel(names)
%     oldname = [d names{n}];
%     newname = sprintf('%s%0*s',d,mLen, names{n});
%     dos(['rename "' oldname '" "' newname '"']); % (1)
% end

%%
% renames all selected movies using cats video format: 
% (e.g. CATSCAM42-yyyymmdd-HHMMSS-fff-###). Takes a start time of a single
% video as input and then assumes all other videos are immediately before
% or after it based on video length

[filename,fileloc] = uigetfile('*.mp4','Grab all files to rename','multiselect','on');
if ~iscell(filename); filename = {filename}; end
cd (fileloc)
f1 = uigetfile('*.mp4','select video for which you know the start time','multiselect','on');
fnum1 = find(cellfun(@(x) strcmp(f1,x),filename));
disp(['Movie selected is number ' num2str(fnum1) ' in the list. Press enter to label with this file number, else enter number you wish to label it with: ']);
fnum2 = input('movie number of selected movie? ');
if ~isempty(fnum2); fnum = fnum2; else; fnum = fnum1; end
dn1 = input('Enter start time of selected movie as date vector: ');
lab = input('What prefix should be used in file names (e.g. ''CATS CAM 3'')? ');
movN = fnum - fnum1+1:length(filename)+fnum-fnum1;
vidDurs = zeros(size(movN));
xlswrite([fileloc 'oldmovienames.xlsx'],filename')
newnames = filename;
for i = [fnum1:-1:1 fnum1+1:length(filename)]
    vid = VideoReader([fileloc filename{i}]);
    vidDurs(movN(i)) = vid.Duration;
    if i <=fnum1
        dn = datenum(dn1) - sum(vidDurs(movN(i:fnum1-1)))/60/60/24;
    else
        dn = datenum(dn1) + sum(vidDurs(movN(fnum1:i-1)))/60/60/24;
    end
    ds = datestr(dn,'yyyymmdd-HHMMSS-fff');
    str = [lab '-' ds '-' sprintf('%05u',movN(i))];
    
    movefile([fileloc filename{i}],[fileloc str '.mp4']);
    disp(['File ' filename{i} ' renamed to ' str '.mp4']);
    newnames{i} = [str '.mp4'];
end
xlswrite([fileloc 'movienamechanges.xlsx'],[{'old names'} {'new names'}; [filename' newnames']])

disp('All files renamed');

%% undo last
for i = 1:length(filename)
    movefile([fileloc newnames{i}],[fileloc filename{i}]);
end
