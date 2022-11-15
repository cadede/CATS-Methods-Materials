% rename video files

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
   if vidT
        if  i==1 ; try v = mmread([fileloc filename{i}],[1 10]); dur = v.totalDuration;
            catch; vid2 = VideoReader([fileloc filename{i}]); dur = vid2.Duration;
            end
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