% This uses the RMS flownoise with the RMS jiggle and tries to line up
% peaks in the signal, treating each audio file separately.  This is for a
% sweet of time when the CATS tags had a video offset time from the data
% that needed to be accounted for.


load([fileloc filename(1:end-4) 'Info.mat'],'tagondec');
DB = flownoise;
JJ = J; JJ(isnan(JJ)) = 0; JJ = runmean(JJ,fs);
D = DB; D(isnan(D)|isinf(D)) = min(D(~isinf(D)));  D = runmean(D,fs);

% figure (200); clf;
dbpks = peakfinder(D,5);
jpks = peakfinder(JJ,5);

% plot(D); hold on; plot(JJ,'g'); plot(dbpks,D(dbpks),'ro'); plot(jpks,JJ(jpks),'rs')
oldVDN = vidDN;
offsets = zeros(size(vidDN));
% if onlyAud; offsets = 
%
for i = 1:length(vidDN)
    if isnan(vidDN(i)) || vidDN(i)>DN(find(tagondec,1,'last')) || (onlyAud && ~isempty(frameTimes{i})); continue; end
    [~,a] = min(abs(DN-vidDN(i))); [~,b] = min(abs(DN-(vidDN(i)+vidDurs(i)/24/60/60)));
    a = max(a,find(tagondec,1)); b = min(b,find(tagondec,1,'last'));
    dbps = dbpks(dbpks>=a & dbpks<=b);% dbps = dbps(2:end-1);
    minioffsets = nan(size(dbps));
%     jpks = peakfinder(JJ(a:b),5)+a-1;
    for ii = 1:length(dbps)
        [~,aa] = min(abs(jpks-dbps(ii)));
        minioffsets(ii) = (jpks(aa)-dbps(ii))/fs;
    end
    offsets(i) = nanmean(minioffsets(abs(minioffsets)<maxoffset));
    if isnan(offsets(i)) && i>1; try  offsets(i) = offsets(find(offsets(1:i-1)~=0,1,'last'));exT = ' (taken from last value)';
        catch; offsets(i) = 0; exT = ''; end;  else exT = ''; end
    figure(200+i); clf;
     set(gcf,'windowStyle','docked');
    s1 = subplot(211); 
    plot(a:b,D(a:b)); hold on; plot(a:b,JJ(a:b),'g'); plot(dbps,D(dbps),'ro'); plot(jpks,JJ(jpks),'rs');
    s2 = subplot(212); 
    plot((a:b)+offsets(i)*fs,D(a:b)); hold on; plot(a:b,JJ(a:b),'g'); plot(dbps+offsets(i)*fs,D(dbps),'ro'); plot(jpks,JJ(jpks),'rs');
    title(['Offset ' num2str(i) ' = ' num2str(offsets(i)) ' ' exT]);
    linkaxes([s1 s2],'x','y');
    set(gca,'xlim',[a b])
    [~,~,button] = ginput(1);
    if isempty(button)
        continue;
    elseif button == 48
        offsets(i) = 0;
    end

end

disp(offsets);
% if okay with offsets, continue, else ctrl c to go back and fix them
disp('If you''re okay with the observed offsets, press enter to continue and save the results, else press ctrl-c to go back and fix the offsets, prior versions saved in "oldvar" variable')

pause;

if ~exist('camondec','var'); camondec = camon; end
if ~exist('audondec','var'); audondec = camon; end
if ~exist('frameTimes','var') && ~nocam; load([fileloc filename(1:end-4) 'movieTimes.mat'],'frameTimes'); end
oldvar = struct('vidDN',vidDN,'camondec',camondec,'audondec',audondec,'DB',DB,'camon',camon,'audon',audon);
camondec = false (size(DN)); audondec = camondec;
for i = 1:length(vidDN)
%     if isnan(offsets(i)) || offsets(i) == 0; continue; end
    if isnan(vidDN(i)); continue; end
    [~,a] = min(abs(DN-vidDN(i))); b = min(round(a+vidDurs(i)*fs),length(DN));
    vidDN(i) = vidDN(i)+offsets(i)/24/60/60;
    [~,a2] = min(abs(DN-vidDN(i))); b2 = min(round(a2+vidDurs(i)*fs),length(DN));
    if ~isempty(frameTimes{i})
        camondec(a:b) = false;
        camondec(a2:b2) = true;
    else
        audondec(a:b) = false;
        audondec(a2:b2) = true;
    end
    oi = DB(a:b);
    DB(a:b) = nan;
    if b2-a2+1==length(oi)
    DB(a2:b2) = oi;
    elseif a2+length(oi)-1<=length(DB)
       DB(a2:a2+length(oi)-1) = oi;
    else
        DB(a2:end) = oi(1:length(DB)-a2+1);
    end

end
if length(camon)~=length(camondec)
camon = false(size(tagon));
for ii = 1:ofs/fs
    camon(ii:ofs/fs:length(camondec)*ofs/fs) = camondec;
end
audon = false(size(tagon));
for iii = 1:ofs/fs
    audon(iii:ofs/fs:length(audondec)*ofs/fs) = audondec;
end   
else camon = camondec; audon = audondec;
end
    
%     Jig = [JX JY JZ J];
disp(['Mean offset: ' num2str(nanmean(offsets))]);
flownoise = DB;
save([fileloc filename(1:end-4) 'Info.mat'],'flownoise','camondec','audondec','camon','audon','vidDN','offsets','oldvar','-append');