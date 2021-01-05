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
%
for i = 1:length(vidDN)
    if isnan(vidDN(i)) || vidDN(i)>DN(find(tagondec,1,'last')); continue; end
    [~,a] = min(abs(DN-vidDN(i))); [~,b] = min(abs(DN-(vidDN(i)+vidDurs(i)/24/60/60)));
    a = max(a,find(tagondec,1)); b = min(b,find(tagondec,1,'last'));
    dbps = dbpks(dbpks>=a & dbpks<=b); dbps = dbps(2:end-1);
    minioffsets = nan(size(dbps));
    for ii = 1:length(dbps)
        [~,aa] = min(abs(jpks-dbps(ii)));
        minioffsets(ii) = (jpks(aa)-dbps(ii))/fs;
    end
    offsets(i) = nanmean(minioffsets(abs(minioffsets)<maxoffset));
    if isnan(offsets(i)) && i>1; offsets(i) = offsets(i-1); exT = ' (taken from last value)'; else exT = ''; end
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
pause;

% if okay with offsets, continue, else ctrl c to go back and fix them
oldvar = struct('vidDN',vidDN,'camondec',camondec,'audondec',audondec,'DB',DB,'camon',camon,'audon',audon);
for i = 1:length(vidDN)
    if isnan(offsets(i)) || offsets(i) == 0; continue; end
    [~,a] = min(abs(DN-vidDN(i))); b = min(round(a+vidDurs(i)*fs),length(Depth));
    vidDN(i) = vidDN(i)+offsets(i)/24/60/60;
    [~,a2] = min(abs(DN-vidDN(i))); b2 = min(round(a2+vidDurs(i)*fs),length(Depth));
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
camon = false(size(tagon));
for ii = 1:ofs/fs
    camon(ii:ofs/fs:length(camondec)*ofs/fs) = camondec;
end
audon = false(size(tagon));
for iii = 1:ofs/fs
    audon(iii:ofs/fs:length(audondec)*ofs/fs) = audondec;
end    
    
    Jig = [JX JY JZ J];
disp(['Mean offset: ' num2str(nanmean(offsets))]);
flownoise = DB;
save([fileloc filename(1:end-4) 'Info.mat'],'flownoise','camondec','audondec','camon','audon','vidDN','offsets','oldvar','Jig','-append');