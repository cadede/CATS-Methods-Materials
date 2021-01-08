function [data,Adata,Atime,Hzs] = importACOUdata()

% try
% a = getdrives;
% for i = 1:length(a)
%     [result,vol]=system(['vol ' a{i}(1) ':']);
%     if strfind(vol,'DTAGS'); vol = a{i}(1); break; end
% end
% catch
% end

cf = pwd; %try cd([vol ':\']); catch; end
[filename,fileloc]=uigetfile('*.MT', 'select MT files','multiselect','on');
cd(cf);
if ischar(filename)
    filename = {filename};
end
%
D = cell(size(filename)); 
for i = 1:length(filename)
    [d, header, info] = MTRead([fileloc filename{i}]);
    if i == 1; H = header; I = info; 
    else H = [H header]; I = [I info];
    end
    D{i} = d; 
        
% %     p2 = p/max(abs(p));
% %     wavwrite(p,info.srate,[fileloc filename{i}(1:end-2) '2.wav']);
%    d2 = d/max(abs(d)); % use these two lines for audio
%    audiowrite([fileloc filename{i}(1:end-2) 'wav'],d2,round(info.srate))
end
ismag = ~cellfun(@isempty,cellfun(@(x) strfind(x,'Mag'),{H.abbrev},'uniformoutput',false));
isacc = ~cellfun(@isempty,cellfun(@(x) strfind(x,'Acc'),{H.abbrev},'uniformoutput',false));
isH = ~cellfun(@isempty,cellfun(@(x) strfind(x,'Pow'),{H.abbrev},'uniformoutput',false)) | ~cellfun(@isempty,cellfun(@(x) strfind(x,'Freq'),{H.abbrev},'uniformoutput',false)); %is hydrophone
isO = ~(ismag | isacc | isH);

Afs = round(I(find(isacc,1,'first')).srate);
srates = [I.srate];
try fs = round(max(srates(srates<Afs))); catch; fs = round(max(srates<100)); end
if isempty(fs) || (~isempty(Afs) && Afs<50); fs = Afs; end
if isempty(fs); error('SampleRate'); end
ii = find(srates == fs,1,'first');

if any(isacc)
    X = find(isacc&~cellfun(@isempty,cellfun(@(x) strfind(x,'X'),{H.abbrev},'uniformoutput',false)));
    Y = find(isacc&~cellfun(@isempty,cellfun(@(x) strfind(x,'Y'),{H.abbrev},'uniformoutput',false)));
    Z = find(isacc&~cellfun(@isempty,cellfun(@(x) strfind(x,'Z'),{H.abbrev},'uniformoutput',false)));
    if length(X)~=length(Y) || length(Y)~=length(Z); error('More of some Accel files than others...'); end
    XYZ = [X; Y; Z];
    Atime = [];
    Adata = nan(0,3);
    DN = []; %acc = Adata;
    for i = 1:length(X)
        [maxL,n] = max([I(X(i)).nsamp I(Y(i)).nsamp I(Z(i)).nsamp]);
        for k = 1:3; D{XYZ(k,i)}(end+1:maxL) = nan; end % sometimes one mag axis appears to be a sample short
        A1 = I(XYZ(n,i)).datenumber;
        A1 = (A1:1/Afs/24/60/60:A1+(I(XYZ(n,i)).nsamp-1)*1/Afs/24/60/60)';
        Atime = [Atime; A1];
        oi = find((isacc | ismag | isO)&(horzcat(I(:).srate)==fs));
        oi2 = find(abs(vertcat(I(oi).datenumber)-I(XYZ(n,i)).datenumber)<fs/3/24/60/60);
        nsamp = max(vertcat(I(oi(oi2)).nsamp));
        DN = [DN; (A1(1):1/fs/24/60/60:A1(1)+(nsamp-1)/fs/24/60/60)'];
        Adata = [Adata; [D{XYZ(:,i)}]];
%         acc = [acc; decdc([D{XYZ(:,i)}],Afs/fs)];
        if size(Atime,1)~=size(Adata,1); error('Acc size error'); end
    end
end
%
% DN = (I(ii).datenumber:1/fs/24/60/60:I(ii).datenumber+(I(ii).nsamp-1)*1/fs/24/60/60)';
data = table(floor(DN),DN-floor(DN),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),nan(size(DN)),'VariableNames',{'Date','Time','Acc1','Acc2','Acc3','Comp1','Comp2','Comp3'});
for i = 1:length(isO)
    if ~isO(i); continue; end
    if strcmp(H(i).abbrev,'Press'); H(i).abbrev = 'Pressure'; end
    eval(['data.' H(i).abbrev ' = nan(size(data(:,1)));']);
end
toAdd = find(isacc | ismag | isO);
for i = 1:length(toAdd)
    k = find(DN>=I(toAdd(i)).datenumber-fs/3/24/60/60,1);
    oi = find((isacc | ismag | isO)&(horzcat(I(:).srate)==fs));
    oi2 = find(abs(vertcat(I(oi).datenumber)-I(toAdd(i)).datenumber)<fs/3/24/60/60);
    nsamp = max(vertcat(I(oi(oi2)).nsamp));
    if ismember(toAdd(i),X); d = decdc(D{toAdd(i)},Afs/fs); d(end+1:nsamp) = nan; data.Acc1(k:k+nsamp-1) = d(1:nsamp); continue; end
    if ismember(toAdd(i),Y); d = decdc(D{toAdd(i)},Afs/fs); d(end+1:nsamp) = nan; data.Acc2(k:k+nsamp-1) = d(1:nsamp);continue; end
    if ismember(toAdd(i),Z); d = decdc(D{toAdd(i)},Afs/fs); d(end+1:nsamp) = nan; data.Acc3(k:k+nsamp-1) = d(1:nsamp); continue; end
    ofs = round(I(toAdd(i)).srate);
    if abs(floor(fs/ofs)-fs/ofs)>.01; error('sample rates not multiples of each other'); end
    if ofs==fs
        d = D{toAdd(i)};
    elseif ofs<fs
        d = [];
        for j = 1:ceil(fs/ofs)
            d(j:fs/ofs:length(D{toAdd(i)})*fs/ofs+j-1) = D{toAdd(i)};
        end
    else error('sampling rate error');
    end
    d(end+1:nsamp) = nan;
    d = d(1:nsamp);
    if ismag(toAdd(i))
        if strfind(H(toAdd(i)).abbrev,'X'); data.Comp1(k:k+nsamp-1) = d; end
        if strfind(H(toAdd(i)).abbrev,'Y'); data.Comp2(k:k+nsamp-1) = d; end
        if strfind(H(toAdd(i)).abbrev,'Z'); data.Comp3(k:k+nsamp-1) = d; end
        continue;
    end
%     if strcmp(H(i).abbrev,'Press'); H(i).abbrev = 'Pressure'; end
    eval(['data.' H(toAdd(i)).abbrev '(k:k+nsamp-1) = d;']);
    disp([H(toAdd(i)).abbrev ' read successfullly']);
end
data.Gyr1 = nan(size(data.Comp1)); data.Gyr2 = nan(size(data.Comp1)); data.Gyr3 = nan(size(data.Comp1));

accHz = I(find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Acc'),{H.abbrev},'uniformoutput',false)),1)).srate;
gyrHz = nan;
magHz = I(find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Mag'),{H.abbrev},'uniformoutput',false)),1)).srate;
pHz = I(find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Press'),{H.abbrev},'uniformoutput',false)),1)).srate;
THz = I(find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Temp'),{H.abbrev},'uniformoutput',false)),1)).srate;
datafs = fs;
lHz = nan; T1Hz = nan; UTC = nan; GPSHz = nan;
Hzs = struct('accHz',accHz,'gyrHz',gyrHz,'magHz',magHz,'pHz',pHz,'lHz',lHz,'GPSHz',GPSHz,'UTC',UTC,'THz',THz,'T1Hz',T1Hz,'datafs',datafs);
%
lastwarn('');
try
    save([fileloc H(1).stationcode '.' H(1).title '-' H(1).year(3:4) H(1).month H(1).day '-' H(1).sourcesn(end-3:end) '.mat'],'data','Adata','Atime','Hzs');
    if ~isempty(lastwarn)
        error(lastwarn);
    end
catch %v7.3 allows for bigger files, but makes a freaking huge file if used when you don't need it
    save([fileloc H(1).stationcode '.' H(1).title '-' H(1).year(3:4) H(1).month H(1).day '-' H(1).sourcesn(end-3:end) '.mat'],'data','Adata','Atime','Hzs','-v7.3');
    disp('Made a version 7.3 file in order to include all');
end

