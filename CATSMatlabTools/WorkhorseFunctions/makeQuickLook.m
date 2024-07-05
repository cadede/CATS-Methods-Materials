function makeQuickLook(fileloc,makemetadata)
% see required file organization and needs in the last cell of
% MainCATSprhTool.

dbstop if error
if nargin<2; makemetadata = false; end
co = [0 0 1;
      0 0.5 0;
      1 0 0;
      0 0.75 0.75;
      0.75 0 0.75;
      0.75 0.75 0;
      0.25 0.25 0.25];
set(groot,'defaultAxesColorOrder',co);

if ispc; sl = '\'; else sl = '/'; end

if nargin<1
     [~,fileloc] = uigetfile('*.*','Select prh file in the folder where figures and photos and data live');
end
ql = [fileloc 'QL' sl];
%
D = dir(fileloc); isDir = vertcat(D.isdir);
D = {D.name}'; %D = D(3:end);
Dql = dir(ql); Dql = {Dql.name}';
isimag = ~cellfun(@isempty,cellfun(@(x) regexpi(x,'bmp'),D,'uniformoutput',false)) | ~cellfun(@isempty,cellfun(@(x) regexpi(x,'jpg'),D,'uniformoutput',false)) | ~cellfun(@isempty,cellfun(@(x) regexpi(x,'png'),D,'uniformoutput',false)) | ~cellfun(@isempty,cellfun(@(x) regexpi(x,'jpeg'),D,'uniformoutput',false)) | ~cellfun(@isempty,cellfun(@(x) regexpi(x,'tif'),D,'uniformoutput',false));
isimagql = ~cellfun(@isempty,cellfun(@(x) regexpi(x,'bmp'),Dql,'uniformoutput',false)) | ~cellfun(@isempty,cellfun(@(x) regexpi(x,'jpg'),Dql,'uniformoutput',false)) | ~cellfun(@isempty,cellfun(@(x) regexpi(x,'png'),Dql,'uniformoutput',false)) | ~cellfun(@isempty,cellfun(@(x) regexpi(x,'jpeg'),Dql,'uniformoutput',false)) | ~cellfun(@isempty,cellfun(@(x) regexpi(x,'tif'),Dql,'uniformoutput',false));
prhf = D{~cellfun(@isempty,cellfun(@(x) regexp(x,'prh'),D,'uniformoutput',false))&~isimag};
whaleName = prhf(1:regexp(prhf,' ')-1);
isID = ~cellfun(@isempty,cellfun(@(x) regexpi(x,whaleName),D,'uniformoutput',false));
cf = pwd; cd(fileloc);
   blankfileloc = mfilename('fullpath');
        i = strfind( blankfileloc,sl);
         blankfileloc =  blankfileloc(1:i(end));
         bfl = [ blankfileloc 'Structures//'];
         bf = 'NoDataImage.jpg';
        
try
    try map = Dql{~cellfun(@isempty,cellfun(@(x) regexpi(x,'map'),Dql,'uniformoutput',false))&isimagql};
        map = ['QL//' map];
    catch; map = D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'map'),D,'uniformoutput',false))&isimag};
    end
    mapfl = fileloc; 
catch; warning('Could not find map image, select map file or press cancel to use blank file.');
    cd(fileloc);
    [map,mapfl] = uigetfile('*.*','Select map image file');
    cd(cf);
    if isempty(mapfl)||all(mapfl == 0); map = bf; mapfl = bfl; end
end
try prhi =  D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'prh'),D,'uniformoutput',false))&isimag}; catch; end
try prhi =  Dql{~cellfun(@isempty,cellfun(@(x) regexpi(x,'prh'),Dql,'uniformoutput',false))&isimagql}; prhi = ['QL//' prhi]; catch; end
try try tdr = Dql{~cellfun(@isempty,cellfun(@(x) regexpi(x,'tdr'),Dql,'uniformoutput',false))&isimagql};
        tdr = ['QL//' tdr];
    catch; tdr = D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'tdr'),D,'uniformoutput',false))&isimag&isID};
    end
catch; try tdr = D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'tdr'),D,'uniformoutput',false))&isimag}; catch; tdr = []; end
end
try
    try kml = Dql{~cellfun(@isempty,cellfun(@(x) regexpi(x,'kml'),Dql,'uniformoutput',false))&isimagql};
        kml = ['QL//' kml];
    catch; kml = D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'kml'),D,'uniformoutput',false))&isimag};
    end
    kmlfl = fileloc;
catch; warning('Could not find kml image, select kml file or press cancel to use blank file.');
    cd(fileloc);
    [kml,kmlfl] = uigetfile('*.*','Select kml image file');
    cd(cf);
    if isempty(kmlfl)||all(kmlfl == 0); kmlfl = bfl; kml = bf; end
end
try
    try gt = Dql{~cellfun(@isempty,cellfun(@(x) regexpi(x,'ptrack'),Dql,'uniformoutput',false))&isimagql};
        gt = ['QL//' gt];
    catch; gt = D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'ptrack'),D,'uniformoutput',false))&isimag};
    end
    gtfl = fileloc;
catch; warning('Could not find ptrack image, select ptrack file or press cancel to use blank file.');
    cd(fileloc);
    [gt,gtfl] = uigetfile('*.*','Select ptrack image file');
    cd(cf);
    if isempty(gtfl)|| all(gtfl == 0); gt = bf; gtfl = bfl; end
end
try
    try pt = Dql{~cellfun(@isempty,cellfun(@(x) regexpi(x,'geotrack'),Dql,'uniformoutput',false))&isimagql};
        pt = ['QL//' pt];
    catch; pt = D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'geotrack'),D,'uniformoutput',false))&isimag};
    end
    ptfl = fileloc;
catch; warning('Could not find geotrack image, select ptrack or press cancel to use blank file.');
    cd(fileloc);
    [pt,ptfl] = uigetfile('*.*','Select geotrack image file');
    cd(cf);
    if isempty(ptfl)||all(ptfl==0); pt = bf; ptfl = bfl; end
end
try try cam = Dql{~cellfun(@isempty,cellfun(@(x) regexpi(x,'cam'),Dql,'uniformoutput',false))&isimagql};
        cam = ['QL//' cam];
        try cam2 = Dql{~cellfun(@isempty,cellfun(@(x) regexpi(x,'cam2'),Dql,'uniformoutput',false))&isimagql}; cam2 = ['QL//' cam2]; catch; cam2 = []; end
    catch; cam = D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'cam'),D,'uniformoutput',false))&isimag};
        try cam2 = D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'cam2'),D,'uniformoutput',false))&isimag}; catch; cam2 = []; end
    end
    nocam = false;
catch; nocam = true; warning('No Camera images found (can put in QL folder labeled "cam"');
end
pfold = D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'Pic'),D,'uniformoutput',false))&isDir};
D = dir([fileloc pfold]); D = {D.name}';
try
ID = [pfold '//' D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'ID_'),D,'uniformoutput',false))}];
IDfl = fileloc;
[im,imc] = imread([IDfl ID]); try im = ind2rgb(im,imc); catch; end
catch; warning('Could not find ID photo, select ID photo or press cancel to use blank file.');
    cd(fileloc);
    [ID,IDfl] = uigetfile('*.*','Select ID photo');
    cd(cf);
    if isempty(IDfl)||all(IDfl == 0); ID = bf; IDfl = bfl; end
end
try
TAG = [pfold '//' D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'TAG_'),D,'uniformoutput',false))}];
TAGfl = fileloc;
[im,imc] = imread([TAGfl TAG]); try im = ind2rgb(im,imc); catch; end
catch; warning('Could not find photo of tag on animal, select photo or press cancel to use blank file.');
    cd(fileloc);
    [TAG,TAGfl] = uigetfile('*.*','Select tag on photo');
    cd(cf);
    if isempty(TAGfl)||all(TAGfl==0); TAG = bf; TAGfl = bfl; end
    
end
if sum(~cellfun(@isempty,cellfun(@(x) regexpi(x,'drone_'),D,'uniformoutput',false))) > 0
    DRONE = [pfold '//' D{~cellfun(@isempty,cellfun(@(x) regexpi(x,'drone_'),D,'uniformoutput',false))}];
    [im,imc] = imread([fileloc DRONE]); try im = ind2rgb(im,imc); catch; end
    DRONE = im;
    drone = true;
else drone = false; warning('No drone photo found, can put drone folder in "pics" folder');
end
[im,imc] = imread([mapfl map]); try im = ind2rgb(im,imc); catch; end
map = im;
if ~nocam
    [im,imc] = imread([fileloc cam]); try im = ind2rgb(im,imc); catch; end
    cam = im;
    if ~isempty(cam2); [im2,imc] = imread([fileloc cam2]); try im2 = ind2rgb(im2,imc); catch; end;
        cam2 = im2;
        c1 = size(cam); c2 = size(cam2);
        if ~issame(c1,c2); w = max([c1; c2]); cam3 = 255*ones(w(1),c1(2)+c2(2),w(3)); cam3 = uint8(cam3);
            cam3(1:c1(1),1:c1(2),1:c1(3)) = cam;
            cam3(1:c2(1),c1(2)+1:end,1:c2(3)) = cam2;
            cam = cam3;
        else
            cam = [cam cam2];
        end
    end
end
[im,imc] = imread([kmlfl kml]); try im = ind2rgb(im,imc); catch; end
kml = im;
[im,imc] = imread([gtfl gt]); try im = ind2rgb(im,imc); catch; end
gt = im;
[im,imc] = imread([ptfl pt]); try im = ind2rgb(im,imc); catch; end
pt = im;
%note: fix this to be more flexible if no ID file or tag file exists (make
%blank image)
try
[im,imc] = imread([IDfl ID]); try im = ind2rgb(im,imc); catch; end
ID = im;
catch
    ID = ones(size(pt));
end
try
[im,imc] = imread([TAGfl TAG]); try im = ind2rgb(im,imc); catch; end
TAG = im;
catch
    TAG = ones(size(pt));
end

try
    fileloc2 = fileloc(1:regexp(fileloc,'tag_data')-1);
    [~,~,txt] = xlsread([fileloc2 'TAG GUIDE.xlsx']);
catch
    try
        DD = dir(fileloc2); DD = {DD.name}';
        TG = [fileloc2 DD{~cellfun(@isempty,cellfun(@(x) regexpi(x,'TAG GUIDE'),DD,'uniformoutput',false))}];
        [~,~,txt] = xlsread(TG);
        disp(['Read ' TG ' as TAG GUIDE']);
    catch
        [filename2,fileloc2] = uigetfile('*.*','Select TAG GUIDE');
        [~,~,txt] = xlsread([fileloc2 filename2]);
    end
end
row = find(cellfun(@(x) strcmp(x,whaleName),txt(:,1)));
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Animal'),txt(3,:),'uniformoutput',false)));
name = txt{row,col};
if strcmp(name,'U'); name = 'Unident.'; end
col = find(cellfun(@any,cellfun(@(x) strcmp(x,'Notes'),txt(3,:),'uniformoutput',false)));
note = txt{row,col};
try
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'First seen'),txt(3,:),'uniformoutput',false)));
seen = txt{row,col}; if ~isempty(seen) && ~any(isnan(seen)); note = ['First seen ' num2str(seen(1:min(length(seen),4))) ', ' note]; end
catch
   warning('No "first seen" column, item not included'); 
end
try
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Genetic Sex'),txt(3,:),'uniformoutput',false)));
sex = txt{row,col}; if ~isempty(sex) && ~any(isnan(sex)); note = ['Sex: ' sex ', ' note]; end
catch
   warning('No "genetic sex" column, item not included'); 
end
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Study_Area'),txt(3,:),'uniformoutput',false)));
place = txt{row,col};
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Total Data Time'),txt(3,:),'uniformoutput',false)));
datatime = txt{row,col};
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Total Video Time'),txt(3,:),'uniformoutput',false)));
vidtime = txt{row,col};
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Drone'),txt(3,:),'uniformoutput',false)));
whalelength = txt{row,col};

col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Lat_On'),txt(3,:),'uniformoutput',false)));
deplat = txt{row,col};
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Long_On'),txt(3,:),'uniformoutput',false)));
deplong = txt{row,col};
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Tag Type'),txt(3,:),'uniformoutput',false)));
tagtype = txt{row,col};
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Tag #'),txt(3,:),'uniformoutput',false)));
try
tagnum = txt{row,col};
catch
    warning('No tag # column, using 0'); tagnum = 0;
end
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Spec'),txt(3,:),'uniformoutput',false)));
spec = txt{row,col}; switch spec; case 'bw'; genus = 'Balaenoptera'; spec = 'musculus'; case 'bs'; genus = 'Balaenoptera'; spec = 'borealis'; case 'bp'; genus = 'Balaenoptera'; spec = 'physalus'; case 'mn'; genus = 'Megaptera'; spec = 'Novaeangliae'; case 'oo'; genus = 'Orcinus'; spec = 'orca'; case 'er'; genus = 'Escrichtius'; spec = 'robustus'; case 'bb'; genus = 'Balaenoptera'; spec = 'bonaerensis'; case 'be'; genus = 'Balaenoptera'; spec = 'edeni'; case 'rt'; genus = 'Rhincodon'; spec = 'typus'; otherwise; genus = 'Unk'; spec = 'sp'; end
col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'PI Contact'),txt(3,:),'uniformoutput',false)));
PI = txt{row,col};

% s = input('Make metadata file in ATN format? 1 = yes, 2 = no ');
if makemetadata
    try
        filelocATN = fileloc(1:regexp(fileloc,'tag_data')-1);
        copyfile([filelocATN 'ATN_Metadata.xls'],[fileloc 'ATN_Metadata.xls'],'f');
    catch
        try
            copyfile([fileloc2 'ATN_Metadata.xls'],[fileloc 'ATN_Metadata.xls']);
        catch
            try
                loc = mfilename('fullpath');
                iiii = strfind(loc,'CATSMatlabTools');
                loc = loc(1:iiii+15);
                temploc = [loc 'templates//'];
                copyfile([temploc 'ATN_Metadata.xls'],[fileloc 'ATN_Metadata.xls']);
            catch
            [filenameATN,filelocATN] = uigetfile('*.*','Select ATN_Metadata Template');
            copyfile([filelocATN filenameATN],[fileloc 'ATN_Metadata.xls'],'f');
            end
        end
    end
    xlswrite([fileloc 'ATN_Metadata.xls'],{whaleName,'CATS',tagtype,tagnum,deplong,deplat,datestr(datenum(whaleName(3:8),'YYMMDD'),'YYYY-MM-DD'),genus,spec,'2021-01-01','',name,PI,'davecade@alumni.stanford.edu'},'A2:N2');
end

try
    [im,imc] = imread([fileloc tdr]); try im = ind2rgb(im,imc); catch; end
    tdr = im;
    [im,imc] = imread([fileloc prhi]); try im = ind2rgb(im,imc); catch; end
    prhi = im;
catch %from vert2sxsprh
    load([fileloc prhf]);
    pitch = load([fileloc prhf],'pitch'); pitch = pitch.pitch;
    load([fileloc prhf],'roll');
    DN2= 1:length(DN); %changes the DN graphs back to regular indices the lazy way (instead of deleting the DN wrapper)
    clear pat;
    fig = figure(2); clf; %set(fig2,'color',[.8 .8 .8]);
    boxsize = 300; %bo
    ysize = 260;%215; % vertical pixel size of graph
    gap = 3; %pixels between video and graph
    lowrat = .75; % the
    set(fig,'windowstyle','normal');
    a = find(p>1,1,'first');
%     a = I29 - 400; % I set the I29 to be the index of the first video I wanted to see (i.e. I29 = find(T.DV(:,4)>=16&T.DV(:,5)>=24&T.DV(:,6)>=04,1,'first')), then it only graphs from there on out
    b = find(p>1,1,'last');
    if ~exist('viddeploy','var'); viddeploy = []; end; if ~exist('vidDN','var'); vidDN = []; end; if ~exist('vidDurs','var'); vidDurs = []; end
    nonpat = nan(size(viddeploy));
    labs = num2str(viddeploy); %labs = labs(:,5:6);
    if ~isempty(viddeploy); combinesmall = true; else combinesmall = false; end
    if combinesmall; combos = getcombos(vidDN,vidDurs,viddeploy,2700,1200); else combos = num2cell(viddeploy); end
    for i = 1:length(combos); for j = 2:length(combos{i}); labs(viddeploy==combos{i}(j),:) = ' '; end; end
    for i = 1:length(viddeploy)
        starttime = max(0,round((DN(find(tagon,1)-15*fs)-vidDN(viddeploy(i)))*24*60*60)); 
        et = min(vidDurs(viddeploy(i)),(DN(end)-vidDN(viddeploy(i)))*24*60*60);
        [~,a1] = min(abs(DN-(vidDN(viddeploy(i))+starttime/24/60/60)));
%         a1 = round(a1+starttime(i)*fs);
        plot(DN2([a1 a1]),[-10 1000],'r','linewidth',2); hold on;
        ys = get(gca,'ylim');
        b1 = min(length(p),a1 + min(round((et-starttime)*fs),round((vidDurs(viddeploy(i))-starttime) * fs)) -1);
        nonpat(i) = rectangle('position',[DN2(a1) ys(1) DN2(b1)-DN2(a1) ys(2)],'facecolor',[255 192 203]/255);
        oi = get(gca,'children');
        set(gca,'children',[oi(2:end); oi(1)]);
        text(DN2(a1),0,labs(i,:),'verticalalignment','bottom','fontsize',10,'color','r')
    end
    plot(DN2(a:b),p(a:b),'b','linewidth',1.5);
    scrn = get(0,'screensize');
    set(fig,'units','pixels','Position',[1,50,scrn(3),round(ysize*lowrat)]); % assumes 2560 frame size, which it should be
    set(gca,'xlim',[DN2(a) DN2(b)],'ylim',[0 max(p(a:b))],'fontsize',16);
    oi = datestr(DN(get(gca,'xtick')),'HH:MM:SS');
    set(gca,'xticklabel',oi,'ydir','rev');
    xlabel('Local Time'); ylabel('Depth');
    set(gca,'units','points');
    pos = get(gca,'position'); pos(1) = 50;
    set(gca,'position',pos); set(gca,'units','normalized'); pos = get(gca,'position');
    set(gca,'position',[pos(1:2) 1-pos(1)-.005 pos(4)]);
    xs = get(gca,'xlim');
    ys = get(gca,'ylim');
    text(xs(1),ys(2),'Full Deployment Record','fontsize',16,'verticalalignment','bottom');
    try saveas(fig,[ql whaleName ' TDR' '.bmp']); catch; saveas(fig,[fileloc whaleName ' TDR' '.bmp']); end
    try [im,imc] = imread([ql whaleName ' TDR' '.bmp']); catch; [im,imc] = imread([fileloc whaleName ' TDR' '.bmp']); end
    try im = ind2rgb(im,imc); catch; end
    tdr = im;
    
    fig5 = figure(5); clf; set(fig5,'color','w','windowstyle','normal');
    set(gca,'units','points');
    pos = get(gca,'position'); pos(1) = 50;
    set(gca,'position',pos); set(gca,'units','normalized'); pos = get(gca,'position');
    set(gca,'position',[pos(1:2) 1-2.3*pos(1) pos(4)]);
    xs = xs(1):xs(end); %xs from above
    [axp5, hp51,hp52] = plotyy(xs,pitch(xs)*180/pi,xs,head(xs)*180/pi);
    hold(axp5(1),'on'); hold(axp5(2),'on');
    hp53 = plot(xs,roll(xs)/pi*180,'r.','parent',axp5(2),'markersize',3);
    hp54 = plot(-10,5,'r-','linewidth',3); hp55 = plot(-10,5,'b-','linewidth',3); hp56 = plot(-10,5,'k--','linewidth',2);
    set(hp51,'color','g','linewidth',2); set(hp52,'color','b','marker','.','markersize',3,'linestyle','none');
    set(axp5,'xlim',[xs(1) xs(end)]);
    set(axp5(1),'ylim',[-90 90],'ytick',-90:30:90,'box','off','ycolor','g'); set(axp5(2),'ycolor','k','ylim',[-181 181],'ytick',-120:60:120);
    ylabel('degrees (pitch)');
    %         if sum([headgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))>60/180*pi; pitchgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))>30/180*pi]) > sum([headgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))<-60/180*pi; pitchgy(xs(1):round(.08*(xs(2)-xs(1))+xs(1)))<-30/180*pi])
    %             legloc = 'southwest'; else legloc = 'northwest'; end
    %         leg = legend(hp6,'Gyros','location',legloc,'orientation','horizontal');
    set(get(axp5(2),'ylabel'),'string',{'degrees' '(roll and heading)'},'fontsize',16);
    set(get(axp5(1),'ylabel'),'fontsize',16);
    set(fig5,'units','pixels','Position',[1,70,scrn(3),round(ysize*lowrat)]);%fix this
    set(gca,'units','points'); % these four lines repeated to ensure it stretches out right
    pos = get(gca,'position'); pos(1) = 50;
    set(gca,'position',pos); set(gca,'units','normalized'); pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2)+.025 1-2.5*pos(1)-.005 pos(4)]);
    set(axp5,'xtick',round(xs(1)+(xs(end)-xs(1))/10:(xs(end)-xs(1))/10:xs(end)));
    oi = datestr(DN(get(axp5(2),'xtick')),'HH:MM:SS');
    set(axp5,'xticklabel',oi,'fontsize',16);
    
%     plot(axp5(1),xs,skipgaps(pitchgy(xs)*180/pi,100),'g--','linewidth',2);
%     plot(axp5(2),xs,skipgaps(rollgy(xs)*180/pi,100),'r--','linewidth',2);
%     plot(axp5(2),xs,skipgaps(headgy(xs)*180/pi,100),'b--','linewidth',2);
    ys = get(axp5(1),'ylim');
    xs = get(axp5(1),'xlim');
    ltext = text(xs(1)-(xs(2)-xs(1))*.028,ys(1)-(ys(2)-ys(1))/40,'Local Time: ','parent',axp5(1),'verticalalignment','top','fontname',get(axp5(1),'fontname'),'fontsize',get(axp5(1),'fontsize'),'horizontalalignment','left');
    try saveas(fig5,[ql whaleName ' prh' '.bmp']); catch; saveas(fig5,[fileloc whaleName ' prh' '.bmp']); end
    try [im,imc] = imread([ql whaleName ' prh' '.bmp']); catch; [im,imc] = imread([fileloc whaleName ' prh' '.bmp']);  end
    try im = ind2rgb(im,imc); catch; end
    prhi = im;
end
%
clear ql;
ql.cdata = ones(1000,1920,3);
%name
tx = 400; ty = 100;
t = insertText(ones(ty,tx,3),[0 0],whaleName,'fontsize',50,'boxcolor','white');%,'font','freesansbold','fontstyle','bold','boxcolor','white');

tx = 220; ty = 60;
[t, tmap] = rgb2ind(t,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(1:ty,1:tx,:) = t;

%name2
tx = 400; ty = 100;
t = insertText(ones(ty,tx,3),[0 0],['"' name '"'],'fontsize',50,'boxcolor','white');%,'font','freesansbold','fontstyle','bold','boxcolor','white');

tx = 220; ty = 55;
[t, tmap] = rgb2ind(t,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(41:40+ty,1:tx,:) = t;

%place
tx = 400; ty = 100;
t = insertText(ones(ty,tx,3),[0 0],place,'fontsize',50,'boxcolor','white');%,'font','freesansbold','fontstyle','bold','boxcolor','white');

tx = 220; ty = 55;
[t, tmap] = rgb2ind(t,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(76:75+ty,1:tx,:) = t;

%datatime
HH = floor(datatime*24);
MM = floor((datatime*24-HH)*60);
SS = floor(((datatime*24-HH)*60-MM)*60);
oiD = ['data: ' num2str(HH) ':' sprintf('%02.f', MM) ':' sprintf('%02.f', SS)];
tx = 400; ty = 100;
t = insertText(ones(ty,tx,3),[0 0],oiD,'fontsize',50,'boxcolor','white');%,'font','freesansbold','fontstyle','bold','boxcolor','white');

tx = 220; ty = 55;
[t, tmap] = rgb2ind(t,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(111:110+ty,1:tx,:) = t;

%vidtime
if isfloat(vidtime)
    HH = floor(vidtime*24);
    MM = floor((vidtime*24-HH)*60);
    SS = floor(((vidtime*24-HH)*60-MM)*60);
else HH = 0; MM = 0; SS = 0;
end
oiD = ['video: ' num2str(HH) ':' sprintf('%02.f', MM) ':' sprintf('%02.f', SS)];
tx = 400; ty = 100;
t = insertText(ones(ty,tx,3),[0 0],oiD,'fontsize',50,'boxcolor','white');%,'font','freesansbold','fontstyle','bold','boxcolor','white');

tx = 220; ty = 55;
[t, tmap] = rgb2ind(t,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(146:145+ty,1:tx,:) = t;

%whalelength
if ischar(whalelength) && ~strcmp(whalelength,'n/a')
    oiD = ['length: ' whalelength];
    tx = 400; ty = 100;
    t = insertText(ones(ty,tx,3),[0 0],oiD,'fontsize',50,'boxcolor','white');%,'font','freesansbold','fontstyle','bold','boxcolor','white');
    
    tx = 220; ty = 55;
    [t, tmap] = rgb2ind(t,255);
    [t, tmap] = imresize(t,tmap,[ty tx]);
    t = ind2rgb(t,tmap);
    ql.cdata(181:180+ty,1:tx,:) = t;
end

% TDR
tx = 1200; ty = 140;
[t, tmap] = rgb2ind(tdr,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(1:ty,230:230+tx-1,:) = t;

% PRH
tx = 1200; ty = 120;
[t, tmap] = rgb2ind(prhi,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(131:131+ty-1,230:230+tx-1,:) = t;

% TAG
ty = 280; % tx = 1920 - tx - 230 - 5; 
tx = round(ty/size(TAG,1)*size(TAG,2));%round(1000/1920 * );
[t, tmap] = rgb2ind(TAG,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(1:ty,1431:1430+tx,:) = t;

% ID
ty = 280; % tx = 1920 - tx - 230 - 5; 
tx = round(ty/size(ID,1)*size(ID,2));
[t, tmap] = rgb2ind(ID,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(291:290+ty,1431:1430+tx,:) = t;
ql.cdata(:,1921:end,:) = [];

% kml
ty = 1000-2*280-20; tx = 1920 - 1200 - 230; 
[t, tmap] = rgb2ind(kml,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(1000-ty+1:1000,1920-tx+1:1920,:) = t;

%map
ty = 1000-2*280-20; tx = 1920 - 1200 - 230; 
map2 = map(round(45/504*size(map,1)):round(490/504*size(map,1)),1:round(635/672*size(map,2)),:);
[t, tmap] = rgb2ind(map2,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(1000-ty+1:1000,1200+222-tx+1:1200+222,:) = t;

%ptrack
ty = 1000-2*280-20; tx = 1920 - 1200 - 230; 
if ~drone pt2 = pt(round(65/976*size(pt,1)):round(940/976*size(pt,1)),round(120/1902*size(pt,2)):round(1750/1902*size(pt,2)),:);
else pt2 = DRONE(round(65/976*size(DRONE,1)):round(940/976*size(DRONE,1)),round(120/1902*size(DRONE,2)):round(1750/1902*size(DRONE,2)),:);
end
[t, tmap] = rgb2ind(pt2,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(1000-ty+1:1000,1:tx,:) = t;


%gtrack
ty = 1000-2*280-20; tx = 1920 - 1200 - 230; 
gt2 = gt(round(65/976*size(gt,1)):round(940/976*size(gt,1)),round(120/1902*size(gt,2)):round(1750/1902*size(gt,2)),:);
[t, tmap] = rgb2ind(gt2,255);
[t, tmap] = imresize(t,tmap,[ty tx]);
t = ind2rgb(t,tmap);
ql.cdata(1000-ty+1:1000,933-tx+1:933,:) = t;

% cam
if nocam
    tx = 0;
else
    ty = 580-255; % tx = 1920 - tx - 230 - 5;
    tx = round(ty/size(cam,1)*size(cam,2));%round(1000/1920 * );
    [t, tmap] = rgb2ind(cam,255);
    [t, tmap] = imresize(t,tmap,[ty tx]);
    t = ind2rgb(t,tmap);
    ql.cdata(256:255+ty,1200+222-tx+1:1200+222,:) = t;
end
%note
w = 1200+222-tx;
tx = 1500; ty = 45;
try ll = 0; sp = [regexp(note,' ') length(note)];
catch warning('No note detected for display in quicklook, rerun in tag guide or enter here (use single quotes around statement): ');
    note = input('note? ');
    ll = 0; sp = [regexp(note,' ') length(note)];
end
sy = 256;
%
if nocam
    w = 870;
end
while ll<length(note)
    tx = 1500; ty = 45;
    oi = note(ll+1:sp(find(sp<ll+w/10,1,'last')));
    ll = length(oi)+ll;
    t = insertText(ones(ty,tx,3),[0 0],oi,'fontsize',35,'boxcolor','white');%,'font','freesansbold','fontstyle','bold','boxcolor','white');
    tx = 870; ty = 25;
    [t, tmap] = rgb2ind(t,255);
    [t, tmap] = imresize(t,tmap,[ty tx]);
    t = ind2rgb(t,tmap);
    t(:,w+1:end,:) = [];
    ql.cdata(sy+1:sy+ty,1:min(w,tx),:) = t;
%     image(ql)
    sy = sy+ty;
    if sy+ty>580; break; end
end
%
figure; 

set(gcf,'Units','pixels', 'outerposition',[0 0 1920 1000], 'menubar','none');
set(gca,'Visible', 'Off', 'Position',[0 0 1 1],'xlim',[0 1920],'ylim',[0 1000]);
 set(gca,'ydir','rev');

image(ql);

imwrite(ql.cdata,[fileloc '_' whaleName 'Quicklook.jpg']);

