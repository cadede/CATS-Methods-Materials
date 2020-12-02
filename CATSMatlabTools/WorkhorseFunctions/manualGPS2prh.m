function [GPS,GPSerr] = manualGPS2prh (prhloc, prhfile)
% Choose a prhfile and a GPS file (with date/time column, lat, long)

if nargin<2
    try cd(prhloc);
        prhfile = uigetfile('*prh.mat','Choose prh file');
    catch
    end
end

if nargin<1
    [prhfile,prhloc] = uigetfile('*prh.mat','Choose prh file');
end

D = dir(prhloc);
files = {D.name};
I = ~cellfun(@isempty,cellfun(@(x) regexp(x,'GPS'),files,'uniformoutput',false)) & ~cellfun(@isempty,cellfun(@(x) regexp(x,'.xls'),files,'uniformoutput',false));
if sum(I)>1 || sum(I) == 0
    cd(prhloc);
    [GPSfile,GPSloc] = uigetfile('*.*','Choose file with GPS hits, press cancel to only use deploy and recover locations from tag guide');
else
    GPSfile = files{I}; GPSloc = prhloc;
    disp(['Using GPS file: ' GPSfile]);
end
disp('Values should be in local time,GPS should be in decimal degrees');

load([prhloc prhfile]);
% GPS = GPS(1,:);
%
if INFO.tagnum>47 && ~ismember(INFO.tagnum,[70 71]) % GPS in back so get floating positions as soon as tag pops off (but not for fastGPS tags)
    II = (1:length(p))';
    try    tagoffGPS = GPS(find(~isnan(GPS(:,1))&II>find(tagon,1,'last')&II<find(tagon,1,'last')+20*fs*60,30),:);
        if ~isempty(tagoffGPS) && sum(isnan(tagoffGPS(:,1)))~=length(tagoffGPS(:,1))
            figure(105); clf; plot(tagoffGPS(:,2),tagoffGPS(:,1),'bs');
            strs = cellstr(num2str((1:size(tagoffGPS,1))'));
            text(tagoffGPS(:,2),tagoffGPS(:,1),strs);
            title('To use one of these GPS points (from after the tag pops off and is floating at the surface), click the point (numbered in chronological order), else press any other key to use the recover time from the tag guide');
            [x,y,button] = ginput(1);
            if button == 1
                dist = sqrt(sum((tagoffGPS-repmat([y,x],length(tagoffGPS),1)).^2,2));
                [~,b] = min(dist);
                tagoffGPS = tagoffGPS(b,:);
            else
                clear tagoffGPS;
            end
        else
            clear tagoffGPS
        end
    catch; clear tagoffGPS;
    end
end



if ~ischar(GPSfile);
    hits = cell(3,3);
    hits(1,:) = [{'Time'} {'Lat'} {'Long'}];
    [~,~,txt] = xlsread([prhloc(1:regexp(prhloc,'CATS\')+4) 'TAG GUIDE.xlsx']);
    row = find(cellfun(@(x) strcmp(x,getWhaleID(prhfile)),txt(:,1)));
    col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Tag_On'),txt(3,:),'uniformoutput',false)));
    hits(2,1) = txt(row,col);
    col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Recover_Time'),txt(3,:),'uniformoutput',false)));
    hits(3,1) = txt(row,col);
    col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Lat_On'),txt(3,:),'uniformoutput',false)));
    hits(2,2) = txt(row,col);
    col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Long_On'),txt(3,:),'uniformoutput',false)));
    hits(2,3) = txt(row,col);
    col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Recover_Lat'),txt(3,:),'uniformoutput',false)));
    hits(3,2) = txt(row,col);
    col = find(~cellfun(@isempty,cellfun(@(x) strfind(x,'Recover_Long'),txt(3,:),'uniformoutput',false)));
    hits(3,3) = txt(row,col);
    try
        if any(cellfun(@(x) any(isnan(x)),hits(3,:))) % if there is not a recovery location
            disp('No Recovery Time or Location.  Do you have a tagoff location (e.g. high quality Argos close to tagoff time?)');
            t = input('1 = yes 2 = no  ');
            if t == 1;
                tagoffGPS(1) = input('Lat?  ');
                tagoffGPS(2) = input('Long?  ');
            end
        elseif (datenum(hits(end,1))-DN(find(tagon,1,'last'))) > 1/24 && INFO.tagnum<=47; % if it was more than an hour between tag off and tag recovery
            disp(['Tag recovery was ' num2str((datenum(hits(end,1))-DN(find(tagon,1,'last')))*24) ' hours after tag off, do you have a better tag off location (e.g. high quality Argos close to tagoff time?).  ']);
            t = input('1 = yes 2 = no  ');
            if t == 1;
                tagoffGPS(1) = input('Lat?  ');
                tagoffGPS(2) = input('Long?  ');
            end
        end
    catch
    end
    if exist('tagoffGPS','var')
       hits{3,1} = datestr(DN(find(tagon,1,'last')-fs),'mm/dd/yyyy HH:MM:SS');  hits{3,2} = tagoffGPS(1);  hits{3,3} = tagoffGPS(2); 
    end
else
    [~,~,hits] = xlsread([GPSloc GPSfile]);
end
timecol = find(~cellfun(@isempty,cellfun(@(x) regexpi(x,'Time'),hits(1,:),'uniformoutput',false)),1);
latcol = find(~cellfun(@isempty,cellfun(@(x) regexpi(x,'Lat'),hits(1,:),'uniformoutput',false)),1);
longcol =  find(~cellfun(@isempty,cellfun(@(x) regexpi(x,'Long'),hits(1,:),'uniformoutput',false)),1);
hits = hits(2:end,:);
%

try gDN = datenum(hits(:,timecol));
catch gDN = datenum(vertcat(hits{:,timecol}));
end
lat = vertcat(hits{:,latcol}); long = vertcat(hits{:,longcol});
nums = 1:length(lat);
button = 1;
while button~=113
    if isempty(gDN); button = 113; continue; end
    figure(42); clf;
    s1 = subplot(211);
    plot(DN,p,gDN,0,'rs'); set(gca,'ydir','rev');
    xlim([min(gDN)-5/60/24 max(gDN)+6/60/24]);
    ylim([-5 20]);
    for i = 1:length(gDN)
        text(gDN(i),0,num2str(nums(i)),'verticalalignment','bottom','horizontalalignment','center','fontsize',8);
    end
    set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
    s2 = subplot(223);
    [x,y,uz] = deg2utm(lat,long); %need utmzone
%     xp = utm2m(x,y,uz);
    if size(GPS,1)>2
        gI = ~isnan(GPS(:,1))&tagon;
        oGPS = GPS(gI,:);
        try [xo,yo,uzo] = deg2utm(oGPS(:,1),oGPS(:,2)); xo = utm2m(xo,yo,uzo); catch; xo = []; yo = []; uzo = []; end
        p1 = plot(xo,yo,'bx-'); hold on;
    else
        p1 = plot(0,0,'bx'); hold on;
        xo = []; yo = []; uzo = [];
    end
    xp = utm2m([x; xo],[y; yo],[uz;uzo]);
    xp = xp(1:length(x));
    p2 = plot(xp,y,'rx-');
    for i = 1:length(gDN)
        text(xp(i),y(i),num2str(nums(i)),'verticalalignment','bottom','horizontalalignment','center','fontsize',6);
    end
    legend([p1,p2],'Tag','Other','location','northeastoutside');
        xlim([min([xp;xo])-100 max([xp;xo])+100]); ylim([min([y;yo])-100 max([y;yo])+100]); axis equal;
    
    title(s1,'Can Zoom in and out on dive profile to examine points, press enter when ready to adjust points.'); 
    z1 = zoom(s1); set(z1,'enable','on'); set(z1,'Motion','horizontal');
    z2 = zoom(s2); set(z2,'enable','on','Motion','both');
    pause;
    set(s1,'xticklabel',datestr(get(s1,'xtick'),'HH:MM:SS'));
    title(s1,{'Now click next to a new GPS point on the dive profile to remove it, or press m next to a point to move it in time, then click on the point to move it to'; 'Press q to finish'});
    set(gcf,'CurrentAxes',s1)
    [x1,~,button] = ginput(1);
    [~,b] = min(abs(gDN-x1));
    if button == 109
        [x2,~] = ginput(1);
        gDN(b) = x2;
    elseif button~=113
        gDN(b) = [];
        lat(b) = [];
        long(b) = [];
        nums(b) = [];
    end
end
if size(GPS,1)>2
    button = 3;
    while isempty(button) || button~=113
        figure(42); clf;
        s1 = subplot(211);
        if ~isempty(gDN); plot(DN,p,gDN,0,'rs'); set(gca,'ydir','rev');
        xlim([min(gDN)-5/60/24 max(gDN)+6/60/24]);
        ylim([-5 20]);
        end
        for i = 1:length(gDN)
            text(gDN(i),0,num2str(nums(i)),'verticalalignment','bottom','horizontalalignment','center','fontsize',8);
        end
        set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
        s2 = subplot(223);
        [x,y,uz] = deg2utm(lat,long);
%         xp = utm2m(x,y,uz);
        if size(GPS,1)>1
            gI = ~isnan(GPS(:,1))&tagon;
            oGPS = GPS(gI,:);
            gI = find(gI);
            try [xo,yo,uzo] = deg2utm(oGPS(:,1),oGPS(:,2)); xo = utm2m(xo,yo,uzo); catch; xo = []; yo = []; uzo = []; end
            gspeeds = sqrt(sum(diff([xo yo]).^2,2))./(diff(DN(gI))*24*60*60);
            bigspeeds = find([gspeeds;0]>6&[0;gspeeds]>6); % find the spots where speed must average >6 m/s to get to and from the point
            p0 = plot(xo(bigspeeds),yo(bigspeeds),'gs','markerfacecolor','g'); hold on;
            p1 = plot(xo,yo,'bx-'); 
        else
            xo = []; yo = []; uzo = [];
            p1 = plot(0,0,'bx-'); hold on;
        end
        xp = utm2m([x; xo],[y; yo],[uz;uzo]);
        xp = xp(1:length(x));
        p2 = plot(xp,y,'rx-');
        for i = 1:length(gDN)
            text(xp(i),y(i),num2str(nums(i)),'verticalalignment','bottom','horizontalalignment','center','fontsize',6);
        end
        legend([p1,p2,p0],{'Tag','Other','Speed > 6 m/s'},'location','northeastoutside');
        if isempty(button) || button == 3
            xlim([min([xp;xo])-100 max([xp;xo])+100]); ylim([min([y;yo])-100 max([y;yo])+100]); axis equal;
            oldxs = get(gca,'xlim');
            oldys = get(gca,'ylim');
            title(s1,{'Goal is to erase any bad tag hits. on the lower left plot'; 'Can Zoom in and out, press enter when ready to adjust the points.'});
            z1 = zoom(s1); set(z1,'enable','on'); set(z1,'Motion','horizontal');
            z2 = zoom(s2); set(z2,'enable','on','Motion','both');
            pause;
        else; xlim(oldxs); ylim(oldys);
        end
        set(s1,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
        title(s1,{'Now click next to a tag GPS point on the map to remove it'; 'Press "del" to delete all green points'; 'press enter to change the zoom' ;  'or press q to finish'});
        set(gcf,'CurrentAxes',s2)
        [x1,y1,button] = ginput(1);
        if button == 127
            for i = 1:length(bigspeeds)
                x1 = xo(bigspeeds(i)); y1 = yo(bigspeeds(i));
                [~,b] = min(arrayfun(@(x,y) sqrt((x1-x)^2+(y1-y)^2),xo,yo));
                GPS(gI(b),:) = nan;
                GPSerr(gI(b),:) = nan;
                oldxs = get(gca,'xlim');
                oldys = get(gca,'ylim');
            end
        elseif isempty(button)
            continue;
        elseif button~=113
            [~,b] = min(arrayfun(@(x,y) sqrt((x1-x)^2+(y1-y)^2),xo,yo));
            GPS(gI(b),:) = nan;
            GPSerr(gI(b),:) = nan;
            oldxs = get(gca,'xlim');
            oldys = get(gca,'ylim');
        end
    end
end
oGPS = GPS(1,:); if size(oGPS,2) == 1; oGPS = GPS(1:2)'; end
if size(GPS,1)<=2
    GPS = nan(length(p),2);
    GPSerr = nan(length(p),3);
end
GPS(1,:) = oGPS;
for i = 1:length(gDN)
    [~,b] = min(abs(DN-gDN(i)));
%     if sum(isnan(GPS(max(3,b-60*fs):min(length(p),b+60*fs),1))) == length(GPS(max(3,b-60*fs):min(length(p),b+60*fs),1)) %if there are no tag hits within 60 seconds
        GPS(b,:) = [lat(i) long(i)];
%     end
end
save([prhloc prhfile],'GPS','GPSerr','-append');
disp(['Number of total GPS hits: ' num2str(sum(~isnan(GPS(:,1)))-1)]);




