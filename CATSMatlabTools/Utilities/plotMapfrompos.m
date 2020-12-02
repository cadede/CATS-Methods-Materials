function [fig,ax] = plotMapfrompos(GPS,DN,tagon,p,fs,mapfileloc,edges)
   global fileloc
%plots a coastline map and the GPS points in GPS up to the boundary of
%edges.  Default is 10 km.
if nargin<7; edges = 10; end
%
try
%     curdir = pwd;
%     oi = strfind(curdir,'MATLAB');
%     [D2,F2] = subdir(curdir(1:oi+5));
%     oi = find(cellfun(@isempty,cellfun(@(x) strfind(x,'CATSMatlabTools\Utilities\map files'),D2,'uniformoutput',false))==0);
%     [D2,F2] = mapfileloc;
    F2 = dir(mapfileloc); F2 = {F2.name};
    Fnum = find(cellfun(@(y) y==1,cellfun(@(x) strcmp(x,'GSHHS_f_L1.shp'),F2,'uniformoutput',false)));
    FnumB = find(cellfun(@(y) y==1,cellfun(@(x) strcmp(x,'bathy_l.shp'),F2,'uniformoutput',false)));
%     DIR = [D2{oi} '\']; 
    DIR = mapfileloc;
%     filename =  F2{oi}{Fnum}; bathy = F2{oi}{FnumB};
    filename = F2{Fnum}; bathy = F2{FnumB};
    river = cell(9,1);
    for k = 1:9
        Fnum(k) = find(cellfun(@(y) y==1,cellfun(@(x) strcmp(x,['WDBII_river_f_L0' num2str(k) '.shp']),F2,'uniformoutput',false))); 
        river{k} = F2{Fnum(k)};
    end
catch
    error('likely could not find mapfiles, check mapfileloc variable');
end
    
% catch
%     [~,DIR] = uigetfile('*.shp','Find GSHHS shape file (in CATSMatlabTools\Utilities\map files), other files should be in the same folder');
%     cd (DIR)
% %     [bathy,DIR] = uigetfile('*.shp','Find bathy_l shape file (in CATStools\Scripts\Tools), other files should be in the same folder');
%     curdir = pwd;
%     oi = strfind(curdir,'MATLAB');
%     [D2,F2] = subdir(curdir(1:oi+5));
%     oi = find(cellfun(@isempty,cellfun(@(x) strfind(x,'CATSMatlabTools\Utilities\map files'),D2,'uniformoutput',false))==0);
%     Fnum = find(cellfun(@(y) y==1,cellfun(@(x) strcmp(x,'GSHHS_f_L1.shp'),F2{oi},'uniformoutput',false)));
%     FnumB = find(cellfun(@(y) y==1,cellfun(@(x) strcmp(x,'bathy_l.shp'),F2{oi},'uniformoutput',false)));
%     DIR = [D2{oi} '\']; filename =  F2{oi}{Fnum}; bathy = F2{oi}{FnumB};
%     river = cell(9,1);
%     for k = 1:9
%         Fnum(k) = find(cellfun(@(y) y==1,cellfun(@(x) strcmp(x,['WDBII_river_f_L0' num2str(k) '.shp']),F2{oi},'uniformoutput',false))); 
%         river{k} = F2{oi}{Fnum(k)};
%     end
% end
%
mins = min(GPS(tagon,:));
maxs = max(GPS(tagon,:));
if any(isnan(mins))||any(mins-maxs==0)
    mins = GPS(find(~isnan(GPS(:,1))&DN<DN(find(tagon,1)),2,'last'),:);
    maxs = max(mins); mins = min(mins);
end
if any(isnan(mins))||any(mins-maxs==0)
     mins = GPS(find(~isnan(GPS(:,1))&DN>DN(find(tagon,1,'last')),2,'first'),:);
     maxs = max(mins); mins = min(mins);
end
DIST = distdim(edges * 1000,'m','degrees');
longD = distance(mins(1),mins(2),mins(1),maxs(2));
latD = distance(mins(1),mins(2),maxs(1),mins(2));
edgeDist = (maxs-mins).*(DIST./[latD longD]); %the degrees of lat and long that make the "edges" distance
maxs2 = maxs+(maxs-mins).*(DIST./[latD longD]);
mins = mins-(maxs-mins).*(DIST./[latD longD]);
maxs = maxs2;
S = shaperead([DIR filename],'BoundingBox',fliplr([mins;maxs]));
B = shaperead([DIR bathy],'BoundingBox',fliplr([mins;maxs]));
fig = figure(642); clf;
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 

h = nan(6,1);
mapshow(S);
set(gca,'xlim',[mins(2) maxs(2)],'ylim',[mins(1) maxs(1)]);
hold on;
I = ~isnan(GPS(:,1));
h(1) = plot(GPS(I,2),GPS(I,1),'ko','markersize',3,'markerfacecolor','k');
t1 = DN(find(tagon,1)); t2 = DN(find(tagon,1,'last'));
c = colormap;
% plot(GPS(I,2),GPS(I,1),'s');
% surface([GPS(I,2)';GPS(I,2)'],[GPS(I,1)';GPS(I,1)'],zeros(size(GPS(I,:)))',[DN(I)';DN(I)],'facecol','no','edgecol','interp','linew',2);
ts= t1:(t2-t1)/size(colormap,1):t2;
for i = 1:length(ts)-1
    ii = find(DN>=ts(i)&DN<=ts(i+1)&~isnan(GPS(:,1)));
    oi = plot(GPS(ii,2),GPS(ii,1),'o','color',c(i,:),'markerfacecolor',c(i,:));
    if ~isempty(oi); oi = oi(1); h(2) = oi; end
end
    
for k = 1:length(B)
    switch B(k).DEPTH
        case 1; continue;
        case 200
            h(4) = plot(B(k).X,B(k).Y,'color',[190 190 190]/255);
        case 500
            h(5) = plot(B(k).X,B(k).Y,'color',[105 105 105]/255);
        case 2500
            h(6) = plot(B(k).X,B(k).Y,'k');
    end
end
for k = 1:length(river)
    R = shaperead([DIR river{k}],'BoundingBox',fliplr([mins;maxs]));
    oi = mapshow(R);
    if ~isempty(oi); oi = get(oi,'children'); h(3) = oi(1); end
end

h(1) = plot(190,190,'ko','markersize',3,'markerfacecolor','k');
h(2) = plot(190,190,'o','color','r','markerfacecolor','r');
h(3) = plot([189 190],[189 190],'b-');
h(4) = plot(190,190,'color',[190 190 190]/255);
h(5) = plot(190,190,'color',[105 105 105]/255);
h(6) = plot(190,190,'color','k');

cb = colorbar;
% h(3) = nan;
legend(h,'Tag Off','Tag On','River','200 m','500 m','2500 m');
ax = gca;
%
ylabel('Latitude','fontsize',16);
xlabel('Longitude','fontsize',16);
i = imread([DIR 'scalebar5-50.png']);
scaled = (1.9+5.35*5)/(5.35*5);
xs = [mins(2) mins(2)+5/edges*edgeDist(2)*scaled];
ys = [mins(1) mins(1) + diff(xs)*103/740];
im = image(xs,ys,flipud(i));
text(xs(1) + diff(xs)/2,mins(1) + diff(xs)*103/740,'5 km','fontsize',16,'verticalalignment','bottom','horizontalalignment','center');
ys = get(cb,'ytick'); set(cb,'ytick',min(ys):(max(ys)-min(ys))/10:max(ys)); %sets 10 y ticks
ylabs = datestr(t1:(t2-t1)/10:t2,'HH:MM');
set(cb,'yticklabel',ylabs)
% plot([mins(2) mins(2)+5/edges*edgeDist(2)],ys(2)

            surfs = findsurfacings(p,fs,tagon,60);
            xs = get(ax,'xlim'); ys = get(ax,'ylim');
            text(xs(1)+diff(xs)/20,ys(2)-diff(ys)/20,{['Number of surfacings: ' num2str(length(surfs))]; ['Number of GPS hits on surfacings: ' num2str(size(GPS(tagon&~isnan(GPS(:,1))),1))]},'parent',ax);

set(fig,'units','normalized','outerposition',[0 0 1 1]);
 

