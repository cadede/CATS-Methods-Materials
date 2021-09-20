function [h,barh,barlabsh] = plotspeed(ax,x,speed,colors)

% plots the composite speed column in a speed table using colors indicating
% where each value came from
% [h,barh,barlabsh] = plotspeed(ax,x,speed,colors)
% ax (optional)- if indicated, plots speed in the indicated axis, if not
% indicated "plotspeed(x,speed...)", ax = gca;
% x- the x values corresponding to the speed you want to plot
% speed- a speed table with at least the speed.comp and speed.type variables
% colors- a 4x3 vector of colors.  Optional to include (default is black, red, blue, green
% h - handles for each of the plotted lines
% barh - handle for the colorbar.
% barlabsh - handles for the text labels in the colorbar

if ~isempty(ax) && ~all(ishandle(ax))
    if nargin == 3; colors = speed;  end
    speed = x;
    x = ax;
    ax = gca;
end
if ~exist('colors','var'); colors = [0 0 0; 1 0 0; [255 140 0]/255; 0 1 0]; end
if size(colors,1)~=4 || size(colors,2)~=3; error('"colors" must be 4 x 3'); end
if max(size(x)) ~= size(speed,1) || ~isvector(x); error('"x" must be a vector of x values the same length as the speed table'); end

I = cellfun(@(x) strcmp(x,'IN'),speed.type);
dI = [diff(I); 0]==1; I(dI) = true; %add a point to close the holes when you plot
FNI = cellfun(@(x) strcmp(x,'FN'),speed.type);
dI = [diff(FNI); 0]==1; FNI(dI) = true;
JJI = cellfun(@(x) strcmp(x,'JJ'),speed.type);
dI = [diff(JJI); 0]==1; JJI(dI) = true;
SPI = cellfun(@(x) strcmp(x,'SP'),speed.type);
dI = [diff(SPI); 0]==1; SPI(dI) = true;

segments = cell(4,1); for i = 1:4; segments{i} = nan(size(speed.comp)); end
segments{1}(I) = speed.comp(I);
segments{2}(FNI) = speed.comp(FNI);
segments{3}(JJI) = speed.comp(JJI);
segments{4}(SPI) = speed.comp(SPI);
h = nan(4,1);
oi = get(ax,'nextplot');
for i = 1:4
    h(i) = plot(ax,x,segments{i},'color',colors(i,:)); set(ax,'nextplot','add');
end
set(ax,'nextplot',oi);
colormap(ax,colors);
barh = colorbar('location','south','xticklabel',[]);
xs = get(barh,'xlim'); ticks = xs(1):diff(xs)/8:xs(2); ticks = ticks([2 4 6 8]);
set(barh,'xtick',[]);
pos = get(barh,'position');
oi = get(ax,'position');
pos(1) = pos(1)+0.83*pos(3); pos(3) = 0.2*pos(3); pos(2) = oi(2); %pos(4) = .8*pos(4); 
set(barh,'position',pos);
try barlabsh = text(ticks,mean(get(barh,'ylim'))*ones(4,1),{'Int','Flow','Jig','OCDR'},'parent',barh,'color','w','fontweight','bold','horizontalalignment','center','fontsize',13);
catch
    ax2014b = axes('position',pos,'visible','off');
    barlabsh = text(ticks,mean(get(barh,'ylim'))*ones(4,1),{'Int','Flow','Jig','OCDR'},'parent',ax2014b,'color','w','fontweight','bold','horizontalalignment','center','fontsize',13);
end



