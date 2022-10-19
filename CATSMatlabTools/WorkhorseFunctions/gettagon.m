function tagon = gettagon(Depth,fs,starttime,At)

% get tagon and tagoff times
% id tagon and tagoff times by zooming in and selecting the boundaries of
% time on the whale.
% inputs: Depth variable
%          fs (sampling rate)
%          starttime (matlab datenumber of the starttime- put 0 if unknown)
%          At (another comparable variable.  Set up to be Acceleration, but could use temperature or even depth again just to make the script work if no other data is available)
% output: tagon (an index of values for when the tag was on the whale

DN = starttime:1/24/60/60/fs:starttime+(length(Depth)-1)/24/60/60/fs;


tagon = find(Depth>1,1,'first');
if isempty(tagon); tagon = 1; end
tagoff = find(Depth>1,1,'last');
if isempty(tagoff); tagoff = length(Depth); end
figure(41); clf;
sp2 = subplot(4,1,4);
plot(DN,At);
ylabel('Acceleration'); xlabel('Time');

sp1 = subplot(4,1,1:3);
plot(DN,Depth); hold on; set(gca,'ydir','rev'); try ylim([0 max(Depth)]); catch; ylim([-5 5]); end
ylabel('Depth'); xlabel('Samples');
linkaxes([sp1 sp2],'x');
pat = rectangle('position',[DN(tagon) -10 (tagoff-tagon)/fs/24/60/60 600],'facecolor',[1 1 .4]);
oi = get(gca,'children'); oi=[oi(2:end); oi(1)];
set(gca,'children',oi,'xticklabel',[]);
set(sp2,'xticklabel',datestr(get(gca,'xtick'),'mm/dd HH:MM:SS'));
TEX = text(DN(1),max(Depth),'LEFT click to CHANGE a boundary for tag on/tag off, RIGHT click to ZOOM in (closest boundary will be affected)','verticalalignment','bottom','fontweight','bold','horizontalalignment','left');
[x,~,button] = ginput(1);
zoom reset; fact = (tagoff-tagon)/5/fs/60;
while ~isempty(button)
    % regraph
    delete(pat);
    pat = rectangle('position',[DN(tagon) -10 (tagoff-tagon)/fs/24/60/60 600],'facecolor',[1 1 .4]); %patch([tagon(i) tagoff(i) tagoff(i) tagon(i)],[-10 -10 600 600],[255 255 100]/255);
    oi = get(gca,'children'); oi=[oi(2:end); oi(1)]; set(gca,'children',oi);
    if button == 3
        fact = fact/2;
        if fact>10 && fact < (tagoff-tagon)/5/fs/60/2; fact = 10; end
        xlim([x-fact/24/60 x+fact/24/60]); %surrounds the point by fact minutes
        set(sp2,'xticklabel',datestr(get(gca,'xtick'),'mm/dd HH:MM:SS'));
        delete (TEX); TEX = text(min(get(gca,'xlim')),max(get(gca,'ylim')),'LEFT click to CHANGE a boundary for tag on/tag off, RIGHT click to ZOOM in','verticalalignment','bottom','fontweight','bold','horizontalalignment','left');
        [x,~,button] = ginput(1);
        if button == 3; continue; end
    end
    if button == 1
        [~,e] = min(abs(DN([tagon tagoff])-x));
        if e==2; [~,tagoff] = min(abs(DN-x));
        else [~,tagon] = min(abs(DN-x));
        end
        zoom out; delete(pat); pat = nan(size(tagon));
        pat = rectangle('position',[DN(tagon) -10 (tagoff-tagon)/fs/24/60/60 600],'facecolor',[1 1 .4]); %patch([tagon(i) tagoff(i) tagoff(i) tagon(i)],[-10 -10 600 600],[255 255 100]/255);
        oi = get(gca,'children'); oi=[oi(2:end); oi(1)]; set(gca,'children',oi);
        delete (TEX); TEX = text(min(get(gca,'xlim')),max(get(gca,'ylim')),'LEFT click to CHANGE a boundary for tag on/tag off, RIGHT click to ZOOM in','verticalalignment','bottom','fontweight','bold','horizontalalignment','left');
        set(sp2,'xticklabel',datestr(get(gca,'xtick'),'mm/dd HH:MM:SS'));
        [x,~,button] = ginput(1); fact = (tagoff-tagon)/5/fs/60; continue;
    end
end

oi = zeros(size(Depth));
oi(tagon:tagoff) = 1;
tagon = logical(oi);
disp(['Tagon: ' datestr(DN(find(tagon,1)),'mm/dd/yy HH:MM:SS.fff')]);
disp(['Tagoff: ' datestr(DN(find(tagon,1,'last')),'mm/dd/yy HH:MM:SS.fff')]);