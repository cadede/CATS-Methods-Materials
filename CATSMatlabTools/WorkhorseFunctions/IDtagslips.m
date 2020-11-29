function [tagslip, startsI, endsI] = IDtagslips(DN,At,Depth,fs,tagondec,tagslip,camondec)

if sum(isnan(Depth))>0
    Depth = fixgaps(Depth);
end

tag1 = find(tagondec==1,1,'first');
tag2 = find(tagondec==1,1,'last');
figure(100); clf;
set(100,'windowStyle','docked');
smoothp = runmean(Depth,round(fs/2)); 
surfs = peakfinder(smoothp(tagondec),1,2,-1);
surfs = surfs + tag1 - 1;
Asurf = nan(length(surfs),3);
for i = 1:length(surfs)
    Asurf(i,:) = nanmean(At(max(1,surfs(i)-round(5*fs)):min(surfs(i)+round(5*fs),size(At,1)),:));
end
J = njerk(At,fs);
button = 100;

while ~isempty(button)
    pats = nan(0,3);
    zz = figure(100);  clf;
    s1 = subplot(8,1,1:5);
%     for k = 1:length(pats(:)); try delete(pats(k)); catch; end; end
%     try delete(p1); delete(p3); catch; end
    p1 = plot(Depth); hold on; set(gca,'ydir','rev','xlim',[1 tag2],'ylim',[-7,max(Depth)]);
    if tagslip(1,2)<=tag1; tagslip(tagslip<tag1) = tag1; else tagslip = [tag1 tag1; tagslip]; end
    if tagslip(end,1) >= tag2; tagslip(tagslip>tag2) = tag2; else tagslip(end+1,1:2) = [tag2 tag2]; end
    pats = plotslips(tagslip);
    ys = get(gca,'ylim');
    oi = 1:length(camondec);
    p3 = plot(oi(camondec),camondec(camondec)*-5,'gs','linewidth',3,'markerfacecolor','g','markersize',2);
    if pats(1,2) == 0; pats(1,2) = rectangle('position',[tag1 -50 1 1],'facecolor',[1 1 .4]); set(gca,'ylim',ys); end
    try legend([p1 pats(1,1) pats(1,2) p3], 'Depth','Tag Slips','Tag Slip Duration','Cam On','location','Southeast');
    catch
        legend([p1 pats(1,1) pats(1,2)], 'Depth','Tag Slips','Tag Slip Duration','location','Southeast');
    end
    
    title('Left click on an apparent tag slip location to add a location, click within a tagslip to adjust a spot, press enter when finished');
    s2 = subplot(8,1,6:8);
    
    plot(surfs,Asurf,'s','markersize',2); set(gca,'xlim',[1 tag2]);
    hold on;
    plot(surfs,runmean(Asurf,3),'linewidth',2);
    ylim([-1.1 1.1])
    legend('X @ surface','Y @ surface','Z @ surface','orientation','horizontal','location','southoutside');
    linkaxes([p1 s2],'x');
    %
    [x,~,button] = ginput(1);
    if isempty(button) || button~=1; continue; end
    t1 = find(tagslip(:,2)<x,1,'last');
    t2 = find(tagslip(:,1)>x,1,'first');
    if t2 == t1+2; t = t1+1;
    elseif any(abs(tagslip-x)<30*fs) % if you're within 30 seconds of a tagslip, assume the click was supposed to be on that tag slip
        [~,t] = min(abs(tagslip(:,1)-x));
    elseif t2==t1 + 1; tagslip = [tagslip(1:t1,:); round([x x]); tagslip(t2:end,:)]; t = t1+1;
    else error ('indicated tag slip is outside of boundaries');
    end
    
    b2 = 1000;
    while ~isempty(b2)
        zs = zoomfig(t,tagslip,DN,Depth,J,At,surfs,Asurf,fs,tagondec);
        axes (zs(2))
        [x2,~,b2] = ginput(1);
        if isempty(b2); continue; end
        [~,x2] = min(abs(DN-x2));
        switch b2
            case 1; tagslip(t,:) = x2;
            case 3; tagslip(t,:) = [];
                b2 = [];
            case 49; tagslip(t,1) = x2;
            case 50; tagslip(t,2) = x2;
        end
    end
    todel = [];
    for ii = 1:size(tagslip,1)-1
        if tagslip(ii,2)>tagslip(ii+1,1); tagslip(ii,2) = tagslip(ii+1,2); todel = [todel ii+1]; end
    end
    tagslip(todel,:) = [];
end
startsI = tagslip(1:end-1,2);
endsI = tagslip(2:end,1);
% startsI = round(sort([tag1; x]));
% endsI = round(sort([x; tag2]));
% plot([startsI startsI]',repmat([0 2*max(smoothp)], length(startsI),1)','k','linewidth',3,'parent',s1);
% if nopress; lsurfs = surfs; end


function pats = plotslips(tagslip)
pats = nan(0,3);
ys = get(gca,'ylim');

for i = 1:size(tagslip,1)
    try pats(i,2) = rectangle('position',[tagslip(i,1) ys(1) tagslip(i,2)-tagslip(i,1) ys(2)-ys(1)],'facecolor',[1 1 .4]); catch; end
    pats(i,1) = plot(tagslip(i,1)*ones(1,2),ys,'k','linewidth',2);
    pats(i,3) = plot(tagslip(i,2)*ones(1,2),ys,'k','linewidth',2);
end
oi = get(gca,'children'); 
 I = ismember(oi,pats(:,2));
 oi=[oi(~I); oi(I)];
set(gca,'children',oi);

function zs = zoomfig(t,tagslip,DN,Depth,J,At,surfs,Asurf,fs,tagondec)
figure(100+t-1); clf
set(100+t-1,'windowStyle','docked');
ss1 = subplot(3,1,1);
slipI = tagslip(t,1) - 2*60*fs:tagslip(t,1) + 2*60*fs;
[pp1,h1,h2] = plotyy(DN(slipI),Depth(slipI),DN(slipI),J(slipI)); hold on; set(pp1(1),'ydir','rev');%,'ylim',[min(Depth(slipI)),max(Depth(slipI))]);
set(pp1(2),'ycolor','m','ylim',[0 max(J(tagondec))]);
set(h2,'color','m');
pat1 = plotslips(DN(tagslip(t,:))');
title({'press enter to accept, else left click to locate tagslip at a single point, right click to remove tag slip'; 'press "1" to set a left boundary of a slip period, press "2" to set a right boundary'});
ss2 = subplot(3,1,2);
pp2 = plot(DN(surfs),Asurf,'s','markersize',2);% set(gca,'xlim',[1 tag2]);
hold on;
plot(DN(surfs),runmean(Asurf,3),'linewidth',2);
plot(DN,At,':');
pat2 = plotslips(DN(tagslip(t,:))');
ylim([-1.1 1.1])
legend('X @ surface','Y @ surface','Z @ surface','orientation','horizontal','location','southoutside');
set([ss1 pp1 ss2],'xlim',DN([min(slipI) max(slipI)]));
set([ss1 pp1 ss2],'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
ss3 = subplot(3,1,3);
plot(DN(surfs),Asurf,'s','markersize',2);% set(gca,'xlim',[1 tag2]);
ylim([-1.1 1.1])
hold on;
plot(DN(surfs),runmean(Asurf,3),'linewidth',2);
pat4 = plotslips(DN(tagslip));
pat3 = rectangle('position',[DN(slipI(1)) -1.1 DN(slipI(end))-DN(slipI(1)) 2.2],'facecolor',[.5 0 0]);
oi = get(gca,'children');  I = ismember(oi,pat3);  oi=[oi(~I); oi(I)]; set(gca,'children',oi);
title ('whole deployment: red is zoomed region');
set(ss3,'xticklabel',datestr(get(gca,'xtick'),'HH:MM'));
zs = [ss1 ss2];