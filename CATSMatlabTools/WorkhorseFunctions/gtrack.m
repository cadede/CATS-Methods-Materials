function [T,pT,newspeed,newhead] = gtrack(pitch, head, p, fs,speed,tagon,DN,GPS,GPSerr,errThresh,adjspeed,excludeloops)
%%
% outputs: T is the georeferenced pseudotrack with three columns of [m East, m North, m depth];
% outputs: pT is the uncorrected pseudotrack with three columns of [m East, m North, m depth]; THIS IS DIFFERENT THAN JOHNSON ptrack [north, east,depth];
% newspeed and newhead are the corrections made to speed and head to fit the pseudotrack to the given points.
tagonO = tagon; % tagon original (to get rid of GPS points in the situation where the tag came off and then was redeployed.
tagon(find(tagon,1):find(tagon,1,'last')) = true;
% sp = speed.JJ; nhead = head + 30.10*pi/180;
nhead = fixgaps(head);
nhead(isnan(nhead)) = 0;
pthresh = 2;
% gets rid of spikes in speed at the surfacing
sp = speed;
sp2 = sp; sp2(isnan(sp2)) = min(sp2);
sp2 = runmean(sp2,2*fs);
sp2(p<=pthresh & p>0.5) = nan; % but only want to get rid of high speeds in case low speeds have been entered for sleeping whales
sp2 = fixgaps(sp2);
sp(p<=pthresh& p>0.5) = sp2(p<=pthresh& p>0.5); sp = fixgaps(sp); sp(isnan(sp)) = min(sp);
speed = sp;
figure(40); clf; ax = plotyy(DN,p,DN,sp);

Gi = find(~isnan(GPS(:,1)));
G0 = find(Gi<=find(tagon,1),1,'last');
if isempty(G0) && Gi(1)-10*fs<=find(tagon,1); G0 = 1; end
try [~,i] = min([abs(Gi(find(Gi<=find(tagon,1,'last'),1,'last'))-find(tagon,1,'last')),abs(Gi(find(Gi>=find(tagon,1,'last'),1,'first'))-find(tagon,1,'last'))]);
    if i == 1; Ginf = find(Gi<=find(tagon,1,'last'),1,'last'); else Ginf = find(Gi>=find(tagon,1,'last'),1,'first'); end
catch; Ginf = find(Gi<=find(tagon,1,'last'),1,'last');
end

T = ptrack(pitch(tagon), nhead(tagon), p(tagon), fs, [], sp(tagon));
T0 = nan(length(p),3);
T0(tagon,:) = T;
t = T0(:,1);
T0(:,1) = T0(:,2); T0(:,2) = t; clear t;
pT = T0;
figure(102); clf;
set(gcf,'windowstyle','normal');
set(gcf,'units','normalized','outerposition',[0 0 .9 1]);
set(gcf,'units','pixels');
oi = get(gcf,'outerposition');
set(gcf,'outerposition',[oi(1:2) oi(4)*1.1 oi(4)])
p1 = plot(T0(tagon,1),T0(tagon,2)); hold on;
set(gca,'Position',[0.13 0.11 .775 .815]);

if isempty(G0)
    % if there is not a tagon GPS location
    [x1,y1,z1] = deg2utm(GPS(Gi(1),1),GPS(Gi(1),2));
    x1 = x1 - pT(Gi(1),1);
    y1 = y1 - pT(Gi(1),2);
    [GPS(find(tagon,1)-1,1),GPS(find(tagon,1)-1,2)] = utm2deg(x1,y1,z1);
    Gi = [find(tagon,1)-1; Gi];
    G0 = 1;
    disp(['No tagon location, location calculated from ptrack: ' num2str(GPS(Gi(1),:))]);
else
    GPS(find(tagon,1)-1,:) = GPS(Gi(G0),:);
    GPSerr(find(tagon,1)-1,:) = GPSerr(Gi(G0),:);
    Gi(G0) = find(tagon,1)-1;
end

if nargin<10 || isempty(errThresh)
    errThresh = max([std(abs(GPSerr(Gi(G0:Ginf),:)),[],1)/8; [10 10 10]]);
end
if nargin<12  % an attempt to exlude bubble nets from adjustment (leave circles as the pseudo tracks).  Mixed success
    excludeloops = false;
end
if nargin<11 || isempty(adjspeed)
    adjspeed = true;
end
todel = [];
for i = G0+1:Ginf
    if any(abs(GPSerr(Gi(i),:))>errThresh) || ~ismember(Gi(i),find(tagonO)) % new: get rid of Gi that is not when the tag is on (good for if there is a gap where the tag was off and then redeployed).
        todel = [todel; i];
    end
end
Gi(todel) = [];
Ginf = Ginf - length(todel);





[x,y,uzone] = deg2utm(GPS(Gi,1),GPS(Gi,2));
[xp,mainUZ] = utm2m(x,y,uzone);
% if any(~cellfun(@(x) strcmp(x,uzone(1,:)),mat2cell(uzone,ones(size(uzone,1),1),4)))
%     disp('GPS hits span UTM zones, this may cause problems on axes but all else should be fine');
% end

x0 = xp(G0); y0 = y(G0);
xp = xp-xp(G0);
y = y-y(G0);
G = plot(xp(G0:Ginf),y(G0:Ginf),'ks');
set(G,'markerfacecolor','k');
axis equal

if excludeloops
loops = isnet(head,pitch,fs);
else loops = false(size(p));
end

%         sp(p<1) = 0.05;
%         speed(p<1) = 0.05;

T = T0;
T(G0,[1 2]) = [0 0];
for i = G0:Ginf-1
    I1 = Gi(i); I2 = Gi(i+1);
    if I2>find(~isnan(T(:,1)),1,'last'); I2 = find(~isnan(T(:,1)),1,'last'); end
    e = [xp(i+1) y(i+1)]-T(I2,1:2);
    %find min speed threshold
    if adjspeed
%         sp(p<1) = 0.05;
%         speed(p<1) = 0.05;
        spT = prctile(sp(I1:I2),10)+.1;
        
        %     spI = sp<=spT|p<=pthresh; % speeds to try adjusting
        for ii = .4:-0.1:0
            sp2 = sp;
            sp2(sp<=spT) = sp2(sp<=spT)*ii;
            sp2(p<=pthresh) = spT*ii;
            sp2([1:I1-1 I2+1:end]) = sp([1:I1-1 I2+1:end]);
            Tadj = T;
            Tadj(I1:I2,:) = ptrack(pitch(I1:I2) , nhead(I1:I2), p(I1:I2) , fs, [], sp2(I1:I2));
            t = Tadj(I1:I2,1); Tadj(I1:I2,1) = Tadj(I1:I2,2); Tadj(I1:I2,2) = t;
            if I1 == 1 || any(isnan(T(I1-1,[1 2]))); oi = [0 0]; else oi = T(I1-1,1:2); end
            oi = repmat(oi,I2-I1+1,1);
            Tadj(I1:I2,1:2) = Tadj(I1:I2,1:2) + oi;
            oi = repmat(T(I2,1:2)-Tadj(I2,1:2),size(Tadj,1)-I2,1);
            Tadj(I2+1:end,1:2) = Tadj(I2+1:end,1:2) -oi;
            e2 = [x(i+1) y(i+1)]-Tadj(I2,1:2);
            if sqrt(e2(1)^2+e2(2)^2)<sqrt(e(1)^2+e(2)^2)
                e = e2; T = Tadj; %speed = sp2;
            end
        end
    end
    %sp = speed;    
    
    E = zeros(size(T(:,[1 2])));
    E(I2+1:end,1) = e(1);
    E(I2+1:end,2) = e(2);
    inc = repmat(e/(diff([I1 I2])-sum(loops(I1:I2))),I2-I1,1); %don't add errors to the loops
    inc(loops(I1:I2)) = 0;
%     inc = repmat(e/diff([I1 I2]),I2-I1,1); %use this line instead of two above to not do the loop thing
    inc = cumsum(inc);
    E(I1+1:I2,:) = inc;
    T(:,[1 2]) = T(:,[1 2]) + E;
end
plot(T(tagon,1),T(tagon,2),'m');
I1 = find(tagon,1); I2 = find(tagon,1,'last');
t1 = plot(T(I1:fs*60*30:I2,1),T(I1:fs*60*30:I2,2),'mo','markerfacecolor','m','markersize',4);
t2 = plot(T0(I1:fs*60*30:I2,1),T0(I1:fs*60*30:I2,2),'bo','markerfacecolor','b','markersize',4);
plot(T(I1:fs*60*5:I2,1),T(I1:fs*60*5:I2,2),'mo','markerfacecolor','m','markersize',2);
p3 = plot(T0(I1:fs*60*5:I2,1),T0(I1:fs*60*5:I2,2),'bo','markerfacecolor','b','markersize',2);
p4 = nan(0,0);
for jj = I1:fs*60*30:I2
    text(T(jj,1),T(jj,2),datestr(DN(jj),'HH:MM'),'verticalalignment','bottom','horizontalalignment','center','fontsize',6);
    p4(end+1) = text(T0(jj,1),T0(jj,2),datestr(DN(jj),'HH:MM'),'verticalalignment','bottom','horizontalalignment','center','fontsize',6);
end
    text(T(I1,1),T(I1,2),'START','verticalalignment','top','horizontalalignment','center','fontsize',22);

set(ax(2),'nextplot','add');
plot(ax(2),DN,speed,'r--');
set(ax(2),'xticklabel',[]); set(ax(1),'xticklabel',datestr(get(ax(1),'xtick'),'HH:MM:SS'));
linkaxes(ax,'x');

figure(102);
xs = get(gca,'xtick');
if length(xs)<7; xs = xs(1):diff(xs(1:2))/2:xs(end); end; set(gca,'xtick',xs);
ys = get(gca,'ytick');
if length(ys)<7; ys = ys(1):diff(ys(1:2))/2:ys(end); end; set(gca,'ytick',ys);

xs = get(gca,'xtick')';
xs = xs(1:2:end);
set(gca,'xtick',xs);
ys = get(gca,'ytick')';
% if length(ys)<size(uzone,1); UZ = uzone(1:size(ys,1),:);
% else
    UZ = repmat(mainUZ,length(ys),1); %uzone(1,:)
% end
dys = utm2deg(zeros(size(ys,1),1)+x0,ys+y0,UZ);
% if length(xs)<size(uzone,1); UZ = uzone(1:size(xs,1),:);
% else
    UZ = repmat(mainUZ,length(xs),1); %uzone(1,:)
% end
[~,dxs] = utm2deg(xs+x0,zeros(size(xs,1),1)+y0,UZ);
oi = sprintf('%0.3f',dxs); oi = reshape(oi,size(oi,2)/size(dxs,1),size(dxs,1))';
xlabs = [oi repmat([char(176) ' ('],size(xs,1),1) num2str(xs) repmat(' m)',size(xs,1),1)];
oi = sprintf('%0.3f',dys); oi = reshape(oi,size(oi,2)/size(dys,1),size(dys,1))';
ylabs = [oi repmat([char(176) ' ('],size(ys,1),1) num2str(ys) repmat(' m)',size(ys,1),1)];
set(gca,'xticklabel',xlabs);
set(gca,'yticklabel',ylabs);
axis auto 
legend([t2 t1],'ptrack (5 min)','geo-ptrack (5 min)');
newf = copyfig(102);
set(newf,'windowstyle','normal');
figure(102);
delete ([p1 t2 p3 p4]);
set(gca,'xtickmode','auto','ytickmode','auto');
axis equal
xs = get(gca,'xtick');
if length(xs)<7; xs = xs(1):diff(xs(1:2))/2:xs(end); end; set(gca,'xtick',xs);
ys = get(gca,'ytick');
if length(ys)<7; ys = ys(1):diff(ys(1:2))/2:ys(end); end; set(gca,'ytick',ys);

xs = get(gca,'xtick')';
xs = xs(1:2:end);
set(gca,'xtick',xs);
ys = get(gca,'ytick')';
% if length(ys)<size(uzone,1); UZ = uzone(1:size(ys,1),:);
% else UZ = repmat(uzone(1,:),length(ys),1);
% end
UZ = repmat(mainUZ,length(ys),1); %uzone(1,:)
dys = utm2deg(zeros(size(ys,1),1)+x0,ys+y0,UZ);
% if length(xs)<size(uzone,1); UZ = uzone(1:size(xs,1),:);
% else UZ = repmat(uzone(1,:),length(xs),1);
% end
UZ = repmat(mainUZ,length(xs),1); %uzone(1,:)
[~,dxs] = utm2deg(xs+x0,zeros(size(xs,1),1)+y0,UZ);
oi = sprintf('%0.3f',dxs); oi = reshape(oi,size(oi,2)/size(dxs,1),size(dxs,1))';
xlabs = [oi repmat([char(176) ' ('],size(xs,1),1) num2str(xs) repmat(' m)',size(xs,1),1)];
oi = sprintf('%0.3f',dys); oi = reshape(oi,size(oi,2)/size(dys,1),size(dys,1))';
ylabs = [oi repmat([char(176) ' ('],size(ys,1),1) num2str(ys) repmat(' m)',size(ys,1),1)];
set(gca,'xticklabel',xlabs);
set(gca,'yticklabel',ylabs);
axis auto 


% back out implied speed and heading and pitch
newhead = head;
newspeed = speed;
newhead(2:end) = atan2(diff(T(:,1)),diff(T(:,2)));
newspeed(2:end) = diff(T(:,1))./sin(newhead(2:end))./cos(min(80*pi/180,pitch(2:end)))*fs; %puts a threshold in there to divide by to get speed since at vertical speeds you can divide by nearly 0
% for i = I1:I2;
%     newhead(i) = atan2((T(i+1,1)-T(i,1))/(T(i+1,2)-T(i,2)));
%     newspeed(i) = (T(i+1,1)-T(i,1))/sin(nhead(i))*fs;
% end
figure(41); clf;
s1 = subplot(2,1,1); plot(DN,speed,DN,newspeed,'g--');set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));legend('Original Speed','Adjusted Speed'); ylim([0 10]);
s2 = subplot(2,1,2); plot(DN,head*180/pi,DN,newhead*180/pi,'g--'); set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));legend('Original heading','Adjusted heading');
linkaxes([s1 s2],'x');
figure(newf);

 [x1,y1,z1] = deg2utm(GPS(Gi(G0),1),GPS(Gi(G0),2));
    x1 = x1 + T(find(tagon,1,'last'),1);
    y1 = y1 + T(find(tagon,1,'last'),2);
    [l1,l2] = utm2deg(x1,y1,z1);
%     Gi = [find(tagon,1)-1; Gi];
%     G0 = 1;
    disp(['End of track: ' num2str([l1 l2])]);

