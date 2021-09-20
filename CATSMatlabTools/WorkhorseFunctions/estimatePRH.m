
function [Aw,Mw,Gw,W,Wchange,Wchangeend,tagprh,pitch,roll,head,calperiodI,tagslip,speedper] = estimatePRH(At,Mt,Gt,fs,DN,Depth,tagon,dec,tagslip,calperiodI,W)
% calperiodI and W are optional inputs as they are calculated here, but can be inputted from previous runs
% dec is the declination in radians at the deployment location.  Enter 0 if you want heading calculated to magnetic north
% tagslip is a m x 2 matrix, where the first row and last rows are the tag on and tag off indices (same value in each column).  Column 1 is the start of a tag slip, column 2 is the end of a tagslip.

% this function walks the user through identifying periods of time useful
% to convert tag frame sensor data to animal frame.  See description around
% Figs 7 & 8 in the manuscript accompanying this code.

% Johnson method based on a tutorial from: Fine-scale animal
% movement workshop, Hobart, March 2011.  Page 10 gives a lot of good information and
% summarizes the processing steps required to estimate an animal's
% orientation.
% The biggest change is that for consistency in rotation, we use a NED
% reference frame, meaning gravity is -1 normally.  When you face an axis
% up, the reading will be positive, when you face down, negative.  The
% rotations are counter clockwise when looking from the positive side of
% the third axis.  The resulting rotation matrices are these:
% H = [cos(h) -sin(h) 0; sin(h) cos(h) 0; 0 0 1];
% P = [cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)];
% R = [1 0 0; 0 cos(r) -sin(r); 0 sin(r) cos(r)];
% so P is P' from Johnson, but consistent with other right-hand rotations
% from other sources.


global fileloc filename

if ~exist('calperiodI','var') || isempty(calperiodI); calperiodI = cell(size(tagslip,1)-1,1); end
if ~exist('W','var') || isempty(W); W = cell(size(tagslip,1)-1,1);
    [W,calperiodI,tagslip,tocal,startsI,endsI,Wchange,Wchangeend] = reconcileSlips(W,calperiodI,tagslip);
    Aw = nan(size(At)); Mw = Aw; Gw = Aw;
else
    if length(W)~=length(calperiodI); error('check length of W'); end
    [W,calperiodI,tagslip,tocal,startsI,endsI,Wchange,Wchangeend] = reconcileSlips(W,calperiodI,tagslip);
    [Aw,Mw,Gw] = applyW(W,startsI,endsI,At,Mt,Gt,true);
end

if any(endsI<startsI); error('check startsI and endsI of calibration periods'); end
if size(tagslip,1)~=length(startsI)+1; error('tag slips should include tag on and tag off and must be start and end of each tag slip period'); end
if size(tagslip,2)~= 2; error('tag slip must be two columns with start and end of tag slips in each column'); end

tag1 = find(tagon,1); tag2 = find(tagon,1,'last');
instructions = sprintf('Controls:\nLeftClick: Calibrate clicked section\nRight Click: rotate dolphin model\n#1-9: zoom to x minutes on either side\nt: enter adjust tag slip mode \nEnter: Finish (accept calibrations)');
inst1 = sprintf('Controls:\nLeftClick: zoom to 10 minutes on either side\n#1-9: zoom to x minutes on either side\nRight Click: do not zoom\nEnter: Accept calibration');
inst2 = sprintf('Controls:\nLeftClick L and R boundaries of surfacing (neutral position)\nRightClick L and R boundaries of ascent to surface or \n      start of dive (where animal rotates around y-axis)\nt: enter adjust tag slip mode\nEnter: zoom to main screen');
inst3 = sprintf('Controls:\nt: enter adjust tag slip mode\nRightClick: rotate dolphin model\nPress any key to return to main screen');
button = 1;
figure(100); clf;
tagprh = nan(length(startsI),3);
surfs = []; dives = [];
while ~isempty(button) || any(cellfun(@isempty,W))
    [x,button,ax,pp,rr] = mainGraph(At,Aw,Mw,DN,Depth,tagslip,tocal,tag1,tag2,instructions);
    % pp & rr are current pitch and roll
    if isempty(button) && any(cellfun(@isempty,W)); title('Not all periods are calibrated, are you sure you wish to quit? Press enter to confirm or any key to continue.','fontweight','bold','fontsize',14,'parent',ax)
        [~,~,button] = ginput(1); if isempty(button); break; else continue; end;
    elseif isempty(button); continue;
    end
    [~,xx] = min(abs(DN-x));
    switch button
        case 1 % calibrate section
            i = find(startsI<xx,1,'last');
            accepted = false;
            if ~isempty(calperiodI{i}); surfs = calperiodI{i}(1:end-2); surfs = reshape(surfs,2,length(surfs)/2)'; dives = calperiodI{i}(end-1:end);
            else surfs = nan(0,2); dives = surfs;
            end
            while ~accepted
                [x,button,ax] = mainGraph(At,Aw,Mw,DN,Depth,tagslip,tocal,startsI(i),endsI(i),inst1,surfs,dives);
                if isempty(button); if ~isempty(surfs)&&~isempty(dives); accepted = true; continue;
                    else title('Not all cal periods have been determined, press any key to select again','fontweight','bold','fontsize',14,'parent',ax); ginput(1); continue;
                    end;
                end
                [~,xx] = min(abs(DN-x));
                switch button
                    case 1; a = xx-10*60*fs; b = xx+10*60*fs;
                    case num2cell(49:57);  a = xx-(button-48)*60*fs; b = xx+(button-48)*60*fs;
                    otherwise; a = startsI(i); b = endsI(i);
                end
                while ~isempty(button)
                     titl = 'Can have multiple surfacings (orange) but only one dive (blue) selected.  Press "del" to remove a section.';
                    [x,button,ax] = mainGraph(At,Aw,Mw,DN,Depth,tagslip,tocal,a,b,inst2,surfs,dives,titl);
                    if isempty(button); continue; end
                    [~,xx] = min(abs(DN-x));
                    if button == 127 || button == 8; % if delete is pressed ( or backspace)
                        if isempty(surfs) && isempty(dives); continue; end
                        [~,ii] = min(abs(xx - [surfs(:); dives(:)]));
                        if ii>length(surfs(:)); dives = []; continue; end
                        surfs(ceil(ii/2),:) = [];
                    elseif button == 116 % add tag slip info
                        [W,calperiodI,tagslip,tocal,startsI,endsI,Wchange,Wchangeend]= adjslips(tagslip,At,Aw,Mw,dec,DN,Depth,tocal,surfs,dives,W,calperiodI);
                        [Aw,Mw,Gw] = applyW(W,startsI,endsI,At,Mt,Gt,true);
                    else
                        [x2,~] = ginput(1); % get a second click
                        [~,xx(1,2)] = min(abs(DN-x2));
                    end
                    if button == 1
                        if isempty(surfs); ii = 1; else ii = find(xx(1)>surfs(:,1),1,'last'); if isempty(ii); ii = 0; end; end
                        if isempty(surfs); surfs = sort(xx);%~isempty(surfs) && ~isempty(ii) && xx(2)<surfs(ii,2); surfs(ii,:) = []; continue;
                        else surfs = [surfs(1:ii,:); sort(xx); surfs(ii+1:end,:)];
                        end
                    elseif button == 3
                        dives = sort(xx);
                    end
                    if ~isempty(surfs)&&~isempty(dives) ; % you have a full set of calibration periods
                        if any([surfs(:); dives']>endsI(i) | [surfs(:); dives'] <startsI(i)); % if the surfs and dives are outside of the calibration period
                            calperiodI{i} = [];
                        else
                            I = startsI(i):endsI(i);
                            s1 = surfs';
                            calperiodI{i} = [s1(:)' dives];
                            [Aw(I,:),Mw(I,:),Gw(I,:),~,~,~,W{i},tagprh(i,:)] = tagframe2whaleframe(At(I,:),Mt(I,:),Gt(I,:),Depth(I),[],calperiodI{i}-startsI(i)+1);
                            [~,~,~,tocal] = reconcileSlips(W,calperiodI,tagslip);
                        end
                        %                          tocal(i,:) = [nan nan];
                    end
                    if isempty(surfs) || isempty(dives) && (~isempty(surfs)&&~isempty(dives)); % if one is empty but not both, reset W. If both are empty, could be just an examination, so don't necessarily 0 out W, but if only one is empty, assume W needs ot be adjusted
                        W{i} = []; tagprh(i,:) = nan(1,3); calperiodI{i} = [];
                         I = startsI(i):endsI(i);
                        Aw(I,:) = nan(length(I),3); Mw(I,:) = nan(length(I),3); Gw(I,:) = nan(length(I),3); 
                        [~,~,~,tocal] = reconcileSlips(W,calperiodI,tagslip);
                    end
                end
            end
            if i~=1 && ~isempty(W{i-1}) && endsI(i-1)~=startsI(i)
                I = endsI(i-1):startsI(i);
                [Aw(I,:),Mw(I,:),Gw(I,:)] = applyWduringslip(At(I,:),Mt(I,:),Gt(I,:),W{i-1},W{i});
            end
            if i~=length(W) && ~isempty(W{i+1}) && endsI(i)~=startsI(i+1)
                I = endsI(i):startsI(i+1);
                [Aw(I,:),Mw(I,:),Gw(I,:)] = applyWduringslip(At(I,:),Mt(I,:),Gt(I,:),W{i},W{i+1});
            end
            save([fileloc filename(1:end-4) 'Info.mat'],'W','calperiodI','-append');
            button =1;
            surfs = []; dives = [];
        case 3; makeDolphin(pp(xx),rr(xx),0); title('heading not yet calculated'); continue
        case num2cell(49:57)
            button2 =3;
            while button2 == 3
                [x,button2,~,pp,rr] = mainGraph(At,Aw,Mw,DN,Depth,tagslip,tocal,xx-(button-48)*60*fs,xx+(button-48)*60*fs,inst3); %[~,xx] = min(abs(DN-x));
                if button2 == 3; makeDolphin(pp(xx),rr(xx),0); title('heading not yet calculated'); 
                elseif button2 == 116 % add tag slip info
                        [W,calperiodI,tagslip,tocal,startsI,endsI,Wchange,Wchangeend]= adjslips(tagslip,At,Aw,Mw,dec,DN,Depth,tocal,surfs,dives,W,calperiodI);
                        [Aw,Mw,Gw] = applyW(W,startsI,endsI,At,Mt,Gt,true);
                end
            end
               continue;
        case 116 % add tag slip info
            [W,calperiodI,tagslip,tocal,startsI,endsI,Wchange,Wchangeend]= adjslips(tagslip,At,Aw,Mw,dec,DN,Depth,tocal,surfs,dives,W,calperiodI);
            [Aw,Mw,Gw] = applyW(W,startsI,endsI,At,Mt,Gt,true);
        otherwise; continue
    end
    
    
end

% Wchange = endsI; Wchangeend = startsI;
speedper = [startsI [startsI(2:end); endsI(end)]]; 
todel = [];
for i = 2:size(speedper,1)
    if ~any(abs(circ_dist(tagprh(i,:),tagprh(i-1,:))*180/pi)>15) %15 is arbitrary.  If the slippage in p,r and h of the tag is < 15 degrees, just consider the tag to have been stable
        todel = [todel; i-1];
        speedper(i,1) = speedper(i-1,1);
    end
end
speedper(todel,:)= [];
[pitch,roll,head] = calcprh(Aw,Mw,dec);
 
function [W,calperiodI,tagslip,tocal,startsI,endsI,Wchange,Wchangeend]= adjslips(tagslip,At,Aw,Mw,dec,DN,Depth,tocal,surfs,dives,W,calperiodI)
titl = 'Tag Slip adjustment';
global fileloc filename
inst = sprintf('Controls:\nLeft Click: create a new tag slip (and cal period)\n1: adjust nearest tag slip start to here\n2: adjust nearest tag slip end to here\n"del" (or backspace): delete nearest tag slip\nEnter: go back to prior mode');
xs = get(gca,'xlim');
for ii = 1:2; [~,xs(ii)] = min(abs(DN-xs(ii))); end
button = 1;
while ~isempty(button)
    [x,button] = mainGraph(At,Aw,Mw,DN,Depth,tagslip,tocal,xs(1),xs(2),inst,surfs,dives,titl,dec);
    titl = 'Tag Slip adjustment';
    if isempty(button); continue; end
    [~,xx] = min(abs(DN-x));
    switch button
        case 1
            ii = find(xx>tagslip(:,2),1,'last');
            if xx>tagslip(ii+1,1); titl = {'Can''t add a tagslip in the middle of another tag slip, try again.'}; continue; end
            tagslip = [tagslip(1:ii,:); [xx xx]; tagslip(ii+1:end,:)];
            [W,calperiodI,tagslip,tocal,startsI,endsI,Wchange,Wchangeend] = reconcileSlips(W,calperiodI,tagslip);
            Aw = applyW(W,startsI,endsI,At,At,At,true);
        case 49 % adj tag slip start
            ii = find(xx>tagslip(:,1),1,'last');
            if xx>tagslip(ii,2); ii = ii+1; end % if you click to the right of a tagslip end, adj the beginning of the next tagslip
            if ii == size(tagslip,1); titl = {'Can''t adjust tagoff slip!'}; continue; end
            tagslip(ii,1) = xx;
            [W,calperiodI,tagslip,tocal,startsI,endsI,Wchange,Wchangeend] = reconcileSlips(W,calperiodI,tagslip);
            Aw = applyW(W,startsI,endsI,At,At,At,true);
        case 50
            ii = find(xx<tagslip(:,2),1,'first');
            if xx<tagslip(ii,1); ii = ii-1; end % if you click to the right of a tagslip end, adj the beginning of the next tagslip
            if ii == 1; titl = {'Can''t adjust tagon slip!'}; continue; end
            tagslip(ii,2) = xx;
            [W,calperiodI,tagslip,tocal,startsI,endsI,Wchange,Wchangeend] = reconcileSlips(W,calperiodI,tagslip);
            Aw = applyW(W,startsI,endsI,At,At,At,true);
        case {8, 127} % if delete is pressed ( or backspace)
            if size(tagslip,1)<3; titl = ('No tagslips to delete!  Select again'); continue; end
            [~,ii] = min(abs(median(tagslip,2)-xx));
            if ii == 1; ii = 2; elseif ii == size(tagslip,1); ii = ii-1; end
            tagslip(ii,:) = []; 
            [W,calperiodI,tagslip,tocal,startsI,endsI,Wchange,Wchangeend] = reconcileSlips(W,calperiodI,tagslip);
            Aw = applyW(W,startsI,endsI,At,At,At,true);
    end
end
tempslips = tagslip;
if ~exist('startsI','var'); startsI = tagslip(1:end-1,2); %if you exit without going through the process of changing anything
    endsI = tagslip(2:end,1);
    Wchange = endsI(1:end-1); Wchangeend = startsI(2:end);
end
save([fileloc filename(1:end-4) 'Info.mat'],'tempslips','-append');
disp('New tag slip locations saved as "tempslips" variable in "INFO" file, load manually if process is interrupted and you wish to use these values.');
warning('"slips" variable will be overwritten with new slips upon successful function termination');
            


function [x,button,s1,pitch,roll] = mainGraph(At,Aw,Mw,DN,Depth,tagslip,tocal,a,b,instructions,surfs,dives,titl,dec) %a,b are the zoom factors
if a<1; a = 1; end; if b>length(Depth); b = length(Depth); end
figure(100); clf
s1 = subplot(3,6,1:5);
p1 = plot(DN,Depth); hold on; set(gca,'ydir','rev','ylim',[-7,max(Depth(a:b))]);
ylabel('Depth');
if exist('surfs','var')&&~isempty(surfs); plotpatches([DN(surfs(:,1)) DN(surfs(:,2))],[255 165 0]/255); end
if exist('dives','var')&&~isempty(dives); plotpatches(DN(dives)',[0 255 255]/255); end
plotpatches([DN(tagslip(:,1)) DN(tagslip(:,2))],[1 1 .4]);
tx = ~isnan(tocal(:,1));
plotpatches([DN(tocal(tx,1)) DN(tocal(tx,2))],[205 92 92]/255);
for i = 1:size(tagslip,1)-1; t = text(DN(tagslip(i,2)),-7,num2str(i), 'verticalalignment','bottom','horizontalalignment','left'); end
if exist('titl','var'); title(titl,'fontsize',14,'fontweight','bold'); end
s2(1) = subplot(6,6,13:17);
plot(DN,At(:,1),'b:',DN,Aw(:,1),'b'); legend('At_{X}','Aw_{X}','orientation','horizontal');
s2(2) = subplot(6,6,19:23);
p2 = plot(DN,At(:,2),':',DN,Aw(:,2)); set(p2,'color',[0 .5 0]); legend('At_{Y}','Aw_{Y}','orientation','horizontal');
s2(3) = subplot(6,6,25:29);
plot(DN,At(:,3),'r:',DN,Aw(:,3),'r'); legend('At_{Z}','Aw_{Z}','orientation','horizontal');
At_mag = sqrt(sum(At.^2,2));
% pitch=asin(Aw(:,1)./At_mag);  % in radians
% roll=atan2(-Aw(:,2),-Aw(:,3));
[pitch,roll] = calcprh(Aw,Mw);
set(s2,'ylim',[-1.5 1.5])
s3 = subplot(6,6,31:35);
[ax,h1,h2] = plotyy(DN,pitch*180/pi,DN,roll*180/pi);
set(ax(1),'ycolor','g','ylim',[-90 90]); set(h1,'color','g');
set(ax(2),'ycolor','r','ylim',[-180 180]); set(h2,'color','r');
if exist('dec','var') && ~isempty(dec)
    [~,~,head] = calcprh(Aw(a:b,:),Mw(a:b,:),dec);
    set(ax(2),'ycolor','k','nextplot','add');
    head(abs(diff(head))*180/pi>340) = nan;
    h3 = plot(ax(2),DN(a:b),head*180/pi,'b');
    legend([h2,h3],'roll','head');
end
ylabel('pitch','parent',ax(1)); ylabel('roll','parent',ax(2));
set([s1 s2 ax],'xlim',DN([a b]));
set([s1, ax(1)],'xticklabel',datestr(get(ax(1),'xtick'),'HH:MM:SS')); set([s2 ax(2)],'xticklabel',[]);
% subplot(3,4,4);
annotation('textbox', [.8, .9, 0, 0], 'string', instructions,'FitBoxToText','on')
axes(s1);
n = 1;
% if strcmp(inst2,instructions); n = 2; end
[x,~,button] = ginput(n);


function pats = plotpatches(tagslip,colors)
pats = nan(0,3);
ys = get(gca,'ylim');

for i = 1:size(tagslip,1)
    try pats(i,2) = rectangle('position',[tagslip(i,1) ys(1) tagslip(i,2)-tagslip(i,1) ys(2)-ys(1)],'facecolor',colors); catch; end
    pats(i,1) = plot(tagslip(i,1)*ones(1,2),ys,'k','linewidth',2);
    pats(i,3) = plot(tagslip(i,2)*ones(1,2),ys,'k','linewidth',2);
end
oi = get(gca,'children'); 
 I = ismember(oi,pats(:,2));
 oi=[oi(~I); oi(I)];
set(gca,'children',oi);


function makeDolphin(p,r,h)%,fig)
figure(6); clf;
F = plot_3d_model;
rot_3d_model_NED(F,[p r h]);%prh(find(tagon,1),:));
% oi = get(gca,'children');
% cm = get(gcf,'colormap');
% v = get(gca,'view');
% % figure(fig);
% subplot(3,5,15);
% copyobj(oi,gca)
% set(gcf,'colormap',cm)
% set(gca,'view',v,'xtick',[],'ytick',[],'ztick',[])
