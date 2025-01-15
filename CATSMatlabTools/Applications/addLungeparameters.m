function Ltable = addLungeparameters(prhfile,lungefile,lungetotal)
dbstop if error
global L
% input a prh file and a lunge file. Saves a table in the lunge file that
% has two approach parameters, mouth opening, mouth closing and filter time
% options.
% lungetotal provides a limit to stop auditing after a certain number of
% lunges (in case sample size is huge and you are only looking for
% biomechanical averages). Default is inf;

if nargin<3; lungetotal = inf; end

%
%
% clear
% addpath(genpath('C:\Users\Ron\Dropbox\CATS Tools\CATSMatlabTools'))
% cd 'E:\Minke filtration times\'
% load('E:\Minke filtration times\MOfiltnewmat')
% load('E:\Minke filtration times\bblungeT.mat')
% outfold = 'E:\Minke filtration times\';
% prhfold = 'E:\Minke filtration times\prh\';
%
% try load([outfold 'bbfiltTday.mat']);filtT = filtTday;
% catch; filtT = MO(:,[1 3:6]);filtTday = filtT;
%     save([outfold 'bbfiltTday.mat'],'filtTday');
% end

load(prhfile);
load(lungefile);
if ~exist('Ltable','var'); Ltable = table(LungeI,LungeDN); end
if isempty(LungeI); warning('No lunges in lunge file'); return; end
if abs(DN(LungeI(1)) - LungeDN(1)) > 1/24/60/60; error('LungeI does not match prh times. Lunge file may have been created from a prh file with a different start index'); end
try A1 = Ltable.Approach1; catch; A1 = nan(size(LungeI)); end
try A2 = Ltable.Approach2; catch; A2 = nan(size(LungeI)); end
try MO = Ltable.MO; catch; MO = LungeI; end
try MC = Ltable.MC; catch; MC = nan(size(LungeI)); end
try filt = Ltable.filt; catch; filt = nan(size(LungeI)); end
try EX = Ltable.exclude; catch; EX = false(size(LungeI)); end

%        MOI = cellfun(@(x) strcmp(x,minkes{i}),filtT.ID);
%     load([prhfold(1:end-1) 'lunges\' minkes{i} '_lunges.mat']);
%     LI = filtT.MO(MOI); MC = filtT.MC(MOI); MO = LI;
%     lI = cellfun(@(x) strcmp(x,minkes{i}),lungeT.ID);
closeI = @(x) find(min(abs(DN-x)) == abs(DN-x));
%     lungeI = arrayfun(@(x) closeI(x),lungeT.DN(lI));
if length(LungeI)<lungetotal
    lungetotal = length(LungeI);
end
%     lungeInoMO = lungeI; xx = [];
%     for ii = 1:length(LI);
%         if min(abs(LI(ii)-lungeI))<5*fs; % get rid of lunges already used
%             [~,remI] = min(abs(LI(ii) - lungeI));
%             xx = [xx remI];
%         end
%     end
%     lungeInoMO(xx) = [];
% Can use this section to limit your lunge randomization
%     II = lungeI(find(p(lungeI)>20 & lungeT.Daybin(lI)>2&lungeT.Daybin(lI)<13 & ismember(lungeI,lungeInoMO)));
sI = randperm(length(LungeI)); sI = sI(1:lungetotal); sI = sort(sI);
%     LI = LungeI(sI);
%     if length(LI)<lungetotal;
%         LI = [LI; II(sI(1:lungetotal-length(LI)))]; sI(1:lungetotal-length(LI)) = []; %sI gets smaller when used
%     end
kk = 1;
speedJJ = speed.JJ; speedJJ = fixgaps(speedJJ); speedJJ(isnan(speedJJ)) = 0; speedJJ = runmean(speedJJ,fs/2);
DT = datetime(round(DN*24*60*60*10)/24/60/60/10,'convertfrom','datenum','format','dd-MMM-yyyy HH:mm:ss.S');
while kk <=lungetotal %&& length(sI)>=0;
    figure(1); clf;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    sp1 = subplot(3,1,1);
    k = sI(kk);
    I = LungeI(k)-45*fs:LungeI(k)+120*fs;
    I(I<1) = []; I(I>length(DN)) = []; %ensures windows doesn't ask to go beyond the length of the data file
    [ax,h1,h2] = plotyy(DT(I),p(I),DT(I),speedJJ(I));
    set(ax(1),'ydir','rev');
    set(ax,'nextplot','add','xlim',[DT(I(1)) DT(I(end))]); 
    L = plotmarkers(ax,DT,speedJJ,LungeI,MO,MC,A1,A2,k,filt);
    ylabel('Depth (m)'); ylabel('SpeedJJ (m/s)','parent',ax(2)); % rol
    title([INFO.whaleName ', lunge num ' num2str(k)]);
    legend(L,'Lunges','Current Lunge','Engulfment','Approach 1','Approach 2','Filter end')
    instructions = sprintf('Controls:\nEnter: Move to next lunge\np: go to previous lunge\na: Approach 1\n2: Approach 2\nLeftClick: Mouth open\nRightClick: Mouth close\nf: Filter end\nDelete: clear closest marker\nX: Exclude current lunge from analysis');
     annotation('textbox', [0, 0.5, 0, 0], 'string', instructions,'FitBoxToText','on') % adds annotation box
    
    % in the second graph, plot pitch, roll, heading
    sp2 = subplot(3,1,2);
    [ax2,h1,h2] = plotyy(DT(I),pitch(I)*180/pi,DT(I),roll(I)*180/pi); % convert radian metrics to degrees
    set(h1,'color','g','linewidth',2);
    set(ax2(1),'ycolor','g','ylim',[-90 90]);
    set(h2,'color','r','marker','.','markersize',4,'linestyle','none'); % make dots instead of lines for roll and head so that when either crosses -180 to 180, it does not connect the lines
    set(ax2(2),'nextplot','add','ycolor','r','ylim',[-180 180]);
    plot(ax2(2),DT(I),head(I)*180/pi,'b','marker','.','markersize',4,'linestyle','none'); % plot heading as blue
    ylabel('Pitch ({\circ})'); ylabel('Roll and Head ({\circ})','parent',ax2(2)); % roll and heading are on a -180 to -180 scale while pitch is on a -90 to 90 scale.  The /circ indicator is the latex interpreter for degree
    %
    % plot other data useful for identifying events.  Here we plot y-axis
    % gyroscope (useful for identifying up/down tail beats) and jerk data
    highP = Gw(:,2);
    subplot(3,1,3);
    % first determine strokes and glides by choosing a threshold of detection
    lowpassfilt = .5; % filter out high frequency noise above 2 Hz
    [fluke_rate,q] = dsf(highP(tagon),fs,lowpassfilt); %tagon = tag on 0 or 1,  fs is the sampling rate
    nf = round(4*fs/lowpassfilt) ; % how many blank spaces to make?
    fc = lowpassfilt/(fs/2) ; %specifies the cut-off frequency in Hz of a low-pass filter
    %		 to apply to X before computing the spectra. This prevents high frequency
    %		 transients e.g., in foraging, from dominating the spectra.
    filtGw = fir_nodelay(highP,nf,fc) ;
    filtGw(isnan(filtGw)) = 0;
    filtA = fir_nodelay(Aw,nf,fc) ;
    [pitch,roll] = a2pr([filtA(:,1:2),-filtA(:,3)]);
    roll = -roll;
    
    reallylowpassfilt = 0.5*fluke_rate; % separates body posture (< than this frequency) from rotations (> than this frequency)
    fc = reallylowpassfilt; % high speed filter
    nf = round(4*fs/fc) ;
    fc = fc/(fs/2) ;
    bodyGw = fir_nodelay(filtGw,nf,fc);
    
    
    %       flukeamp = nan(size(p));
    % %       flukeamp(tagon) = wrapToPi(head2(tagon)-bodyhead(tagon));
    % %     flukespeed = nan(size(p));
    % %     flukespeed(tagon) = filtGw-[wrapToPi(diff(bodyhead(tagon))); 0]*fs;
    [glds,stks] = stroke_glide([filtGw(tagon)-bodyGw(tagon) roll(tagon) head(tagon)],fs,1,1*pi/180,15);
    stks(:,1) = round(stks(:,1)*fs+find(tagon,1)-1); % convert to indices
    glds = round(glds*fs+find(tagon,1)-1);
    gldI = false(size(p));
    for j = 1:size(glds,1); gldI(glds(j,1):glds(j,2)) = true; end
    upstks = stks(stks(:,2) == 1,:);
    
    
    sp3 = subplot(3,1,3);
    J = njerk(Aw,fs);
    I2 = false(size(DN)); I2(I) = true;
    I2(find(tagon,fs*10)) = false; I2(find(tagon,fs*10,'last')) = false; % exclude the first and last 10 seconds of jerk to exclude tag deployment and detachment
    [ax3,h1,h2] = plotyy(DT(I2),J(I2),DT(I),(filtGw(I)-bodyGw(I))*180/pi);
    set(h1,'color','m'); % type "help plot" to see some of the color options available
    set(ax3(1),'ycolor','m');
    set(h2,'color','g','linewidth',1);
    set(ax3(2),'ycolor','g','nextplot','add');
    plot(ax3(2),DT(upstks(:,1)),filtGw(upstks(:,1)),'rx')
    ylabel('|Jerk| (m/s^{3})');  ylabel('y-axis gyro ({\circ}/s)','parent',ax3(2));
    
    AX = [ax ax2 ax3]; % put all axes in one variable to adjust all at once.  If you did not use plotyy, could have also used sp, sp2, sp3 as axes handles, but plotyy creates two new axes handles that are plotted into
    set(AX,'xlim',[DT(I(1)) DT(I(end))]); % set all xlimits to
    
    [x,~,button] = ginput(1);
    while ~isempty(button)
        switch button 
            case 112 % p
            kk = kk-2; button = []; continue;
            case 97 %a
                A1(k) = closeI(datenum(num2ruler(x,ax(1).XAxis)));
            case 50 %2
                A2(k) = closeI(datenum(num2ruler(x,ax(1).XAxis)));
            case 1 % left click
                MO(k) = closeI(datenum(num2ruler(x,ax(1).XAxis)));
                LungeI(k) = closeI(datenum(num2ruler(x,ax(1).XAxis)));
            case 3 % right click
                MC(k) = closeI(datenum(num2ruler(x,ax(1).XAxis)));
            case 102 % f
                filt(k) = closeI(datenum(num2ruler(x,ax(1).XAxis)));
            case 127 % delete
                xx = find(min(abs([A1(k) A2(k) MC(k) filt(k)]-closeI(datenum(num2ruler(x,ax(1).XAxis))))) == abs([A1(k) A2(k) MC(k) filt(k)]-closeI(datenum(num2ruler(x,ax(1).XAxis)))), 1 );
                switch xx
                    case 1; A1(k) = nan;
                    case 2; A2(k) = nan;
                    case 3; MC(k) = nan;
                    case 4; filt(k) = nan;
                end
            case 120
                EX(k) = true; button = []; continue;
            otherwise
                [x,~,button] = ginput(1); continue;
        end
        L = plotmarkers(ax,DT,speedJJ,LungeI,MO,MC,A1,A2,k,filt);
        [x,~,button] = ginput(1);
    end
    kk = kk+1;
    Ltable.Approach1 = A1; Ltable.Approach2 = A2; Ltable.MO = MO; Ltable.MC = MC; Ltable.filt = filt; Ltable.exclude = EX;
    save(lungefile,'Ltable','-append');
end
    
end

    function L = plotmarkers(ax,DT,speedJJ,LungeI,MO,MC,A1,A2,k,filt)
    global L
    try delete(L); catch; end
    try oi = get(ax(2),'children'); delete(oi(1:6)); catch; end
    L(1) = plot(ax(2),DT(LungeI),speedJJ(LungeI),'rx');
    L(2) = plot(ax(2),DT(LungeI(k)),speedJJ(LungeI(k)),'rs','markerfacecolor','r');
    if ~isnan(MC(k))
        L(3) = plot(ax(2),DT(MO(k):MC(k)),speedJJ(MO(k):MC(k)),'k','linewidth',4);
    else
        L(3) = plot(ax(2),DT(1:2),[0 0],'k','linewidth',4);
    end
    if ~isnan(A1(k))
        L(4) = plot(ax(2),DT(A1(k)),speedJJ(A1(k)),'ks','markerfacecolor','k');
    else
        L(4) = plot(ax(2),DT(1),0,'ks','markerfacecolor','k');
    end
    if ~isnan(A2(k))
        L(5) = plot(ax(2),DT(A2(k)),speedJJ(A2(k)),'ms','markerfacecolor','m');
    else
        L(5) = plot(ax(2),DT(1),0,'ms','markerfacecolor','m');
    end
     if ~isnan(filt(k))
        L(6) = plot(ax(2),DT(filt(k)),speedJJ(filt(k)),'gs','markerfacecolor','g');
    else
        L(6) = plot(ax(2),DT(1),0,'gs','markerfacecolor','g');
     end
    end
  