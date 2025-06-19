function [ax1,ax2,ax3,H1,H2,H3] = plotprh3 (prh,lunges,FIG,gyroAx,nogyro)
global L
% plots a set of 3 subplots using the data in the prh structure
% LungeI is an index of lunges to plot (optional)
% FIG is an optional argument that specifies the figure number to plot in
% gyroAx should be 2 for finding pitch strokes or 3 for finding heading
% strokes
% nogryo = true if you want to plot strokes using the accelerometer
% OUTPUTS:
% ax1-3, the axis handles of each subplot
% H1-3, the plot handles in each subplot

if nargin<5||isempty(nogyro)
    nogyro = false;
end

if nargin<4 || isempty(gyroAx)
    gyroAx = 2;
end
if nogyro; gyroAx = gyroAx-1; end

if nargin<3 || isempty(FIG)
    FIG = figure; 
else; figure(FIG); clf;
end

if nargin<2 || isempty(lunges)
    LungeI = [];
else
    LungeI = lunges.LungeI;
    try Ltable = lunges.Ltable; catch; end
end
fs = prh.fs;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    subplot(3,1,1); % divide the figure into three rows and one column, and plot in the first plot
%     title(ID);
    I = prh.tagon; % limit the plots to only be when the tag is on the animal
    try % the try/catch switch can be used in case your deployment does not have speed, in that case, the "catch" creates an empty variable
        speedJJ = prh.speed.JJ; % create a dummy variable with your speed data and smooth it with a 1 second running mean
        speedJJ(isnan(speedJJ)) = min(speedJJ); speedJJ = runmean(speedJJ,round(prh.fs)); speedJJ(isnan(prh.speed.JJ)) = nan;
    catch; speedJJ = nan(size(prh.p));
    end
    % in the first graph, plot depth (p) and speed
    DN = datetime(round(prh.DN*24*60*60*10)/24/60/60/10,'convertfrom','datenum','format','dd-MMM-yyyy HH:mm:ss.S');
    [ax1,h1,h2] = plotyy(DN(I),prh.p(I),DN(I),speedJJ(I)); % plotyy allows for plotting time series with different y-axes on the same x-axis
    set(h1,'color','b','linewidth',2);
    a = find(I,1); b = find(I,1,'last');
    set(ax1(1),'ylim',[-max(prh.p(I))/10 max(prh.p(I))],'ydir','rev','nextplot','add'); % set the y-limits of the depth plot and plot so 0 is at the top, and set it so that if you plot the next line, it plots on top of the depth plot instead of replacing it
    set(ax1,'xlim',[DN(a) DN(b)])
    try h3 = plot(ax1(1),DN(LungeI),prh.p(LungeI),'rs','markerfacecolor','r'); % if you loaded a file with events indexed as "LungeI", plot those events as red squares on the depth graph
    catch; h3 = nan;
    end
    ylabel('Depth (m)'); ylabel('Speed (m/s)','parent',ax1(2)); % for plotting a ylabel in the non-active axis, need to identify which axis you are plotting into
    set(ax1(2),'nextplot','add');
    try
        for k = 1:size(Ltable,1)
            try L = plotmarkers(ax1,DN,speedJJ,Ltable.LungeI,Ltable.MO,Ltable.MC,Ltable.Approach1,Ltable.Approach2,k,Ltable.filt); catch; end
        end
    catch
    end
        H1 = [h1 h2 h3];




    % in the second graph, plot pitch, roll, heading
    subplot(3,1,2);
    [ax2,h1,h2] = plotyy(DN(I),prh.pitch(I)*180/pi,DN(I),prh.roll(I)*180/pi); % convert radian metrics to degrees
    set(h1,'color','g','linewidth',2);
    set(ax2(1),'ycolor','g');
    set(h2,'color','r','marker','.','markersize',4,'linestyle','none'); % make dots instead of lines for roll and head so that when either crosses -180 to 180, it does not connect the lines
    set(ax2(2),'nextplot','add','ycolor','r');
    h3 = plot(ax2(2),DN(I),prh.head(I)*180/pi,'b','marker','.','markersize',4,'linestyle','none'); % plot heading as blue
    ylabel('Pitch ({\circ})'); ylabel('Roll and Head ({\circ})','parent',ax2(2)); % roll and heading are on a -180 to -180 scale while pitch is on a -90 to 90 scale.  The /circ indicator is the latex interpreter for degree
    
    H2 = [h1 h2 h3];
    
    % plot other data useful for identifying events.  Here we plot y-axis
    % gyroscope (useful for identifying up/down tail beats) and jerk data
    if ~nogyro; highP = prh.Gw(:,gyroAx); else; highP = prh.Aw(:,gyroAx); end
    subplot(3,1,3);
    % first determine strokes and glides by choosing a threshold of detection
    lowpassfilt = .5; % filter out high frequency noise above 2 Hz
    [fluke_rate,q] = dsf(highP(prh.tagon),fs,lowpassfilt); %tagon = tag on 0 or 1,  fs is the sampling rate
    nf = round(4*fs/lowpassfilt) ; % how many blank spaces to make?
    fc = lowpassfilt/(fs/2) ; %specifies the cut-off frequency in Hz of a low-pass filter
    %		 to apply to X before computing the spectra. This prevents high frequency
    %		 transients e.g., in foraging, from dominating the spectra.
    filtGw = fir_nodelay(highP,nf,fc) ;
    filtGw(isnan(filtGw)) = 0;
    filtA = fir_nodelay(prh.Aw,nf,fc) ;
    [pitch,roll] = a2pr([filtA(:,1:2),-filtA(:,3)]);
    roll = -roll;

    reallylowpassfilt = 0.5*fluke_rate; % separates body posture (< than this frequency) from rotations (> than this frequency)
    fc = reallylowpassfilt; % high speed filter
    nf = round(4*fs/fc) ;
    fc = fc/(fs/2) ;
    bodyGw = fir_nodelay(filtGw,nf,fc);


%       flukeamp = nan(size(prh.p));
% %       flukeamp(prh.tagon) = wrapToPi(head2(tagon)-bodyhead(tagon));
% %     flukespeed = nan(size(p));
% %     flukespeed(tagon) = filtGw-[wrapToPi(diff(bodyhead(tagon))); 0]*fs;
    [glds,stks] = stroke_glide([filtGw(prh.tagon)-bodyGw(prh.tagon) prh.roll(prh.tagon) prh.head(prh.tagon)],fs,1,1*pi/180,15);
    stks(:,1) = round(stks(:,1)*fs+find(prh.tagon,1)-1); % convert to indices
    glds = round(glds*fs+find(prh.tagon,1)-1);
    gldI = false(size(prh.p));
    for j = 1:size(glds,1); gldI(glds(j,1):glds(j,2)) = true; end
    upstks = stks(stks(:,2) == 1,:);


    sp3 = subplot(3,1,3);
    J = njerk(prh.Aw,fs);
    I2 = I; I2(find(I,fs*10)) = false; I2(find(I,fs*10,'last')) = false; % exclude the first and last 10 seconds of jerk to exclude tag deployment and detachment
    [ax3,h1,h2] = plotyy(DN(I2),J(I2),DN,(filtGw-bodyGw)*180/pi);
    set(h1,'color','m'); % type "help plot" to see some of the color options available
    set(ax3(1),'ycolor','m');
    set(h2,'color','g','linewidth',1);
    set(ax3(2),'ycolor','g','nextplot','add','ylim',[min((filtGw-bodyGw)*180/pi) max((filtGw-bodyGw)*180/pi)]);
    h3 = plot(ax3(2),DN(upstks(:,1)),filtGw(upstks(:,1)),'rx');
    ylabel('|Jerk| (m/s^{3})');  ylabel('y-axis gyro ({\circ}/s)','parent',ax3(2));
    H3 = [h1 h2 h3];
    %
    AX = [ax1 ax2 ax3]; % put all axes in one variable to adjust all at once.  If you did not use plotyy, could have also used sp, sp2, sp3 as axes handles, but plotyy creates two new axes handles that are plotted into
    set(AX,'xlim',[DN(find(I,1)) DN(find(I,1,'last'))]); % set all xlimits to match tag on and tag off times

    linkaxes(AX,'x'); % link all the axes so if you zoom in on one, all the x-axes show the same window.
    title(prh.INFO.whaleName,'parent',ax1(1));
end
    
function L = plotmarkers(ax,DT,speedJJ,LungeI,MO,MC,A1,A2,k,filt)
    global L
%     try delete(L); catch; end
%     try oi = get(ax(2),'children'); delete(oi(1:6)); catch; end
%     L(1) = plot(ax(2),DT(LungeI),speedJJ(LungeI),'rx');
    L(1) = plot(ax(2),DT(LungeI(k)),speedJJ(LungeI(k)),'rs','markerfacecolor','r');
    if ~isnan(MC(k))
        L(2) = plot(ax(2),DT(MO(k):MC(k)),speedJJ(MO(k):MC(k)),'k','linewidth',4);
    else
        L(2) = plot(ax(2),DT(1:2),[0 0],'k','linewidth',4);
    end
    if ~isnan(A1(k))
        L(3) = plot(ax(2),DT(A1(k)),speedJJ(A1(k)),'ks','markerfacecolor','k');
    else
        L(3) = plot(ax(2),DT(1),0,'ks','markerfacecolor','k');
    end
    if ~isnan(A2(k))
        L(4) = plot(ax(2),DT(A2(k)),speedJJ(A2(k)),'ms','markerfacecolor','m');
    else
        L(4) = plot(ax(2),DT(1),0,'ms','markerfacecolor','m');
    end
     if ~isnan(filt(k))
        L(5) = plot(ax(2),DT(filt(k)),speedJJ(filt(k)),'gs','markerfacecolor','g');
    else
        L(5) = plot(ax(2),DT(1),0,'gs','markerfacecolor','g');
     end
    end