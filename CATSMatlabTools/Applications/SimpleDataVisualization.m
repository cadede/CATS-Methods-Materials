% simple data visualization

% load prh file of interest. If you have one, also load a lunge or event detection file.
% Running this cell plots the data as a continuous variable.  Use the zoom tool to zoom in to areas of interest.
% To adjust the displayed times to the zoomed in version, run (see explanation below).  Note that this line is also used at the end of this script:
%  set(AX,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
 

figure(1); clf; % chooses the figure to display
sp1 = subplot(3,1,1); % divide the figure into three rows and one column, and plot in the first plot
I = tagon; % limit the plots to only be when the tag is on the animal
try % the try/catch switch can be used in case your deployment does not have speed, in that case, the "catch" creates an empty variable
speedJJ = speed.JJ; % create a dummy variable with your speed data and smooth it with a 1 second running mean
 speedJJ(isnan(speedJJ)) = min(speedJJ); speedJJ = runmean(speedJJ,round(fs)); speedJJ(isnan(speed.JJ)) = nan;
catch; speedJJ = nan(size(p));
end
% in the first graph, plot depth (p) and speed
ax = plotyy(DN(I),p(I),DN(I),speedJJ(I)); % plotyy allows for plotting time series with different y-axes on the same x-axis
set(ax(1),'ylim',[-5 max(p(I))],'ydir','rev','nextplot','add'); % set the y-limits of the depth plot and plot so 0 is at the top, and set it so that if you plot the next line, it plots on top of the depth plot instead of replacing it
try plot(ax(1),DN(LungeI),p(LungeI),'rs','markerfacecolor','r'); % if you loaded a file with events indexed as "LungeI", plot those events as red squares on the depth graph
catch
end
ylabel('Depth (m)'); ylabel('Speed (m/s)','parent',ax(2)); % for plotting a ylabel in the non-active axis, need to identify which axis you are plotting into

% in the second graph, plot pitch, roll, heading
sp2 = subplot(3,1,2);
[ax2,h1,h2] = plotyy(DN(I),pitch(I)*180/pi,DN(I),roll(I)*180/pi); % convert radian metrics to degrees
set(h1,'color','g','linewidth',2);
set(ax2(1),'ycolor','g');
set(h2,'color','r','marker','.','markersize',4,'linestyle','none'); % make dots instead of lines for roll and head so that when either crosses -180 to 180, it does not connect the lines
set(ax2(2),'nextplot','add','ycolor','r');
plot(ax2(2),DN(I),head(I)*180/pi,'b','marker','.','markersize',4,'linestyle','none'); % plot heading as blue
ylabel('Pitch ({/circ})'); ylabel('Roll and Head ({/circ})','parent',ax2(2)); % roll and heading are on a -180 to -180 scale while pitch is on a -90 to 90 scale.  The /circ indicator is the latex interpreter for degree

% plot other data useful for identifying events.  Here we plot y-axis
% gyroscope (useful for identifying up/down tail beats) and jerk data
sp3 = subplot(3,1,3);
J = njerk(Aw,fs);
I2 = I; I2(find(I,fs*10)) = false; I2(find(I,fs*10,'last')) = false; % exclude the first and last 10 seconds of jerk to exclude tag deployment and detachment 
[ax3,h1,h2] = plotyy(DN(I2),J(I2),DN(I),Gw(I,2)*180/pi);
set(h1,'color','m'); % type "help plot" to see some of the color options available
set(ax3(1),'ycolor','m');
set(h2,'color','g','linewidth',1);
set(ax3(2),'ycolor','g');
ylabel('|Jerk| (m/s^{3})');  ylabel('y-axis gyro ({/circ}/s)','parent',ax3(2)); 

%
AX = [ax ax2 ax3]; % put all axes in one variable to adjust all at once.  If you did not use plotyy, could have also used sp, sp2, sp3 as axes handles, but plotyy creates two new axes handles that are plotted into
set(AX,'xlim',[DN(find(I,1)) DN(find(I,1,'last'))]); % set all xlimits to match tag on and tag off times

linkaxes(AX,'x'); % link all the axes so if you zoom in on one, all the x-axes show the same window.
set(AX,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS')); % can play with other time formats as well (e.g. 'mmm-dd HH:MM' will display the month, day, hour and minute).  type help datestr for more