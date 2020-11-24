% makes some graphs so you can see progress.  Don't need to check necessarily unless something looks off,
% At this point everything up to what DTags have is completed (also all what you would need for trackplot (except speed perhaps in new verion of trackplot))
figure(4); clf;
set(gcf,'windowStyle','docked');
s1 = subplot(3,1,1);
plot(DN, At); grid on
title('"TAG" Acceleration (g) via TAG FRAME vs Time');
xlim([DN(1) DN(end)]); ylim([-2 2]);
set(gca,'xticklabel',[]);
s2 = subplot(3,1,2);
plot(DN, Aw); grid on
xlim([DN(1) DN(end)]); ylim([-2 2]);
legend('X','Y','Z','orientation','horizontal');
set(gca,'xticklabel',[]);
title('"ANIMAL" Acceleration (g) via ANIMAL FRAME vs Time')
s3 = subplot(3,1,3);
if ~nopress
    plot(DN, Depth); grid on; title('Animal Depth(m) vs time (s)');set(gca,'YDir','rev'); ylim([0 max(Depth)]);
else plot(DN,Light);
end

xlim([DN(1) DN(end)]);
oi = datestr(get(gca,'xtick')); oi = oi(:,13:end);%[oi(:,1:5) repmat(10,size(oi,1),1) oi(:,6:end)];
set(gca,'xticklabel',oi);
linkaxes([s1 s2 s3],'x');
hold on;  pp1 = plot(repmat(DN(Wchange),1,2)',repmat([-300 300]',1,length(Wchange)),'k','linewidth',3); pp2 = plot(repmat(DN(Wchangeend),1,2)',repmat([-300 300]',1,length(Wchangeend)),'r'); %pp3 = plot(repmat(DN(startsI),1,2)',repmat([-300 300]',1,length(startsI)),'g--');
if ~isempty(Wchange); legend([pp1(1) pp2(1)],'TagslipStart','TagslipStop','location','southeast','orientation','horizontal'); end
% ANIMAL FRAME PITCH, ROLL, and HEADING CALCULATIONS

% 4) PITCH and ROLL
% Plot Animal Frame Pitch and Roll
figure(5); clf;
set(gcf,'windowStyle','docked');
ss1 = subplot(3,1,1);
plot(DN,pitch*180/pi);
xlim([DN(1) DN(end)]);
grid on
title('Pitch Animal (degrees) vs time (s)'); ylabel ('Pitch Animal(degrees)')
ss2 = subplot(3,1,2);
plot(DN,roll*180/pi);%,DN,roll2*180/pi);
title('Roll Animal (degrees) vs. time (s)'); ylabel ('Roll Animal (degrees)');
legend('Roll');
xlim([DN(1) DN(end)]);
grid on
ss3 = subplot(3,1,3);
if ~nopress
    plot(DN, Depth); grid on; title('Depth (m) vs Time (s)'); set(gca,'YDir','rev')
else plot(DN,Light);
end

grid on
xlim([DN(1) DN(end)]);
oi = datestr(get(gca,'xtick')); oi = oi(:,13:end);%[oi(:,1:5) repmat(10,size(oi,1),1) oi(:,6:end)];
set(gca,'xticklabel',oi);
linkaxes([s1 s2 s3 ss1 ss2 ss3],'x');
hold on;  pp1 = plot(repmat(DN(Wchange),1,2)',repmat([-300 300]',1,length(Wchange)),'k','linewidth',3); pp2 = plot(repmat(DN(Wchangeend),1,2)',repmat([-300 300]',1,length(Wchangeend)),'r'); %pp3 = plot(repmat(DN(startsI),1,2)',repmat([-300 300]',1,length(startsI)),'g--');
if ~isempty(Wchange); legend([pp1(1) pp2(1)],'TagslipStart','TagslipStop','location','southeast','orientation','horizontal'); end

% Plot Calibrated Gyro Data
figure(6); clf;
set(gcf,'windowStyle','docked');
subplot(3,1,1)
plot(DN, Gw);
legend('X','Y','Z','orientation','horizontal');
xlim([DN(1) DN(end)]);
oi = datestr(get(gca,'xtick')); oi = oi(:,13:end);%[oi(:,1:5) repmat(10,size(oi,1),1) oi(:,6:end)];
set(gca,'xticklabel',oi);

grid on
title('GYRO Animal vs time (s)'); ylabel ('Gryo Animal')

% Plot the True Heading for Methods 1,2 and Depth vs Time.
figure(7); clf;
set(gcf,'windowStyle','docked');
sss1 = subplot(2,1,1);
head(abs(diff(head)) > 350*pi/180) = nan;
plot(DN,head*180/pi); title({'True Heading: Inspect heading, sometimes rapid heading changes'; 'can be a good indicator of where a tag slip may need a fine scale adjustment'}); ylabel ('True Heading Degrees (degrees)'); xlabel('Time ')
xlim([DN(1) DN(end)]);
oi = datestr(get(gca,'xtick'),'HH:MM:SS'); set(gca,'xticklabel',oi);
grid on
sss2 = subplot(2,1,2);
plot(DN,Depth); title('Depth(m) vs Time (s) '); ylabel('Depth (m)'); xlabel('Time (s) ')
set(gca,'YDir','rev')
grid on
xlim([DN(1) DN(end)]);
oi = datestr(get(gca,'xtick')); oi = oi(:,13:end);%[oi(:,1:5) repmat(10,size(oi,1),1) oi(:,6:end)];
set(gca,'xticklabel',oi);
linkaxes([s1 s2 s3 ss1 ss2 ss3 sss1 sss2],'x');
hold on;  pp1 = plot(repmat(DN(Wchange),1,2)',repmat([-300 300]',1,length(Wchange)),'k','linewidth',3); pp2 = plot(repmat(DN(Wchangeend),1,2)',repmat([-300 300]',1,length(Wchangeend)),'r'); %pp3 = plot(repmat(DN(startsI),1,2)',repmat([-300 300]',1,length(startsI)),'g--');
if ~isempty(Wchange); legend([pp1(1) pp2(1)],'TagslipStart','TagslipStop','location','southeast','orientation','horizontal'); end
