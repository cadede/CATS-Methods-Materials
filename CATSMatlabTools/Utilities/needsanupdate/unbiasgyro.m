function [GwUB, bias, pitchgy,rollgy,headgy] = unbiasgyro(Aw,fs,Gw,Mw,pitch,roll,head,tagondec,camondec)

% Should be mostly automatic, but will take a long time (an hour?).  Each
% line displayed shows the index, so you should know how much longer based
% on how long your files are (i.e. the number of rows in At, Mt, pitch
% etc.)
% output: pitchgy, rollgy, headgy
J = njerk(Aw,fs);
At_mag = sqrt(sum(Aw.^2,2));

njerkM = runmean(J,round(fs));
% next line is untested
bigA = find(njerkM>2*median(njerkM) | abs(pitch)*180/pi>70); %a threshold for times when specific acceleration is likely high or when gimbal lock is likely
% At_magM = runmean(At_mag,round(fs));
% bigA =
tagonI = 1:length(tagondec); tagonI = tagonI(tagondec);
bigA = intersect(bigA,tagonI(round(3*fs):end)); % start looking 3 seconds after tag on
[s,e] = consec(bigA);
I = find(e-s>=fs); %find the spots where the njerk (or gimbal lock) is big for longer than 1 second at a time
s = sort(s(I)); e = sort(e(I)); % may be unnecessary?  Was created for testing by mixing in another value
oi = find(abs(e(1:end-1)-s(2:end))<2*fs); % if there are gaps of <2 seconds where njerkM is low, gyro through them.
if ~isempty(oi); s(oi+1) = []; e(oi) = []; end
%  e(31) = 46300;
Awgy = Aw; pitchgy = pitch; rollgy = roll; headgy = head; pitchgy2 = pitch; rollgy2 = roll; headgy2 = head; head2gy = head;
At_maggy = At_mag; Mhgy = Mh;
% oi = [Gw(:,1) Gw(:,2) Gw(:,3)];
tic
% first go through and add a bias to Gw that would make the gyro turns give you the next location.  When ac is high, this will be overwritten below
GwUB = Gw; % Gw unbiased
bias = zeros(size(Gw));
%we assume that the three gyro rotations at a given time best describe the
%rotation to the next position in time.  This was supported in one example
%with smaller biases doing it that way

for i = 1:length(Gw(:,1))-1
    if any(isnan([head(i:i+1); pitch(i:i+1); roll(i:i+1)])); continue; end
    quat1 = angle2quat(head(i),pitch(i),roll(i));
    quat2 = angle2quat(head(i+1),pitch(i+1),roll(i+1));
    vv=Gw(i,:)'/fs;
    magvv=sqrt(vv'*vv);
    qv=[cos(magvv/2);(1/magvv)*sin(magvv/2)*vv]'; %cos(angle of rotation/2), sin(angle/2)*axis of rotation.  Axis and angle from Tomaži?, S.; Stan?in, S. Simultaneous orthogonal rotation angle. Electrotech. Rev. 2011, 78, 7-11
    qbias = quatmultiply(quatinv(quatmultiply(quat1,qv)),quat2);
    magBias = 2*acos(qbias(1));
    bias(i,:) = qbias(2:4)*magBias/sin(magBias/2)*fs;
    GwUB(i,:) = Gw(i,:)+bias(i,:);
end
bias2 = bias;
for i = 2:length(Gw(:,1))
    if any(isnan([head(i-1:i); pitch(i-1:i); roll(i-1:i)])); continue; end
    quat1 = angle2quat(head(i-1),pitch(i-1),roll(i-1));
    quat2 = angle2quat(head(i),pitch(i),roll(i));
    vv=Gw(i,:)'/fs;
    magvv=sqrt(vv'*vv);
    qv=[cos(magvv/2);(1/magvv)*sin(magvv/2)*vv]'; %cos(angle of rotation/2), sin(angle/2)*axis of rotation.  Axis and angle from Tomaži?, S.; Stan?in, S. Simultaneous orthogonal rotation angle. Electrotech. Rev. 2011, 78, 7-11
    qbias = quatmultiply(quatinv(quatmultiply(quat1,qv)),quat2);
    magBias = 2*acos(qbias(1));
    bias2(i,:) = qbias(2:4)*magBias/sin(magBias/2)*fs;
end
Ise = 1:length(Gw);
for i = 1:length(s)
    Ise(s(i):e(i)) = true;
end
Ise(1) = 0;
Ise(Ise~=1) = false;
sum(sum(bias(tagondec&~Ise',:).^2))
sum(sum(bias2(tagondec&~Ise',:).^2)) %this one should be bigger based on past experience, but could be random

for i = 1:length(s) %[6 29:32 95]%
    [~,I] = min(abs(At_mag(s(i)-round(2*fs):s(i))-1));
    s(i) = s(i)-round(2*fs)-1+I; %find the place up to two seconds before high jerk where acceleration is closest to magnitude of 1.  This will be the baseline
    try [~,I] = min(abs(At_mag(e(i):e(i)+round(2*fs))-1)); catch err; if i==length(s); [~,I] = min(abs(At_mag(e(i):length(At_mag))-1)); else throw (err); end; end
    e(i) = e(i)-1+I;
    I = find(~isnan(head(1:s(i))),1,'last'); %ensure that the heading exists (since it's bad where the cam turns on or off)
    if abs(s(i) - I) > 5*fs; I = find(~isnan(head(s(i):end)),1,'first')+s(i)-1; end
    if abs(s(i)-I) <=5*fs; if ~isnan(head(I-1)); s(i) = I; else s(i) = I+1; end
    else continue; end %if you can't start where a heading exists (within 5 sec), just give up
    I = find(~isnan(head(e(i):end)),1,'first') +e(i) - 1;
    if abs(e(i) - I) > 3*fs; I = find(~isnan(head(1:e(i))),1,'last'); end
    if abs(e(i)-I) <= 3*fs; e(i) = I; end % if the end is too far, just use what you've got and hope the pitch/roll matching is good enough.
    if e(i)<=s(i); continue; end; % if the values you could use don't work, just skip this one (new edit)
    bothways = false; % does the (unadjusted) gyros going backwards as well for half the length of the gyro time
    [rg, pg, hg, R,h] = quatGyros([Gw(:,1) Gw(:,2) Gw(:,3)],s(i),e(i),roll,pitch,head,fs,tagslipdec,njerkM,bothways); %uses quaternions to keep track of the gyros.  adjust them all by a single small rotation for each period so that the gyros end up at the known end points smoothly.  If there is a tagslip or the adjustment is larger than 1/3 degree per turn or you are more than 1 degree off at the end, look for a tag slip region and use two different calibrations.
    pitchgy(s(i):e(i)) = pg; rollgy(s(i):e(i)) = rg; headgy(s(i):e(i)) = hg;
    if ~bothways; endpoint = e(i); else endpoint = round((e(i)-s(i))/2)+s(i); end
    for j = s(i):endpoint; Awgy(j,:) = Awgy(j-1,:)*R{j-s(i)+1}; end
    if bothways
        for j = e(i):-1:endpoint+1; Awgy(j,:) = Awgy(j+1,:)*R{j-s(i)+2}'; end
    end
    
    At_maggy(s(i):e(i)) = sqrt(sum(Awgy(s(i):e(i),:).^2,2));
    pitchgy2(s(i):e(i))=asin(Awgy(s(i):e(i),1)./At_maggy(s(i):e(i)));
    rollgy2(s(i):e(i))=atan2(-Awgy(s(i):e(i),2),-Awgy(s(i):e(i),3));
    for j=s(i):e(i)
        P = [cos(pitchgy2(j)) 0 sin(pitchgy2(j)); 0 1 0; -sin(pitchgy2(j)) 0 cos(pitchgy2(j))];
        R = [1 0 0; 0 cos(rollgy2(j)) -sin(rollgy2(j)); 0 sin(rollgy2(j)) cos(rollgy2(j))];
        Mhgy(j,:) = Mw(j,:)*R'*P';
        headgy2(j) = wrapToPi(atan2(-Mhgy(j,2),Mhgy(j,1)));% add dec back in to get the true north bearing
        if headgy2(j) >0; head2gy(j) = headgy2(j)-pi; else head2gy(j) = headgy2(j)+pi; end
    end
    if length(h)==4; endpoint = e(i); else endpoint = h(end); end
    qbias = quatnormalize(h(1:4));
    magBias = 2*acos(qbias(1));
    bias(s(i):endpoint,:) = repmat(qbias(2:4)*magBias/sin(magBias/2),endpoint-s(i)+1,1)*fs;
    % create an adjusted whale-frame Gyro where the bias is included in each
    % value
    GwUB(s(i):endpoint,:) = Gw(s(i):endpoint,:)+bias(s(i):endpoint,:);
    if length(h)>4 %if length>4, you had two biases on each side of a tagslip
        qbias = quatnormalize(h(5:8));
        magBias = 2*acos(qbias(1));
        bias(endpoint+1:e(i),:) = repmat(qbias(2:4)*magBias/sin(magBias/2),e(i)-endpoint,1)*fs;
        GwUB(endpoint+1:e(i),:) = Gw(endpoint+1:e(i),:)+bias(endpoint+1:e(i),:);
    end
end
toc

highJerk = [s;e];

figure(8); clf;
set(gcf,'windowStyle','docked');
sp(1) = subplot(8,1,1:7); hold on;
plot([pitchgy2 rollgy2 headgy2 pitchgy rollgy headgy]*180/pi,'--')
plot([pitch roll head]*180/pi,'.','markersize',4);

title('Using quaternion vectors from current orientation')
legend('Gyp','Gyr','Gyh','GypAdj','GyrAdj','GyhAdj','p','r','h');
sp(2) = subplot(8,1,8);
ax = plotyy(1:length(Depth),njerk,1:length(Depth),Depth);
set(ax(2),'ydir','rev','ylim',[0 200],'nextplot','add');
plot(ax(2),camondec*100,'k','linewidth',3)
set(ax(1),'ylim',[0 max(njerkM(tagondec))]);
