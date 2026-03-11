[file1,loc1] = uigetfile('*.mat','Select matlab data file with most accurate timestamp');
cf = pwd; cd (loc1);
[file2,loc2] = uigetfile('*.mat','Select matlab data file to synch with IMU1');
[calfile1,calloc1] = uigetfile('*.mat','Select cal file of IMU1');
cd (calloc1);
[calfile2,calloc2] = uigetfile('*.mat','Select cal file of IMU2');
cd(cf);
IMU1 = load([loc1 file1]);
IMU2 = load([loc2 file2]);
CAL1 = load([calloc1 calfile1]);
CAL2 = load([calloc2 calfile2]);

%%
% press = true; % use pressure to synchdia

IMU1DN = IMU1.data.Date+IMU1.data.Time;
IMU2DN = IMU2.data.Date+IMU2.data.Time + IMU1.Hzs.UTC/24 - IMU2.Hzs.UTC/24; % converts to UTC and then to offset of IMU1
IMU1ADN = IMU1.Atime;
IMU2ADN = IMU2.Atime + IMU1.Hzs.UTC/24 - IMU2.Hzs.UTC/24;
try fs2 = IMU2.Hzs.datafs; catch; fs2 = IMU2.Hzs.accHz; end
try fs1 = IMU1.Hzs.datafs; catch; fs1 = IMU1.Hzs.accHz; end
% [~,A1] = applyCal2(IMU1.data,IMU1DN,CAL1,true(size(IMU1DN)),fs1,IMU1.Hzs);
% [~,A2] = applyCal2(IMU2.data,IMU2DN,CAL2,true(size(IMU2DN)),fs2,IMU2.Hzs);

A1 = IMU1.Adata;
A2 = IMU2.Adata;
% A1 = [IMU1.data.Acc1 IMU1.data.Acc2 IMU1.data.Acc3];
% A2 = [IMU2.data.Acc1 IMU2.data.Acc2 IMU2.data.Acc3];


acal = CAL1.acal; aconst = CAL1.aconst; try Acal = CAL1.Acal3d0; catch; end
if exist('Acal','var') && ~isempty(Acal)
    axA = (acal./abs(acal)); axA(isnan(axA)) = 0;
    A1 = A1*axA;
    A1 = (A1*diag(Acal.poly(:,1))+repmat(Acal.poly(:,2)',size(A1,1),1))*Acal.cross;
else
    A1 = (A1-repmat(aconst,size(A1,1),1))*acal;
end
clear Acal;
acal = CAL2.acal; aconst = CAL2.aconst; try Acal = CAL2.Acal3d0; catch; end
if exist('Acal','var') && ~isempty(Acal)
    axA = (acal./abs(acal)); axA(isnan(axA)) = 0;
    A2 = A2*axA;
    A2 = (A2*diag(Acal.poly(:,1))+repmat(Acal.poly(:,2)',size(A2,1),1))*Acal.cross;
else
    A2 = (A2-repmat(aconst,size(A2,1),1))*acal;
end
clear Acal;

% if press
IMU1p = runmean(edgenans(fixgaps(IMU1.data.Pressure)),fs1/2);
IMU2p = runmean(edgenans(fixgaps(IMU2.data.Pressure)),fs2/2);
% else
%     IMU1p = A1(:,3); %IMU1.data.Light2;
%     IMU2p = A2(:,3); %IMU2.data.Light2;
% end
% IMU1T = IMU1.data.Temp; %IMU1T2 = IMU1.data.Temp1;
% IMU1Tf = fir_nodelay(IMU1T,500, .1/50/2);
% % IMU1Tf2 = fir_nodelay(IMU1T2,500, .1/50/2);
% IMU2T = IMU2.data.Temp; %IMU2T2 = IMU2.data.Temp1;
% IMU2Tf = fir_nodelay(IMU2T,500, .1/50/2);
% IMU2Tf2 = fir_nodelay(IMU2T2,500, .1/50/2);
% A2 = IMU2.Adata; %[IMU2.data.Acc1 IMU2.data.Acc2 IMU2.data.Acc3];
% A1 = IMU1.Adata; %[IMU1.data.Acc1 IMU1.data.Acc2 IMU1.data.Acc3];


figure(4); clf;
subplot(211);
axA = plotyy(IMU1ADN,A1,IMU1DN,-IMU1p);
set(axA,'nextplot','add');
plot(axA(1),IMU2ADN,A2-2);
plot(axA(2),IMU2DN,-IMU2p)
legend('IMU1X','IMU1Y','IMU1Z','IMU2X','IMU2Y','IMU2Z','p1','p2')
axT = subplot(212);
plot(IMU1DN,IMU1.data.Temp,IMU2DN,IMU2.data.Temp+1);
linkaxes([axA axT],'x');
legend('IMU1 Temp', 'IMU2 Temp');
set([axA axT],'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS.fff'))
%%
ob = 0; xs = nan(0,2); k = 1;
while ob~=113
    figure(4);
    title('select window to look for synch points, or press q twice if satisfied with line')
    [x,~,ob] = ginput(2);
    if any(ob == 113); break; end
    [~,a] = min(abs(IMU1DN-x(1)));
    [~,b] = min(abs(IMU1DN-x(2)));
    I1 = a:b;
    [~,a] = min(abs(IMU2DN-x(1)));
    [~,b] = min(abs(IMU2DN-x(2)));
    I2 = a:b;
    [~,a] = min(abs(IMU1ADN-x(1)));
    [~,b] = min(abs(IMU1ADN-x(2)));
    AI1 = a:b;
    [~,a] = min(abs(IMU2ADN-x(1)));
    [~,b] = min(abs(IMU2ADN-x(2)));
    AI2 = a:b;
    
    % I1 = 1:3000000;
    figure(5); clf;
    subplot(211);
    [axA,h1,h1p] = plotyy(IMU1ADN(AI1),A1(AI1,:),IMU1DN(I1),-IMU1p(I1));
    set(axA,'nextplot','add');
    h2 = plot(axA(1),IMU2ADN(AI2),A2(AI2,:)-2,'--');
    h2p = plot(axA(2),IMU2DN(I2),-IMU2p(I2));
    c1 = get(h1,'color');
    for ii = 1:3; set(h2(ii),'color',c1{ii}); end
    legend('IMU1X','IMU1Y','IMU1Z','IMU2X','IMU2Y','IMU2Z','p1','p2')
    %  title ('zoom to section of interest')
    % axT = subplot(212);
    % plot(IMU1DN(I1),IMU1.data.Temp(I1),IMU2DN(I2),IMU2.data.Temp(I2)+1);
    % linkaxes([axA axT],'x');
    % legend('IMU1 Temp', 'IMU2 Temp');
    set(axA,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS.fff'))
    afs2 = IMU2.Hzs.accHz;
    afs1 = IMU1.Hzs.accHz;
    % try delete(axT); catch; end;
    % for i = 2
    button = 0;
    %     figure(i);
    %     t = title('Set y-axis zoom for IMU1 (press enter when okay');
    
    
    %     if i == 1; ax = axp; YC = IMU1p; YD = IMU2p; else ax = axA; YC = max(IMU1Tf2,[],2); YD = max(IMU2Tf2,[],2); end
    ax = axA;
    c1 = get(h1,'color');
    c2 = get(h2,'color');
    x2 = get(ax(2),'xlim');
    y1 = get(ax(1),'ylim'); x1 = get(ax(1),'xlim');
    %     linkaxes(ax,'x');    linkaxes(ax,'x');
    axes(ax(1))
    %     pause;
    y2 = get(ax(2),'ylim');% delete(t);
    set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
    %
    while button~=113
        set(ax(1),'xlim',x1); set(ax(1),'ylim',y1);
        try delete(t); catch; end
        try delete(cc); catch; end
        %         set(ax,'xlim',x2); set(ax(2),'ylim',y2); set(ax(1),'ylim',y1);
        for ii = 1:3; set(h2(ii),'color',c2{ii});
            set(h1(ii),'color',c1{ii});
        end
        t = title('zoom to desired level, then press enter');
        pause; delete(t)
        set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
        t = title('Line up cursor with peak on IMU2, then press x, y, or z for that acc axis peak. Press q to quit or r to redo zoom.');
        %         t = title('Click peak for IMU2 then for IMU1 (or q twice if done). Left click maxima, right click minima.');
        [x,~,button] = ginput(1);
        %         if button == 1; ext = 1; else ext = -1; end
        ext = 1; %could switch to look for minima
        if ~ismember(button,120:122); continue; end
        [~,Id] = min(abs(IMU2ADN-x(1)));
        YD = A2(:,button-119);
        [~,peaklocD] = max(YD(Id-.2*afs2:Id+.2*afs2)*ext);
        peaklocD = round(peaklocD+Id-.2*afs2-1);
        cc = plot(ax(1),IMU2ADN(peaklocD),A2(peaklocD,button-119)-2,'o');
        oi = 1:3;
        oi = oi(120:122 ~= button);
        set(h2(oi),'color','w');
        set(h1(oi),'color','w');
        delete(t)
        t = title('Line up cursor with peak on IMU1, then click for that acc axis peak. Press c to cancel point 1.');
        [x,~,button2] = ginput(1);
        if button2 == 99; continue; end
        [~,Ic] = min(abs(IMU1ADN-x(1)));
        YC = A1(:,button-119);
        [~,peaklocC] = max(YC(Ic-.2*afs1:Ic+.2*afs1)*ext);
        peaklocC = round(peaklocC+Ic-.2*afs1-1);
        plot(ax(1),IMU1ADN(peaklocC),A1(peaklocC,button-119),'o');
        xs(k,:) = [IMU2ADN(peaklocD)  IMU1ADN(peaklocC)];
        s2 = subplot(212); hold on;
        %                 plot(IMU2ADN(peaklocD), diff(xs(k,:))*24*60*60,'o');
        plot(xs(:,1), diff(xs,[],2)*24*60*60,'o');
        offsets = diff(xs,[],2)*24*60*60;
        xlim([min(xs(:,1))-5/24/60/60 max(xs(:,1))+5/24/60/60])
        ylim([min(offsets)-2 max(offsets)+2]);
        xlabel('Time');
        ylabel('Offset (s)');
        set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
        k = k+1;
        %         x = [IMU2DN(round(peaklocD+Id-.2*fs2-1)) IMU1DN(round(peaklocC+Ic-.2*fs1-1)) ];
        %         if button ~= 113
        %             xs(k,:) = x'*24*60*60; k = k+1;
        %         end
    end
end
%%
figure(3);% if i == 1; clf; end;
try delete(s3(1)); catch; end
s3(1) = subplot(2,1,1);
%     plot(xs(:,2)/60/60-IMU2DN(1)*24,diff(xs,1,2),'o');
plot(xs(:,1)*24-IMU2ADN(1)*24,diff(xs,1,2)*24*60*60,'o')
%     [p,S] = polyfit(xs(:,2)/60/60-IMU2DN(1)*24,diff(xs,1,2),1);
[p,S] = polyfit(xs(:,1)*24-IMU2ADN(1)*24,diff(xs,1,2)*24*60*60,1);
set(s3(1),'nextplot','add')
plot(s3(1),get(s3(1),'xlim'),get(s3(1),'xlim')*p(1)+p(2),'linewidth',2);
xlabel('Hours since start of IMU2');
ylabel('Seconds of offset from IMU1era');

%     figure(i); set(gca,'nextplot','add');
%     FSnew = round(40000*(1+p(1)/60/60)); %*fs2*750
%     FSold = 40000;
%
%     IMU2pnew = resample(IMU2.data.Pressure,FSnew,FSold); % be careful as old resample in cats tools uses timeseries but that is not necessary any more
%     IMU2DNnew = ((0:length(IMU2pnew)-1)/fs2/24/60/60+IMU2DN(1) + p(2)/24/60/60)';

if abs(p(1)*length(IMU2DN)/fs2/60/60)>1 % if more than 1 second loss over the deployment
    oldDN = 0:length(IMU2DN)-1;
    newDN = 0:(1-p(1)/60/60):max(oldDN);
    p2 = timeseries(IMU2.data.Pressure,oldDN);
    IMU2pnew = resample(p2,newDN);
    IMU2pnew =IMU2pnew.data;
    IMU2DNnew = ((0:length(IMU2pnew)-1)/fs2/24/60/60+IMU2DN(1) + p(2)/24/60/60)';
else
    IMU2pnew = IMU2p;
    IMU2DNnew = IMU2DNnew+p(2)/60/60/24;
    
    
end

%     if ~press
%     IMU1p = IMU1.data.Pressure;
%     IMU2p = nan(size(IMU2DNnew));
% %     IMU2DN = IMU2.data.Date+IMU2.data.Time+TimeOffset_s/24/60/60;
%     [~,k1] = min(abs(IMU1DN-IMU2DNnew(1))); [~,k2] = min(abs(IMU1DN-IMU2DNnew(end)));
%     [~,kd1] = min(abs(IMU2DNnew-IMU1DN(k1))); [~,kd2] = min(abs(IMU2DNnew-IMU1DN(k2)));
%     if size(IMU2p(kd1:IMU2.Hzs.datafs/IMU1.Hzs.datafs:kd2),1)+1 == length(k1:k2)
%         k2 = k2-1;
%     end
%     IMU2p(kd1:IMU2.Hzs.datafs/IMU1.Hzs.datafs:kd2) = IMU1p(k1:k2);
%      IMU2p = fixgaps(IMU2p); IMU2p(isnan(IMU2p)) = 0;
%       IMU2pnew = IMU2p;%resample(IMU2p,FSnew,FSold);
%     end

%     newIMU2DN = (((IMU2DN-IMU2DN(1))*24)*p(1)+p(2))/24+IMU2DN(1);
figure(6); clf;
plot(IMU1DN,IMU1p,IMU2DN,IMU2p,IMU2DNnew,IMU2pnew,'g--');
legend('IMU1 p','IMU2 p','New IMU2 p');

% disp(num2str(diff(x)*24*60*60));
% disp(num2str(diff(x2)*24*60*60));
%%
IMU2New = IMU2;
IMU2New.data = IMU2.data;
IMU2New.data(length(IMU2DNnew)+1:end,:) = [];
IMU2New.data = [IMU2New.data; IMU2New.data(1:length(IMU2DNnew)-size(IMU2.data,1),:)]; % makes it longer if necessary
IMU2New.data.Date = floor(IMU2DNnew); IMU2New.data.Time = IMU2DNnew-floor(IMU2DNnew);
if abs(p(1)*length(IMU2DN)/fs2/60/60)>1 % if more than 1 second loss over the deployment
    
    for k = 3:size(IMU2.data,2)
        % IMU2New.data{:,k} = resample(IMU2.data{:,k},FSnew,FSold);
        p2 = timeseries(IMU2.data{:,k},oldDN);
        IMU2pnew = resample(p2,newDN);
        IMU2New.data{:,k} =IMU2pnew.data;
    end
    
    % IMU2New.Adata = resample(IMU2New.Adata,FSnew,FSold);
    oldDN = 0:length(IMU2ADN)-1;
    newDN = 0:(1-p(1)/60/60):max(oldDN);
    p2 = timeseries(IMU2New.Adata,oldDN);
    IMU2pnew = resample(p2,newDN);
    IMU2New.Adata =IMU2pnew.data;
end
IMU2New.TimeOffset_s = p(2);
IMU2New.TimeSlope_s_per_h = p(1);
% IMU2New.TimeOffsetNotes = 'slope was .1372, but changed to 0 since drift evidence was not strong'
ODN = IMU2New.data.Date(1)+IMU2New.data.Time(1);
IMU2New.ODN = ODN;
IMU2New.Atime = ((0:size(IMU2New.Adata,1)-1)/IMU2New.Hzs.accHz/24/60/60+ODN)';
figure; plotyy(IMU2New.data.Date+IMU2New.data.Time,IMU2New.data.Pressure,IMU2New.data.Date+IMU2New.data.Time,[IMU2New.data.Acc1 IMU2New.data.Acc2 IMU2New.data.Acc3])
%%
copyfile([loc2 file2],[loc2 file2(1:end-4) 'b4resample.mat']);
saveas(3,[loc2 'timeOffsetCal.jpg'])
saveas(6,[loc2 'timeOffsetResult.fig'])

try
   save([loc2 file2],'-struct','IMU2New');
    if ~isempty(lastwarn)
        error(lastwarn);
    end
    disp('File Creation successful!' );
    disp(newfileloc)
catch %v7.3 allows for bigger files, but makes a freaking huge file if used when you don't need it
 save([loc2 file2],'-struct','IMU2New','-v7.3');
       disp('Made a version 7.3 file in order to include all vars');
    disp('File Creation successful!' );
    disp(newfileloc)
end
%% Alternate version if just time synch required (old, untested)
% use this to choose an offset after zooming in to acc tag on periods (for D1 with no drift issue)
% (not relevant since all seem to have drift). Skip this section and next
title('Click on IMU2 tag on spike, then on IMU1 tag on spike');
[x,~] = ginput(2);
TimeOffset_s = round(diff(x)*24*60*60*IMU2.Hzs.datafs)/IMU2.Hzs.datafs;
TimeSlope_s_per_h = 0;
TimeOffsetNotes = 'Only Offset calculated, none applied, do doing prh creation';
IMU1p = IMU1.data.Pressure;
IMU2p = nan(size(IMU2.data.Pressure));
IMU2DN = IMU2.data.Date+IMU2.data.Time+TimeOffset_s/24/60/60;
[~,k1] = min(abs(IMU1DN-IMU2DN(1))); [~,k2] = min(abs(IMU1DN-IMU2DN(end)));
[~,kd1] = min(abs(IMU2DN-IMU1DN(k1))); [~,kd2] = min(abs(IMU2DN-IMU1DN(k2)));
if size(IMU2p(kd1:IMU2.Hzs.datafs/IMU1.Hzs.datafs:kd2),1)+1 == length(k1:k2)
    k2 = k2-1;
end
IMU2p(kd1:IMU2.Hzs.datafs/IMU1.Hzs.datafs:kd2) = IMU1p(k1:k2);

%     for k = k1:k2
%         [~,kk] = min(abs(IMU2DN-IMU1DN(k)));
%         IMU2p(kk) = IMU1p(k);
%     end
IMU2p = fixgaps(IMU2p); IMU2p(isnan(IMU2p)) = 0;
figure(5); clf;
s1 = subplot(311);
plot(IMU2DN,IMU2p,IMU1DN,IMU1p,'--');
s2 = subplot(312);
plot(IMU2DN,IMU2.data.Gyr3);% [IMU2.data.Comp1 IMU2.data.Comp2 IMU2.data.Comp3]);
s3 = subplot(313);
plot(IMU1DN,A1);
linkaxes([s1 s2 s3],'x')
%%
data = IMU2.data; data.Pressure = IMU2p;%nan(size(IMU2.data.Pressure));


save([fold ID '\D' tagnum '\' IMU2file],'data','TimeOffset_s','TimeSlope_s_per_h','TimeOffsetNotes','-append');
