[file1,loc1] = uigetfile('*.mat','Select file with most accurate timestamp');
[file2,loc2] = uigetfile('*.mat','Select file to synch with IMU1');

IMU1 = load([loc1 file1]);
IMU2 = load([loc2 file2]);

%%
press = true; % use pressure to synchdia

IMU1DN = IMU1.data.Date+IMU1.data.Time;
IMU2DN = IMU2.data.Date+IMU2.data.Time + IMU1.Hzs.UTC/24 - IMU2.Hzs.UTC/24;
IMU1ADN = IMU1.Atime;
IMU2ADN = IMU2.Atime + IMU1.Hzs.UTC/24 - IMU2.Hzs.UTC/24;
try fs2 = IMU2.Hzs.datafs; catch; fs2 = IMU2.Hzs.accHz; end
try fs1 = IMU1.Hzs.datafs; catch; fs1 = IMU1.Hzs.accHz; end
% if press
    IMU1p = runmean(edgnans(fixgaps(IMU1.data.Pressure)),fs1*2);
    IMU2p = runmean(edgenans(fixgaps(IMU2.data.Pressure)),fs2*2);
% else
%     IMU1p = IMU1.data.Light2;
%     IMU2p = IMU2.data.Light2;
% end
IMU1T = IMU1.data.Temp; %IMU1T2 = IMU1.data.Temp1;
IMU1Tf = fir_nodelay(IMU1T,500, .1/50/2);
% IMU1Tf2 = fir_nodelay(IMU1T2,500, .1/50/2);
IMU2T = IMU2.data.Temp; %IMU2T2 = IMU2.data.Temp1;
IMU2Tf = fir_nodelay(IMU2T,500, .1/50/2);
% IMU2Tf2 = fir_nodelay(IMU2T2,500, .1/50/2);
A2 = IMU2.Adata; %[IMU2.data.Acc1 IMU2.data.Acc2 IMU2.data.Acc3];
A1 = IMU1.Adata; %[IMU1.data.Acc1 IMU1.data.Acc2 IMU1.data.Acc3];


figure(4); clf;
subplot(211);
axA = plotyy(IMU1ADN,A1,IMU1DN,-IMU1.data.Pressure);
set(axA,'nextplot','add');
plot(axA(1),IMU2ADN,A2-20);
plot(axA(2),IMU2DN,-IMU2.data.Pressure)
 legend('IMU1X','IMU1Y','IMU1Z','IMU2X','IMU2Y','IMU2Z','p1','p2')
axT = subplot(212);
plot(IMU1DN,IMU1.data.Temp,IMU2DN,IMU2.data.Temp+1);
linkaxes([axA axT],'x');
legend('IMU1 Temp', 'IMU2 Temp');
set([axA axT],'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS.fff'))

%% run this regardless
% IMU1p = IMU1Tf2; IMU2p = IMU2Tf2;
figure (1); clf; axp = plotyy(IMU1DN,IMU1p,IMU2DN,IMU2p);
legend('IMU1 p','IMU2 p');
% IMU1Tf = IMU1p; IMU1Tf2 = A1(:,:); IMU2Tf = IMU2p; IMU2Tf2 = A2(:,:);
figure(2); clf; axt = plotyy(IMU1DN,IMU1p,IMU2DN,IMU2p); set(axt,'nextplot','add');
plot(axt(1),IMU1ADN,A1,'linewidth',1); plot(axt(1),IMU2ADN,A2-20,'linewidth',1);
% legend('IMU1 p','IMU1 A3','IMU2 p','IMU2 A3');
 legend('cp','ac1','ac2','ac3','cp','ad1','ad2','ad3')

set(axp,'ylim',[0 max([IMU1p; IMU2p])]);
set(axt,'ylim',[-40 max([IMU1p; IMU2p])])
% figure(3); clf


%% set i = 1 if both IMU1 and IMU2 have pressure, else set i = 2 if only IMU1 has pressure
for i = 1
    xs = nan(0,2); k = 1; button = 0;
    figure(i);
    t = title('Set y-axis zoom for IMU1 (press enter when okay');
    

    if i == 1; ax = axp; YC = IMU1p; YD = IMU2p; else ax = axt; YC = max(IMU1Tf2,[],2); YD = max(IMU2Tf2,[],2); end
    x2 = get(ax(2),'xlim');
    y1 = get(ax(1),'ylim'); x1 = get(ax(1),'xlim');
    linkaxes(ax,'x');    linkaxes(ax,'x');
    axes(ax(2))
    pause; 
    y2 = get(ax(2),'ylim'); delete(t);
    set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
    while button~=113
        delete(t);
%         set(ax,'xlim',x2); set(ax(2),'ylim',y2); set(ax(1),'ylim',y1);
        t = title('zoom to desired level, then press enter');
        pause; delete(t)
        t = title('Click peak for IMU2 then for IMU1 (or q twice if done). Left click maxima, right click minima.');
        [x,~,button] = ginput(2);
        if button == 1; ext = 1; else ext = -1; end
        [~,Id] = min(abs(IMU2DN-x(1)));
        [~,Ic] = min(abs(IMU1DN-x(2)));
%         [peaklocC,~] = peakfinder(YC(Ic-10*fs1:Ic+10*fs1),[],[],ext);
%         [peaklocD,~] = peakfinder(YD(Id-10*fs2:Id+10*fs2),[],[],ext);
        [~,peaklocC] = max(YC(Ic-.2*fs1:Ic+.2*fs1)*ext); % max within .2 s. Not much of a window so click well!
        [~,peaklocD] = max(YD(Id-.2*fs2:Id+.2*fs2)*ext);

        x = [IMU2DN(round(peaklocD+Id-.2*fs2-1)) IMU1DN(round(peaklocC+Ic-.2*fs1-1)) ];
        if button ~= 113
            xs(k,:) = x'*24*60*60; k = k+1;
        end
    end
    figure(3); if i == 1; clf; end; try delete(s3(i)); catch; end
    s3(i) = subplot(2,1,i); 
    plot(xs(:,2)/60/60-IMU2DN(1)*24,diff(xs,1,2),'o');
    [p,S] = polyfit(xs(:,2)/60/60-IMU2DN(1)*24,diff(xs,1,2),1);
    set(s3(i),'nextplot','add')
    plot(get(gca,'xlim'),get(gca,'xlim')*p(1)+p(2),'linewidth',2);
    xlabel('Hours since start of IMU2');
    ylabel('Seconds of offset from IMU1era');
    
    figure(i); set(gca,'nextplot','add');
    FSnew = round(40000*(1+p(1)/60/60)); %*fs2*750
    FSold = 40000;
    
    IMU2pnew = resample(IMU2.data.Pressure,FSnew,FSold); % be careful as old resample in cats tools uses timeseries but that is not necessary any more
    IMU2DNnew = ((0:length(IMU2pnew)-1)/fs2/24/60/60+IMU2DN(1) + p(2)/24/60/60)';
    if ~press 
    IMU1p = IMU1.data.Pressure;
    IMU2p = nan(size(IMU2DNnew));
%     IMU2DN = IMU2.data.Date+IMU2.data.Time+TimeOffset_s/24/60/60;
    [~,k1] = min(abs(IMU1DN-IMU2DNnew(1))); [~,k2] = min(abs(IMU1DN-IMU2DNnew(end)));
    [~,kd1] = min(abs(IMU2DNnew-IMU1DN(k1))); [~,kd2] = min(abs(IMU2DNnew-IMU1DN(k2)));
    if size(IMU2p(kd1:IMU2.Hzs.datafs/IMU1.Hzs.datafs:kd2),1)+1 == length(k1:k2)
        k2 = k2-1;
    end
    IMU2p(kd1:IMU2.Hzs.datafs/IMU1.Hzs.datafs:kd2) = IMU1p(k1:k2);
     IMU2p = fixgaps(IMU2p); IMU2p(isnan(IMU2p)) = 0;
      IMU2pnew = IMU2p;%resample(IMU2p,FSnew,FSold);
    end
    
%     newIMU2DN = (((IMU2DN-IMU2DN(1))*24)*p(1)+p(2))/24+IMU2DN(1);
    plot(IMU2DNnew,IMU2pnew,'g--');
    legend('IMU1 p','IMU2 p','New IMU2 p');
end
% disp(num2str(diff(x)*24*60*60));
% disp(num2str(diff(x2)*24*60*60));
%%
IMU2New = IMU2;
IMU2New.data = IMU2.data;
IMU2New.data(length(IMU2DNnew)+1:end,:) = [];
IMU2New.data = [IMU2New.data; IMU2New.data(1:length(IMU2DNnew)-size(IMU2.data,1),:)];
IMU2New.data.Date = floor(IMU2DNnew); IMU2New.data.Time = IMU2DNnew-floor(IMU2DNnew);
for k = 3:size(IMU2.data,2)
IMU2New.data{:,k} = resample(IMU2.data{:,k},FSnew,FSold);
end
if ~press
    IMU2New.data.Pressure = IMU2pnew;
end
IMU2New.Adata = resample(IMU2New.Adata,FSnew,FSold);
IMU2New.TimeOffset_s = p(2);
IMU2New.TimeSlope_s_per_h = p(1);
% IMU2New.TimeOffsetNotes = 'slope was .1372, but changed to 0 since drift evidence was not strong'
ODN = IMU2New.data.Date(1)+IMU2New.data.Time(1);
IMU2New.ODN = ODN;
IMU2New.Atime = ((0:size(IMU2New.Adata,1)-1)/IMU2New.Hzs.accHz/24/60/60+ODN)';
figure; plotyy(IMU2New.data.Date+IMU2New.data.Time,IMU2New.data.Pressure,IMU2New.data.Date+IMU2New.data.Time,[IMU2New.data.A1c1 IMU2New.data.A1c2 IMU2New.data.A1c3])
%%
copyfile([fold ID '\D' tagnum '\' IMU2file],[fold ID '\D' tagnum '\' IMU2file(1:end-4) 'b4resample.mat']);
 save([fold ID '\D' tagnum '\' IMU2file],'-struct','IMU2New');
 saveas(3,[fold ID '\D' tagnum '\timeOffsetCal.jpg'])
  saveas(i,[fold ID '\D' tagnum '\timeOffsetResult.fig'])
  %% Alternate version if just time synch required
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
