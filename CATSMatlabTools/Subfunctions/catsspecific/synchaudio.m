function audstart = synchaudio(tagon,DN,fs,fileloc,p,Atemp)

timestamp = input('Use tag on animal time as synch point? (1 = yes, 2 = no) ');
if timestamp == 2; timestamp = input('Enter time of synch point in date vector format ([yyyy mm dd HH MM SS]): ');
    timestamp = datenum(timestamp);
else timestamp = find(tagon,1); timestamp = DN(timestamp);
end
try cd([fileloc 'AudioData\']); catch; cd(fileloc); end
[~,I] = min(abs(DN-timestamp)); I = I-fs*10:I+fs*10;
[audfile,audloc] = uigetfile('*.wav','Get first audio file from "AudioData" folder (assumes all audiofiles from here are consecutive with no gaps)');
if exist([audloc audfile(1:end-4) 'audio.mat'],'file'); load([audloc audfile(1:end-4) 'audio.mat']); else
    aud = struct();
    try [aud.data,aud.rate,aud.bits] = wavread([audloc audfile]);
    catch
        audioI = audioinfo([audloc audfile]);
        [aud.data,aud.rate] = audioread([audloc audfile]);
        aud.bits = audioI.BitsPerSample;
    end
end
FF = figure(100); clf;
try
    II = strfind(audfile,'-');
    otime = datenum(audfile(II(end-2)+1:II(end)+3),'yyyymmdd-HHMMSS-fff');
    starttimedif = abs(otime-DN(1));
catch
    warning('starttime of audiofile could not be read from file name, using data start time');
    starttimedif = 0; otime = DN(1);
end
if starttimedif>.001/24/60/60
    warning(['Starttime of data is ' datestr(DN(1),'yyyy-mmm-dd HH:MM:SS.fff') ', but starttime in audio name is ' datestr(otime,'yyyy-mmm-dd HH:MM:SS.fff')]);
    xx = input('Okay to proceed using audio start time from file name? (1 = yes, 2 = no):');
    if xx~=1; error('check audio start time'); end
end


audioI = round((timestamp-otime)*24*60*60*aud.rate); audioI = audioI-aud.rate*10:audioI+aud.rate*10;
sp1 = subplot(2,1,1); hold on;
ax = plotyy(DN(I),p(I),DN(I),Atemp(I,3)); legend('Depth','Az')
set(ax(1),'ydir','rev')
set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
set(ax,'xlim',[DN(I(1)) DN(I(end))]);
t = title('Click on synch location in upper graph (will select peak within 1s)');
set(ax,'nextplot','add')
subplot(2,1,2);
nfft = 2048; 
nov = nfft/2;
w = hanning(nfft) ;
spec = nan(200,1); k = 1; % 1/10 of a second bins
for i = audioI(1):aud.rate/10:audioI(end)
    x = aud.data(i:i+round(aud.rate/10),1);
    [SL,f] = spectrum_level(x,nfft,aud.rate,w,nov);

    [TOL,fc]=spec2tol(SL,f);
    spec(k,1:length(fc)) = TOL';
    k = k+1;
end
pI = (1:size(spec,1))/10;

% pI = (1:length(audioI))/aud.rate;
plot(pI,spec);%aud.data(audioI));
hold on; 
plot(pI,mean(spec,2),'k','linewidth',3)
lgd = legend(num2str(round(fc')));
lgd.Title.String = 'Center Frequency';
xlabel('Seconds')
xlim([pI(1) pI(end)])
[x,~] = ginput(1);
[~,x] = min(abs(DN-x));
[i1,mag] = peakfinder(abs(Atemp(x-fs:x+fs,3)));
[~,i] = max(mag); i1 = i1(i);
x = x-fs-1+i1;
plot(ax(2),DN(x),Atemp(x-fs+i1,3),'rs','markerfacecolor','r')
set(ax,'xlim',[DN(x-10*fs) DN(x+10*fs)])
set(t,'string','');
subplot(2,1,2); hold on;
title('Click on synch location in lower graph (will select peak within 1s). Look for where noise is broad band and high');
[x2,~] = ginput(1);
% x2 = x2*aud.rate;
% x2 = audioI(1)-1+round(x2);
% [i1,mag] = peakfinder(aud.data(x2-aud.rate:x2+aud.rate));
% [~,i] = max(mag); i1 = i1(i);
% x2 = x2-aud.rate-1+i1;
% plot(pI(x2-audioI(1)+1),mag(i),'rs','markerfacecolor','r')
oi = mean(spec,2);
[i1,mag] = peakfinder(oi(x2*10-10:x2*10+10));
[~,i] = max(mag); i1 = round(i1(i)+x2*10-10-1);
plot(pI(i1),mag(i),'rs','markerfacecolor','r')
x2 = pI(i1)*aud.rate+audioI(1);
audioI2 = x2-audioI(1)+1-10*aud.rate:x2-audioI(1)+1+10*aud.rate;
xlim([audioI2(1) audioI2(end)]/aud.rate)
audstart = DN(x)-x2/aud.rate/24/60/60;
disp(['Old Audio start time: ' datestr(DN(1), 'dd-mmm-yyyy HH:MM:SS.fff')]);
disp(['New Audio Start time: ' datestr(audstart,'dd-mmm-yyyy HH:MM:SS.fff')]);% '. Will be asked to rename audio files after prh process is complete.']);
kk = input('Rename audio files with new time adjument (recommended)? 1 = yes, 2 = no : ');
if kk ~=2
    D = dir(audloc); D = {D.name}; D = D(3:end);
    for kk = 1:length(D)
        II = strfind(D{kk},'-');
        otime = datenum(D{kk}(II(end-2)+1:II(end)+3),'yyyymmdd-HHMMSS-fff');
        newtime = otime+audstart-DN(1);
        newfile = [D{kk}(1:II(end-2)) datestr(newtime,'yyyymmdd-HHMMSS-fff') D{kk}(II(end)+4:end)];
        disp(['old file = ' D{kk}]);
        disp(['new file = ' newfile]);
        movefile([audloc D{kk}],[audloc newfile]);
    end
end
%     close(FF);