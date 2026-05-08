function tvidDN = autosynchvideousingaudio(vidNum,tagon,DN,vidDN,vidDurs,fs,p,A,audfold,movaudfold)
cf = pwd;
cd (movaudfold);
movaudfiles = dir('*.wav');
movaudfiles = {movaudfiles.name};
for i = 1:length(movaudfiles)
    if str2num(movaudfiles{i}(end-7:end-4)) == vidNum
        movaudfile = movaudfiles{i};
        break
    end
end
disp(['Synching audio file from movie: ' movaudfile(1:end-4)])

cd (audfold);
audfiles = dir('*.wav');
audfiles = {audfiles.name};
audtimes = nan(length(audfiles),1);
for i = 1:length(audfiles)
    try
    audfile = audfiles{i};
     II = strfind(audfile,'-');
    audtimes(i) = datenum(audfile(II(end-2)+1:II(end)+3),'yyyymmdd-HHMMSS-fff');
    catch
        error(['audio start time could not be read from audio file ' audfile]);
    end
end
t1 = find(tagon,1); t2 = find(tagon,1,'last');
fileA = find(audtimes<max(vidDN,DN(t1)),1,'last');
fileB = find(audtimes<vidDN+vidDurs/24/60/60,1,'last');
disp('Using audio files:');
for i = fileA:fileB
    disp(audfiles{i});
end

cd (cf);

if exist([movaudfold movaudfile(1:end-4) 'audio.mat'],'file'); load([movaudfold movaudfile(1:end-4) 'audio.mat']);
else
    movaud = struct();
%     try [movaud.data,movaud.rate,movaud.bits] = wavread([movaudfold movaudfile]);
%     catch
        audioI = audioinfo([audfold audfile]);
        [movaud.data,movaud.rate] = audioread([movaudfold movaudfile]);
        movaud.bits = audioI.BitsPerSample;
%     end
end
movaud.dur = size(movaud.data,1)/movaud.rate;

aud2 = struct(); aud2.data = nan(0,1);
audstart = max(audtimes(fileA),vidDN);
movaudstart = vidDN;
for i = fileA:fileB % load all audio files that match video file and make a single file
    audfile = audfiles{i};
%     if exist([audfold audfile(1:end-4) 'audio.mat'],'file'); load([audfold audfile(1:end-4) 'audio.mat']);
%     else
        aud = struct();
        %         try [aud.data,aud.rate,aud.bits] = wavread([audfold audfile]);
        %         catch
        audioI = audioinfo([audfold audfile]);
        a = round(max((vidDN-audtimes(i))*24*60*60*audioI.SampleRate+1,1));
        b = round(min(audioI.TotalSamples,(vidDN+vidDurs/24/60/60-audtimes(i))*24*60*60*audioI.SampleRate));
        [aud.data,aud.rate] = audioread([audfold audfile],[a b]);
        aud.bits = audioI.BitsPerSample;
        %         end
%     end
    if i>fileA && (aud2.rate~=aud.rate || aud2.bits~=aud.bits); error('Difference in sample rate between audio files'); end
    aud2.data = [aud2.data; aud.data(:,1)]; aud2.rate = aud.rate; aud2.bits = aud.bits;
    aud2.dur = size(aud2.data,1)/aud2.rate;
end
aud = aud2;
clear aud2;
aud.start = audstart;
movaud.start = movaudstart;
getauda = @(aud,dn) max(1,round((dn-aud.start)*24*60*60*aud.rate));
getaudb = @(aud,dn) min(size(aud.data,1),round((dn-aud.start)*24*60*60*aud.rate+aud.dur*aud.rate));
%%
[pI,pmag] = peakfinder(-p(tagon),2,-6);
pI = pI+t1-1; pmag = -pmag;

figure(vidNum);
set(vidNum,'units','normalized','outerposition',[0 0 1 1]);  clf;
sp = nan(3,1);
sp(1) = subplot(6,1,1); hold on;
[~,pa] = min(abs(DN-movaudstart));
pa = max(pa,t1);
[~,pb] = min(abs(DN-(movaudstart+movaud.dur/24/60/60)));
pb = min(pb,t2);
I = pa:pb;
ax = plotyy(DN(I),p(I),DN(I),A(I,:)); 
set(ax(1),'ydir','rev')
ylabel('Depth (m)','parent',ax(1))
set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
set(ax,'xlim',[DN(I(1)) DN(I(end))]);
t = title({['Movie: ' movaudfile(1:end-4)]; 'Examine data, the press enter twice to accept offset or press enter then "s" to manually synch hydrophone and camera mic for this video'});
set(ax,'nextplot','add')
plot(ax(1),DN(pI),pmag,'rs','markerfacecolor','r')
legend('Depth','Surfacing','Ax','Ay','Az');%,'Location','northeastoutside')

a = getauda(aud,DN(pa));
b = getaudb(aud,DN(pb));
ma = getauda(movaud,DN(pa));
mb = getaudb(movaud,DN(pb));

audioI = a:b;
nfft = 512; 
nov = nfft/2;
w = hanning(nfft) ;
spec = nan(200,1); k = 1; % 1/10 of a second bins
for i = audioI(1):round(aud.rate/25):audioI(end)-aud.rate/25
    x = aud.data(i:i+round(aud.rate/25),1);
    [SL,f] = spectrum_level(x,nfft,aud.rate,w,nov);

    [TOL,fc]=spec2tol(SL,f);
    spec(k,1:length(fc)) = TOL';
    k = k+1;
end

% pI = (1:size(spec,1))/10;
aI = round(25*((1:size(spec,1))/25+(a-1)/aud.rate))/25/24/60/60+aud.start;
% pI = (1:length(audioI))/aud.rate;
sp(2) = subplot(6,1,2:3);
ax(3) = gca;
h = plot(aI,spec);%aud.data(audioI));
hold on; 
audmean = mean(spec,2);
h(end+1) = plot(aI,audmean,'k','linewidth',3);
lgd = legend(num2str(round(fc')));%,'location','eastoutside');
lgd.Title.String = {'Center Frequency'; 'of 1/3 octave band'};
% xlabel('Seconds')
set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
xlim([aI(1) aI(end)])
title('Audio from hydrophone (should be synched with data)')
% audioI = round((movaudstart-otime)*24*60*60*aud.rate); 
audioI = ma:mb;
% sp(3) = subplot(3,1,3);

nfft = 512; 
nov = nfft/2;
w = hanning(nfft) ;
spec = nan(200,1); k = 1; % 1/10 of a second bins
% tic
for i = audioI(1):round(movaud.rate/25):audioI(end)-movaud.rate/25
    x = movaud.data(i:i+round(movaud.rate/25),1);
    [SL,f] = spectrum_level(x,nfft,movaud.rate,w,nov);

    [TOL,fc]=spec2tol(SL,f);
    spec(k,1:length(fc)) = TOL';
    k = k+1;
end
% toc
maI = round(25*((1:size(spec,1))/25+(ma-1)/movaud.rate))/25/24/60/60+movaud.start;
sp(3) = subplot(6,1,4:5);
ax(4) = gca;
% pI = (1:length(audioI))/aud.rate;
plot(maI,spec);%aud.data(audioI));
hold on; 
movaudmean = mean(spec,2);
h1 = plot(maI,movaudmean,'k','linewidth',3);
% lgd = legend(num2str(round(fc')),'location','eastoutside');
% lgd.Title.String = {'Center Frequency'; 'of 1/3 octave band'};
% xlabel('Seconds')
set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
xlim([aI(1) aI(end)])
title ('Audio from movie files')

linkaxes(ax(1:4),'x'); %set(ax(4),'xlim',get(ax(3),'xlim'))
% [ax(3:4),h1,h2] = plotyy(aI,amp,maI,mamp);
% set(h1,'color','k','linewidth',3)
% set(h2,'color','r','linewidth',3)
[c,lags] = xcorr(audmean,movaudmean,250);% 10 seconds on each side;
[~,ci] = max(c);
offset = lags(ci)/25;
set(ax(4),'nextplot','add')
h2 = plot(ax(4),maI+offset/24/60/60,movaudmean,'b--','linewidth',3);
% title('Amplitude of acoustic data, 1/25 s bins');
% [~,amag] = peakfinder(amp);
% set(ax(3),'ylim',[0 1.1*percentile(amag,.9)])
% [~,mmag] = peakfinder(mamp,percentile(mamp,.5));
% set(ax(3),'ylim',[0 1.1*percentile(mmag,.9)])
legend([h1 h2],'camera mic','Camera mic auto aligned ')
subplot(6,1,6);
plot(lags/25,c,lags(ci)/25,c(ci),'ro'); 
xlabel('seconds');
title(['Cross Correlation coefficient at each offest, chosen value = ' num2str(offset) ' s']);

tvidDN = vidDN+offset/24/60/60;
% PI = pI;
%%
% [pI,pmag] = peakfinder(-p(tagon),2,-6);
% pI = pI+t1-1; pmag = -pmag;
% 
% figure(vidNum);
% set(vidNum,'units','normalized','outerposition',[0 0 1 1]);  clf;
% sp = nan(3,1);
% sp(1) = subplot(3,1,1); hold on;
% [~,pa] = min(abs(DN-movaudstart));
% pa = max(pa,t1);
% [~,pb] = min(abs(DN-(movaudstart+movaud.dur/24/60/60)));
% pb = min(pb,t2);
% I = pa:pb;
% ax = plotyy(DN(I),p(I),DN(I),A(I,:)); 
% set(ax(1),'ydir','rev')
% ylabel('Depth (m)','parent',ax(1))
% set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
% set(ax,'xlim',[DN(I(1)) DN(I(end))]);
% t = title({['Movie: ' movaudfile(1:end-4)]; 'Examine data, the press enter to accept offset or press "s" to manually synch hydrophone and camera mic for this video'});
% set(ax,'nextplot','add')
% plot(ax(1),DN(pI),pmag,'rs','markerfacecolor','r')
% legend('Depth','Surfacing','Ax','Ay','Az')
% 
% a = getauda(aud,DN(pa));
% b = getaudb(aud,DN(pb));
% ma = getauda(movaud,DN(pa));
% mb = getaudb(movaud,DN(pb));
% 
% audioI = a:b;
% nfft = 2048; 
% nov = nfft/2;
% w = hanning(nfft) ;
% spec = nan(200,1);
% k = 1; % 1/25 of a second bins
% % amp = nan(0,1);
% % for i = audioI(1)+ceil(aud.rate/4):round(aud.rate/25):audioI(end)-ceil(aud.rate/4)
% %     x = aud.data(i-round(aud.rate/4):i+round(aud.rate/4),1);%audioI(1):round(aud.rate/25):audioI(end)-aud.rate/25
% %      amp(k,1) = sqrt(mean(x.^2));
% %     k = k+1;
% % end
% for i = audioI(1):round(aud.rate/10):audioI(end)-aud.rate/10
%     x = aud.data(i:i+round(aud.rate/10),1);
%     [SL,f] = spectrum_level(x,nfft,aud.rate,w,nov);
% 
%     [TOL,fc]=spec2tol(SL,f);
%     spec(k,1:length(fc)) = TOL';
%     k = k+1;
% end
% % pI = (1:size(spec,1))/10;
% aI = round(25*((1:size(amp,1))/25+(a-1)/aud.rate))/25/24/60/60+aud.start;
% % pI = (1:length(audioI))/aud.rate;
% % h = plot(aI,spec);%aud.data(audioI));
% 
% 
% % audioI = round((movaudstart-otime)*24*60*60*aud.rate); 
% audioI = ma:mb;
% sp(3) = subplot(3,1,3);
% ax(4) = gca;
% % nfft = 2048; 
% % nov = nfft/2;
% % w = hanning(nfft) ;
% % spec = nan(200,1);
% k = 1; % 1/2 second bins at 1/25 s resolution
% mamp = nan(0,1);
% for i = audioI(1)+ceil(movaud.rate/4):round(movaud.rate/25):audioI(end)-ceil(movaud.rate/4)
%     x = movaud.data(i-round(movaud.rate/4):i+round(movaud.rate/4),1);
% %     [SL,f] = spectrum_level(x,nfft,movaud.rate,w,nov);
% % 
% %     [TOL,fc]=spec2tol(SL,f);
% %     spec(k,1:length(fc)) = TOL';
%     mamp(k,1) = sqrt(mean(x.^2));
%     k = k+1;
% end
% amp =10*log10(amp);
% mamp = 10*log10(mamp);
% sp(2) = subplot(3,1,2);
% ax(3) = gca;
% hold on; 
% % audmean = mean(spec,2);
% % h(1) = plot(aI,amp,'k','linewidth',3);
% % lgd = legend(num2str(round(fc')),'location','eastoutside');
% % lgd.Title.String = {'Center Frequency'; 'of 1/3 octave band'};
% % xlabel('Seconds')
% set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
% xlim([aI(1) aI(end)])
% % title('Audio from hydrophone (should be synched with data)')
% 
% maI = round(25*((1:size(mamp,1))/25+(ma-1)/movaud.rate))/25/24/60/60+movaud.start;
% 
% % pI = (1:length(audioI))/aud.rate;
% %plot(maI,spec);%aud.data(audioI));
% hold on; 
% % movaudmean = mean(spec,2);
% [ax(3:4),h1,h2] = plotyy(aI,amp,maI,mamp);
% set(h1,'color','k','linewidth',3)
% set(h2,'color','r','linewidth',3)
% [c,lags] = xcorr(amp,mamp,250);% 10 seconds on each side;
% [~,ci] = max(c);
% offset = lags(ci)/25;
% set(ax(4),'nextplot','add')
% plot(ax(4),maI+offset/24/60/60,mamp,'b--','linewidth',2);
% title('Amplitude of acoustic data, 1/25 s bins');
% % [~,amag] = peakfinder(amp);
% % set(ax(3),'ylim',[0 1.1*percentile(amag,.9)])
% % [~,mmag] = peakfinder(mamp,percentile(mamp,.5));
% % set(ax(3),'ylim',[0 1.1*percentile(mmag,.9)])
% legend('hydrophone', 'camera mic','Camera mic auto aligned ')
% % lgd = legend(num2str(round(fc')),'location','eastoutside');
% % lgd.Title.String = {'Center Frequency'; 'of 1/3 octave band'};
% % xlabel('Seconds')
% set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
% xlim([aI(1) aI(end)])
% % title ('Audio from movie files')
% linkaxes(ax(1:4),'x'); %set(ax(4),'xlim',get(ax(3),'xlim'))
% 
% subplot(3,1,3);
% plot(lags/25,c,lags(ci)/25,c(ci),'ro'); 
% xlabel('seconds');
% title(['Cross Correlation coefficient at each offest, chosen value = ' num2str(offset) ' s']);
% 
% PI = pI;
% pause
% %%
% button = [];
% k = 1;
% omovstart = movaud.start;
% adj = nan(0,2); m = adj;
% pI = find(DN(PI)>movaudstart & DN(PI)<movaudstart+movaud.dur/24/60/60);
% pI = PI(pI);
% curadj = 0;
% if isempty(pI)
%     set(t,'string',{'No surfacings detected! Use controls to zoom to time window with a suitable synch point'; 'Then press enter'});
%     pause;
%     [~,pI] = min(abs(DN-diff(get(ax(1),'xlim'))));
% end
% 
% while (isempty(button) || button ~= 113) && k<= length(pI)
%     set(ax(1:3),'xlim',[DN(pI(k)-10*fs) DN(pI(k)+10*fs)])
%     set(ax(4),'xlim',get(ax(3),'xlim')+curadj)
%     set(ax(1:3),'xticklabel',datestr(get(ax(1),'xtick'),'HH:MM:SS'));
%      set(ax(4),'xticklabel',datestr(get(ax(4),'xtick'),'HH:MM:SS'));
%     set(t,'string',{'Left click peak in (middle) audio graph, then peak in (bottom) movie graph. '; 'Press "enter" to move to next surfacing, "b" to move back one surfacing, "u" to undo last synch point.'; 'Press q when synch is satisfactory for this movie.'})
%     [x,~,button]=ginput(1);
%     if isempty(button)
%         k = k+1; 
%         if k>length(pI)
%             if curadj == 0; warning('No time adjustment made!'); end
%             tvidDN = vidDN-curadj;
%         end
%         continue;
%     elseif button == 98 % b
%         k = max(k-1,1); continue;
%     elseif button == 113 % q
%         if curadj == 0; warning('No time adjustment made!'); end
%         tvidDN = vidDN-curadj;
%         break;
%     elseif button == 1 % left click
%         adj(end+1,1) = x;
%         a = (x-aI(1))*24*60*60; % a is seconds since start time
%         [i1,mag] = peakfinder(audmean(a*10-5:a*10+5)); % audmean is in 1/10 s bins
%         [~,i] = max(mag); i1 = round(i1(i)+a*10-5-1);
%         m(size(adj,1),1) = plot(ax(3),aI(i1),mag(i),'rs','markerfacecolor','r');
%         [x,~,button]=ginput(1);
%         adj(end,2) = x;
%         a = (x-maI(1))*24*60*60; % a is seconds since start time
%         [i1,mag] = peakfinder(movaudmean(a*10-5:a*10+5)); % audmean is in 1/10 s bins
%         [~,i] = max(mag); i1 = round(i1(i)+a*10-5-1);
%         m(size(adj,1),2) = plot(ax(4),maI(i1),mag(i),'rs','markerfacecolor','r');
%         curadj = mean(diff(adj,[],2));
%         disp('Current differences (s):')
%         disp(diff(adj,[],2)*24*60*60)
%         disp('Mean difference (s):')
%         disp(curadj*24*60*60)
%         set(ax(4),'xlim',get(ax(3),'xlim')+curadj)
%         set(ax(4),'xticklabel',datestr(get(ax(4),'xtick'),'HH:MM:SS'));
%     elseif button == 117 % u
%         try adj(end,:) = []; delete(m(size(adj,1)+1,1)); delete(m(size(adj,1)+1,2)); 
%         catch; warning('No synchpoints left to undo'); 
%         end
%     end
% end
% clf;

