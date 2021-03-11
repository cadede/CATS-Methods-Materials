function accwav_audit(M,wsize,pctoverlap)
% ACCWAV_AUDIT - audit and manually detect vocalizations from accelerometry
% William Oestreich
% Goldbogen Lab
% Stanford University
% Building on softward developed by David Cade and James Fahlbusch for
% feeding event (lunge) detection.
%
% INPUTS:
% M - number of minutes to display on interactive scrolling spectrogram
% wsize - window size for FFT
% pctoverlap - percent overlap for sliding window
% 
% OUTPUTS: 
% Saved .mat file containing datenumber, depth, index, and type for
% manually identified vocalizations.
%
% Note:
% Calling this function will prompt the user to select (1) a .wav file of
% the accelerometer data and (2) the PRH .mat file for the same deployment.
% To generate the .wav file from accelerometer data, see "accwav.m".


% select files (acc.wav and prh)
% acc.wav
disp('Select acc.wav file')
[filename,fileloc]=uigetfile('*.*', 'select acc.wav'); 
[y,Fs] = audioread([fileloc filename]);
disp('Loading acc.wav file'); 
% PRH
disp('Select corresponding PRH file')
[filename,fileloc]=uigetfile('*.*', 'select PRH file of the tag deployment');
cd(fileloc);
disp('Loading PRH file'); 
load([fileloc filename(1:end-3) 'mat'],'DN','INFO','tagon','fs','p'); 
whaleName = INFO.whaleName;

% Check to see if we should start at beginning or at a saved index (i)
if ~exist('progressIndex','var')
    i = find(tagon,1);
else
    i = progressIndex;
end

% create arrays for storing call DN, index, and type
C = nan(0,0);
CI = nan(0,0);
Ctype = nan(0,0);
% for syncing indices in prh (usually 10 Hz) and acc (usually 400 Hz) data
ii = i*(Fs/fs);

% scrolling interactive spectrogram for user manual audit
while i<find(tagon,1,'last')
    figure(101); clf; set(gcf,'position',[10,10,1400,700])
    fprintf('Controls:\nLeft click: call (unspecified)\nA: A call (blue whale)\nB: B call  (blue whale)\nC: C call  (blue whale)\nD:D call  (blue whale)\nO: Other\nB: Move Back\nEnter: move forward\nX: move backward\nS: Save');
    
    % end index for display window
    ee=ii+M*60*Fs;
    if ee>length(y); ee = length(y); end
    
    % for spectrogram calculation below
    window=hamming(wsize); 
    noverlap=round(size(window,1)*(pctoverlap/100)); 
    
    % acc x axis
    subplot(311)
    [S,F,T,P] = spectrogram(y(ii:ee,1),window,noverlap,[],Fs,'yaxis');
    Plog=10*log10(P);
    surf(T+(ii/Fs),F,Plog,'edgecolor','none'); axis tight;view(0,90);
    colormap(hot); 
    caxis([-90 max(Plog(:))]);
    ylabel('Frequency (Hz)')
    title('x axis')
    sgtitle(whaleName)
    
    % acc y axis
    subplot(312)
    [S,F,T,P] = spectrogram(y(ii:ee,2),window,noverlap,[],Fs,'yaxis');
    Plog=10*log10(P);
    surf(T+(ii/Fs),F,Plog,'edgecolor','none'); axis tight;view(0,90);
    colormap(hot); 
    caxis([-90 max(Plog(:))]);
    ylabel('Frequency (Hz)')
    title('y axis')
    
    % acc z axis
    subplot(313)
    [S,F,T,P] = spectrogram(y(ii:ee,3),window,noverlap,[],Fs,'yaxis');
    Plog=10*log10(P);
    surf(T+(ii/Fs),F,Plog,'edgecolor','none'); axis tight;view(0,90);
    colormap(hot); 
    caxis([-90 max(Plog(:))]);
    ylabel('Frequency (Hz)')
    title('z axis')
    xlabel('Time into deployment (s)');
       
    % interactive audit buttons 
    button = 1;
    redraw = false;
    while ~isempty(button)
        redraw = false;
            [x,~,button] = ginput(1);
            if isempty(button); continue; end
            switch button
                case 1  %left click = unidentified call
                    dn = DN(1) + x/60/60/24;
                    [~,xI] = min(abs(DN-dn));
                    [C,II] = sort([C;DN(xI)]);
                    CI = sort([CI;xI]);
                    Ctype = [Ctype;1];
                case 97 %a = blue whale A call
                    dn = DN(1) + x/60/60/24;    
                    [~,xI] = min(abs(DN-dn));
                    [C,II] = sort([C;DN(xI)]);
                    CI = sort([CI;xI]);
                    Ctype = [Ctype;2];
                case 98 %b = blue whale B call
                    dn = DN(1) + x/60/60/24;    
                    [~,xI] = min(abs(DN-dn));
                    [C,II] = sort([C;DN(xI)]);
                    CI = sort([CI;xI]);
                    Ctype = [Ctype;3];
                case 99 %c = blue whale C call
                    dn = DN(1) + x/60/60/24;    
                    [~,xI] = min(abs(DN-dn));
                    [C,II] = sort([C;DN(xI)]);
                    CI = sort([CI;xI]);
                    Ctype = [Ctype;4];
                case 100 %d = blue whale D call
                    dn = DN(1) + x/60/60/24;    
                    [~,xI] = min(abs(DN-dn));
                    [C,II] = sort([C;DN(xI)]);
                    CI = sort([CI;xI]);
                    Ctype = [Ctype;5];
                case 111 %o = Other
                    dn = DN(1) + x/60/60/24;    
                    [~,xI] = min(abs(DN-dn));
                    [C,II] = sort([C;DN(xI)]);
                    CI = sort([CI;xI]);
                    Ctype = [Ctype;6];
                case 120 %if x, go backwards
                    i = max(find(tagon,1),i-M*60*fs);
                    redraw = true; button = [];
            end
    end
    if redraw
        continue;   % advance to next chunk of the deployment for audit
    else
        ii = ee;    % update start index for displaying next chunk of the deployment
    end
    
    % save variables in calls.mat file 
    starttime = DN(1);
    prh_fs = fs;
    acc_fs = Fs;
    CallDN = C;
    CallDepth = p(CI);
    CallI = CI;
    CallType = Ctype;
    codes = 'Call types: 1 - unidentified call, 2 - blue whale A call, 3 - blue whale B call, 4 - blue whale C call, 5 - blue whale D call, 6 - other';
    progressIndex = i;
    save([fileloc strtrim(filename(1:end-11)) ' acc_calls.mat'],'codes','CallDN','CallI','CallType','CallDepth','prh_fs','acc_fs','starttime', 'progressIndex');
    
end