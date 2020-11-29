function [data,Adata,Atime,datagaps,ODN,fs,Afs] = truncatedata(data,Adata,Atime,Hzs,fileloc,filename)

% takes data, Adata, Atime, truncates them to get rid of unneeded data
% (i.e. most (or all) of the data before or after deployment)
% looks for any gaps in data or potential bad data sections
% saves truncated file in fileloc with [filename 'truncate.mat'];

DN = data.Date+data.Time;
% DV = datevec(data.Date+data.Time); %makes a datevec of the date and time
fs = round(1./mean((data.Time(50:60)-data.Time(49:59))*24*60*60));
if abs(round(fs)-fs)<.01; fs = round(fs); end
disp(['sample rate: fs = ' num2str(fs) ' Hz']);
Afs = round(1./mean((Atime(50:60)-Atime(49:59))*24*60*60));
if abs(round(Afs)-Afs)<.01; Afs = round(Afs); end
disp(['Acc sample rate: fs = ' num2str(Afs) ' Hz']);

%
ODN = data.Date(1)+data.Time(1);
if sum(data.Pressure) == 0; disp('No pressure sensor'); data.Pressure(:) = 30; p = data.Pressure; nopress = true;
else
    p = data.Pressure; nopress = false;
end


pmin = median(p(1:fs*60)); % min p = median of the first minute;
p = p-pmin; % adjust for an approximate 0 threshold
pmax = prctile(p,97); % 97th percentile of depth (not max in case of eroneous spikes)
p1 = max([1 find(p>.05*pmax,1)-10*60*fs]); % ten minutes before the first dive greater than 10% max
p2 = min(length(p),find(p>.05*pmax,1,'last')+10*60*fs);
%
button = 3;

while ~isempty(button);
    figure (2); clf;
    s1 = subplot(2,1,1);
    plot(1:length(p),p); set(gca,'ydir','rev');
    xlabel('Sample number')
    ylabel('Pressure (raw units)');
    title(s1,{'Can use zoom tool to examine one side or the other. Press enter when ready to bring up a cursor to select boundary of data to maintain'; 'Data outside of yellow boundaries will be deleted'});
    hold on;
    % rectangle('position',[p1 min(p) p2-p1 2*pmax-min(p)],'facecolor','y')
    %  oi = get(s1,'children'); delete(oi(strcmp(get(get(ax1(1),'children'),'type'),'patch')));
    pat(1) = patch([p1 p2 p2 p1],[min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))],'y');
    oi = get(s1,'children'); oi=[oi(2:end); oi(1)];
    set(s1,'children',oi);
    % set(gca,'ylim',[min(p)-5 max(p) + 5]);
    s2 = subplot(2,1,2); plot(1:length(p),[data.Acc1 data.Acc2 data.Acc3]);
    set([s1 s2],'xlim',[1 length(p)])
    set(s2,'xticklabel',datestr(DN(get(s2,'xtick')),'HH:MM:SS'));
    xlabel('Time (labels reset after pressing enter)');
    ylabel('Acc (raw units)');
    pat(2) = patch([p1 p2 p2 p1],[min(get(gca,'ylim')) min(get(gca,'ylim')) max(get(gca,'ylim')) max(get(gca,'ylim'))],'y');
    oi = get(s2,'children'); oi=[oi(2:end); oi(1)];
    set(s2,'children',oi);
    linkaxes([s1 s2],'x');
    z1 = zoom(s1); set(z1,'enable','on'); set(z1,'Motion','horizontal');
    z2 = zoom(s2); set(z2,'enable','on','Motion','both');
    pause;
    set(s2,'xticklabel',datestr(DN(get(s2,'xtick')),'HH:MM:SS'));
    set(gcf,'CurrentAxes',s1)
    title('Now press enter to accept truncated boundaries or align cursor and press 1 to set beginning boundary or 2 to press end boundary');
    [x1,~,button] = ginput(1);
    if button == 49; p1 = round(x1); elseif button == 50; p2 = round(x1); end
    
end

% truncate data
data = data(p1:p2,:);
DN = data.Date+data.Time; [~,aa] = min(abs(DN(1)-Atime)); [~,ba] = min(abs(DN(end)-Atime));
Adata = Adata(aa:ba,:); Atime = Atime(aa:ba);
disp(['New data start time:' datestr(data.Date(1)+data.Time(1),'mm/dd/yyyy HH:MM:SS')]);
disp(['New data end time:' datestr(data.Date(end)+data.Time(end),'mm/dd/yyyy HH:MM:SS')]);

skippeddata = find(diff(DN*24*60*60)>1.5*1/fs); % find spots where gaps are longer than 1.5*1/fs with missing data and replace with nans.  
baddata = nan(size(skippeddata)); % keeps track of length of baddata
for i = length(skippeddata):-1:1
    numnew = length(DN);
    DN = [DN(1:skippeddata(i)); (DN(skippeddata(i))+1/fs/24/60/60:1/fs/24/60/60:DN(skippeddata(i)+1)-.95/fs/24/60/60)'; DN(skippeddata(i)+1:end)];% -.95 because of occasional rounding errors
    numnew = length(DN)-numnew;
    oi = array2table(nan(numnew,size(data,2)));
    oi.Properties.VariableNames = data.Properties.VariableNames;
    data = [data(1:skippeddata(i),:); oi; data(skippeddata(i)+1:end,:)];
    baddata(i) = numnew;
    skippeddata(i+1:end) = skippeddata(i+1:end) + numnew;
end
skippeddataI = skippeddata;
Length = baddata;
datagaps = table(skippeddataI,Length);
disp([num2str(length(skippeddata)) ' data gaps replaced with nans.  Total # of missed points = ' num2str(sum(isnan(data.Pressure)))]);
data.Date = floor(DN);
data.Time = DN-floor(DN);

if any(diff(DN)<0); 
    disp(datestr(DN(diff(DN)<0),'mmm-dd HH:MM:SS.fff'));
    error('Bad time stamps! Check data at above times'); 
end

skippeddata = find(diff(Atime*24*60*60)>1.5*1/Afs); % find spots with missing data in full Adata.  
for i = length(skippeddata):-1:1
    numnew = length(Atime);
    Atime = [Atime(1:skippeddata(i)); (Atime(skippeddata(i))+1/Afs/24/60/60:1/Afs/24/60/60:Atime(skippeddata(i)+1)-.95/Afs/24/60/60)'; Atime(skippeddata(i)+1:end)];% -.95 because of occasional rounding errors
    numnew = length(Atime)-numnew;
    oi = nan(numnew,size(Adata,2));
    Adata = [Adata(1:skippeddata(i),:); oi; Adata(skippeddata(i)+1:end,:)];
end

if ~isempty(baddata)
    f = input('Fix gaps in bad data (may be bad idea if more than a few seconds of data)? 1 = yes, 2 = no   ');
    if f == 1
        disp('Fixing bad data');
        headers = data.Properties.VariableNames;
        fixcols = true(size(headers));
        fixcols(~cellfun(@isempty, cellfun(@(x) strfind(x,'GPS'),headers,'uniformoutput',false))) = false;
        fixcols(cellfun(@(x) strcmp(x,'Date'),headers) | cellfun(@(x) strcmp(x,'Time'),headers)) = false;
        for j = find(fixcols)
            if sum(isnan(data{j,:})) ~= length(data{j,:})
                data{:,j} = fixgaps(data{:,j});
            end
        end
    else
        disp('Nans left in data, recommend checking raw and imported data to identify the problem');
    end
end

newdatafile = [filename(1:end-4) 'truncate.mat'];

lastwarn('');
try
    save([fileloc newdatafile],'data','Adata','Atime','ODN','Hzs','datagaps');
    if ~isempty(lastwarn)
        error(lastwarn);
    end
catch %v7.3 allows for bigger files, but makes a freaking huge file if used when you don't need it
    save([fileloc newdatafile],'data','Adata','Atime','ODN','Hzs','datagaps','-v7.3');
    disp('Made a version 7.3 file in order to include all data (files were large)');
end

