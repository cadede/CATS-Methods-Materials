function [report, data] = addGPSfrompos(fileloc,data,Hzs,UTC) 
% finds all pos files in a filelocation, and adds GPS info to the matching mat files and the prh file (if they exis.  Also generates a "quick look" with depth and GPS locations
global fixweird
fixweird = [];
%%
if nargin<1
    [~,fileloc] = uigetfile('*.*','Choose any file in the directory of your .pos file (will examine all subdirectories as well');
end
%%
if strcmp(fileloc(end),'\'); fileloc = fileloc(1:end-1); end

% clearvars -except fileloc data Hzs
[D,F] = subdir(fileloc);
D = [{fileloc} D];
oi = dir(fileloc); oi2 = {oi.name};
F = [{oi2(~vertcat(oi.isdir))} F];
posmatch = cellfun(@(x) strfind(x,'.pos'),F,'uniformoutput',false);
posstat = cellfun(@(x) strfind(x,'.pos.stat'),F,'uniformoutput',false);
posfastloc = cellfun(@(x) strfind(x,'Fastloc'),D,'uniformoutput',false);

for j = 1:length(posmatch); posmatch{j} = ~cellfun(@isempty,posmatch{j})&cellfun(@isempty,posstat{j}); end
pos = cell(size(posmatch));
Dmatch = false(size(posmatch));
for j = 1:length(posmatch)
    Dmatch(j) = sum(posmatch{j});
    if ~isempty(posfastloc{j}); Dmatch(j) = false; for jj = 1:length(posmatch{j}); posmatch{j}(jj) = false; end; continue; end
    for jj = 1:sum(posmatch{j}) %shouldn't be more than 1...
        if jj == 2; error(['more than one pos file in ' D{j}]); end
        posfile = F{j}(posmatch{j}); posfile = posfile{jj};
        done = false;
        f = fopen([D{j} '\' posfile]); numheads = 0; skiprest = false;
        while ~done
            try tline = fgets(f); catch; disp([posfile ' has no data']); Dmatch(j) = 0; posmatch{j}(jj) = 0; skiprest = true; break; end
            if strcmp(tline(1),'%')||strcmp(tline(2),'%'); numheads = numheads+1; else done = true; end; %disp(numheads); disp(done);
        end
        if skiprest; continue; end
        fclose(f);
        posT = importdata([D{j} '\' posfile],' ',numheads);
        posT.textdata(1:size(posT.textdata,1)-size(posT.data,1),:) = [];
        DN = datenum(posT.textdata(:,1),'yyyy/mm/dd') + datenum(posT.textdata(:,2),'HH:MM:SS.fff')-floor(datenum(posT.textdata(:,2),'HH:MM:SS.fff')); %datenum in GMT
        pos{j} = [DN posT.data];
    end
end
if nargin == 2
    Dmatch = Dmatch(1);
    DATA= {data}; clear data;
else DATA = cell(0);
end
report = nan(0,2);
%%
for j = 1:length(Dmatch)
        if j~=1; DATA = cell(0); end
    if Dmatch(j)
% j = 1;
        %add GPS to data files
        posfile = F{j}{posmatch{j}}; %posfile = posfile{jj};
        data1 = cellfun(@(x) strfind(x,[posfile(1:end-4) '.mat']),F{j},'uniformoutput',false);
        data2 = cellfun(@(x) strfind(x,[posfile(1:end-4) 'truncate.mat']),F{j},'uniformoutput',false);
        notdata = cellfun(@(x) strfind(x,'ALLDATA'),F{j},'uniformoutput',false);
        datafile = F{j}((~cellfun(@isempty,data1)|~cellfun(@isempty,data2))&cellfun(@isempty,notdata));
        if ~isempty(datafile)
            for k = 1:length(datafile); data1 = load([D{j} '\' datafile{k}],'data'); DATA{end+1} = data1.data; end
            [~,data1] = max(cellfun(@(x) size(x,1), DATA)); % which one is longer
            %             oi =find(~isnan(DATA{data1}.GPSTime),1); %first GPStime recorded
            %             GPSoffset = round((DATA{data1}.GPSTime(oi)+DATA{data1}.GPSDate(oi)-pos{j}(1))*24*2)/2; % GPSoffset to nearest half hour from when data started, most tags this will be 0, but some older had GPStime in local and pos file in UTC
%             tagnum = gettagnum(posfile);
            % get UTC offset (only needed in special cases and for the kml file
            lat1 = nanmedian(pos{j}(:,2)); long1 = nanmedian(pos{j}(:,3));
%             if ~exist('cityGPS','var')
%                 [UTC,cityGPS,citytext,UTCs,TZs] = getUTC(lat1,long1,pos{j}(1));
%             else UTC = getUTC(lat1,long1,pos{j}(1),cityGPS,citytext,UTCs,TZs);
%             end
%             if isempty(tagnum) || (tagnum~=6&&tagnum~=7&&tagnum~=36&&tagnum~=37)
%                 GPSoffset = 0;
%                 if isempty(tagnum); disp('Couldn''t get tag num from file name, but okay as long as it''s not 6 or 7'); end
%             else 
%                 GPSoffset = UTC;
%             end
            % test offset
           GPSoffset = 0;
            Tdata = DATA{data1};
            GPSDN = Tdata.GPSDate+Tdata.GPSTime-GPSoffset/24;
            try load([D{j} '\' F{j}{~cellfun(@isempty,notdata)}],'timedif'); 
                catch; try load([fileloc '\' F{j}{~cellfun(@isempty,cellfun(@(x) regexp(x,'prh'),F{j},'UniformOutput',false))}],'INFO'); timedif = INFO.timedif;
                catch; timedif = input('timedif?  (from header xls file)'); end; end
            DNt = Tdata.Date+Tdata.Time-GPSoffset/24 + timedif/24;  %not sure if timedif is needed?  
            [GPSDN,fixjump] = fixGPSdate(DNt,GPSDN); %looks for overly repeated times (for some tags when GPS is underwater)
%             DNt = DNt + nanmean(GPSDN-DNt);
            Tp = Tdata.Pressure;
            if sum(Tp) == 0; Tp = load([fileloc '\' F{j}{~cellfun(@isempty,cellfun(@(x) regexp(x,'prh'),F{j},'UniformOutput',false))}],'p','DN','INFO','fs'); Tp.timedif = Tp.INFO.timedif;
                Tp.DN = Tp.DN-Tp.timedif/24; [~,a1] = min(abs((Tdata.Date+Tdata.Time)-Tp.DN(1))); 
                dfs = round(1./mean((Tdata.Time(50:60)-Tdata.Time(49:59))*24*60*60)); a2 = a1+dfs/Tp.fs*(length(Tp.p))-1;
                tTp = nan(size(Tdata.Time)); tTp(a1:dfs/Tp.fs:a2) = Tp.p; tTp = fixgaps(tTp); tTp([1:a1-1 a2+1:end]) = 0;
                Tp = tTp; clear tTp;
            end
                
            figure(641); clf; plot(GPSDN,Tp,pos{j}(:,1),pos{j}(:,2)./pos{j}(:,2)*min(Tp),'rs')
            set(gca,'ydir','rev');
            title('red squares are the timestamps of gps points');
            set(641,'units','pixels');
            set(641,'Position', get(0,'Screensize')); % Maximize figure. 
%             xs = get(gca,'xlim'); ys = get(gca,'ylim'); %tt = text(xs(1),ys(2),'Left click just after the last GPS hit before deployment, Right click to zoom','fontsize',14,'verticalalignment','top');
            button = [];
            status = 1;
%             [x,~,button] = ginput(1);
%             while ~isempty(button)
%                 switch button
%                     case 3
%                         set(gca,'xlim',[max(xs(1),x-diff(xs)/4),min(xs(2),x+diff(xs)/4)]); xs = get(gca,'xlim'); set(gca,'ylim',[min(Tp) max(Tp(GPSDN>xs(1)&GPSDN<xs(2)))]); ys = get(gca,'ylim');
%                         posit = get(tt,'position'); posit(1) = xs(1); posit(2) = ys(2); set(tt,'position',posit);
%                         [x,~,button] = ginput(1); continue;
%                     case 1
%                         if status == 1; delete(tt); tt = text(xs(1),ys(2),'Now, Left click just after the last depth before the first dive, Right click to zoom','fontsize',14,'verticalalignment','top');
%                             GPSi = find(pos{j}(:,1)<x,1,'last'); [x,~,button] = ginput(1); status = 2; continue;
%                         elseif status == 2;
%                             delete(tt); 
%                             px = find(GPSDN<x,1,'last');
%                             newOffset = (pos{j}(GPSi,1)-GPSDN(px))*24;
%                             GPSoffset = GPSoffset - newOffset;
%                             GPSDN = Tdata.GPSDate+Tdata.GPSTime-GPSoffset/24;             
%                             GPSDN = fixGPSdate(DNt,GPSDN,[],[],fixjump);
%                             clf; plot(GPSDN,Tp,pos{j}(:,1),pos{j}(:,2)./pos{j}(:,2)*min(Tp),'rs');  xs = get(gca,'xlim'); ys = get(gca,'ylim');
%                             tt = text(xs(1),ys(2),'Press Enter if Done, else press ''R'' to redo','fontsize',14,'verticalalignment','top');
%                             [x,~,button] = ginput(1); continue;
%                         end
%                     case 114
%                         status = 1;
%                         delete(tt); tt = text(xs(1),ys(2),'Left click just after the last GPS hit before deployment, Right click to zoom','fontsize',14,'verticalalignment','top');
%                         [x,~,button] = ginput(1); continue;
%                     otherwise
%                         [x,~,button] = ginput(1); continue;
%                 end
%             end
       end
        %

        for jj = 1:length(DATA)
            GPSDN = DATA{jj}.GPSDate+DATA{jj}.GPSTime-GPSoffset/24;
            GPSDN = fixGPSdate(DATA{jj}.Date+DATA{jj}.Time,GPSDN,[],[],fixjump);
            DATA{jj}.Lat = nan(size(DATA{jj}.Time));
            DATA{jj}.Long = nan(size(DATA{jj}.Time));
            DATA{jj}.GPSsd = nan(size(DATA{jj}.Time,1),3); % the sdne vector to estimate error
            DATA{jj}.Elev = nan(size(DATA{jj}.Time));
            a = find(diff(GPSDN)>0,1000);
            GPSHz = round(1/(nanmedian(GPSDN(a+1)-GPSDN(a))*24*60*60));
%             posHz = round(1/(nanmedian(diff(pos{j}(1:101,1)))*24*60*60));
            posHz = Hzs.GPSHz;%10; % should be 10
%             disp('Used 10 Hz as GPS file.  If not correct change line 138 of addGPS script');
            dataHz = round(1/(nanmedian(diff(DATA{jj}.Time(1:101,1)))*24*60*60));
            databreak = [0; find(diff(DATA{jj}.Time+DATA{jj}.Date)>1.5/dataHz/24/60/60); length(DATA{jj}.Date)];
            %             oi = find(diff(GPSDN)>0.9*GPSHz/24/60/60);
            %             not50 = find(diff(oi)~=dataHz/GPSHz);
            if GPSHz<posHz %lots of data files were truncated to the floor of the second, so add the GPS sampling rate back in
                %                 ff = zeros(dataHz,1); for k = dataHz/posHz+1:dataHz/posHz:dataHz; ff(k:end) = ff(k:end)+1/posHz; end %times to add on
                %                 for k = 1:length(databreak)-1
                %                     a = find(~isnan(GPSDN(databreak(k)+1:databreak(k+1))),1)+databreak(k);
                %                     b = find(GPSDN(a:end)>GPSDN(a)+0.9/GPSHz/24/60/60,1);
                %                     b = min(a+b-1,databreak(k+1));
                %                     fff = [ff(end-(b-a)+1:end); ff(1:end-(b-a))];
                %                     fff = repmat(fff,ceil((databreak(k+1)-databreak(k))/dataHz),1);
                %                     GPSDN(a:databreak(k+1)) = GPSDN(a:databreak(k+1)) + fff(1:databreak(k+1)-a+1)/24/60/60;
                %                 end
                
                oi = find(diff(GPSDN)~=0);
                oi = [find(~isnan(GPSDN),1)-1; oi; length(GPSDN)];
                asec = floor((0:1/dataHz:1-1/dataHz)*posHz)'/posHz;
                for k = 1:length(oi)-1
                    dif = oi(k+1)-oi(k);
                    if dif >= dataHz
                        onesec = floor((0:1/dif:1-1/dif)*posHz)'/posHz;
                    else
                        onesec = asec(end-dif+1:end);
                    end
                    GPSDN(oi(k)+1:oi(k+1)) = GPSDN(oi(k)+1:oi(k+1))+onesec/24/60/60;
                end
            end
            tic
            POS = nan(length(pos{j}(:,1))*dataHz/posHz,length(pos{j}(1,:)));
            for k = 1:dataHz/posHz
                POS(k:dataHz/posHz:end,:) = pos{j};
            end
            bins = min(GPSDN)-.5/posHz/24/60/60:1/posHz/24/60/60:max(GPSDN)+.5/posHz/24/60/60;
            [~,posI] = histc(pos{j}(:,1),bins);
            [~,dataI] = histc(GPSDN,bins);
            posI(posI==0) = -1; % ensures that if pos times don't fit in they are ignored
            [~,Ia,Ib] = intersect(dataI,posI);
%             Ib(Ib==0) = -1; 
%             if any(Ib==0); error('why do the times from the pos file not fit into the times from data file?'); end
            k =0;
            while sum(Ia)~=0
                DATA{jj}.Lat(Ia) = pos{j}(Ib,2); DATA{jj}.Long(Ia) = pos{j}(Ib,3); DATA{jj}.GPSsd(Ia,:) = pos{j}(Ib,[7 8 10]);  DATA{jj}.Elev(Ia) = pos{j}(Ib,4); 
                dataI(Ia) = 0;
                [~,Ia,Ib] = intersect(dataI,posI);
                k = k+1;
                if k>50; error('seems to be an endless while loop...'); end
%                   [Ia(1:10) dataI(Ia(1:10))]
%                   sum(Ia)
            end
            data = DATA{jj};
            load([D{j} '\' datafile{jj}],'Adata','Atime');
            try load([D{j} '\' datafile{jj}],'ODN'); catch; end
            save([D{j} '\' datafile{jj}],'data','Adata','Atime');
            try save([D{j} '\' datafile{jj}],'ODN','-append'); catch; end
        end
       %             
        prhfile = cellfun(@(x) strfind(x,'prh.mat'),F{j},'uniformoutput',false);
        prhfile = F{j}(~cellfun(@isempty,prhfile));

%         
        % assumes that a data file you found with this pos file will relate
        % to the prh file
        for jj = 1:length(prhfile)
            load([D{j} '\' prhfile{jj}],'p','DN','fs','tagon','GPS','INFO');
            oi = GPS(1,:);
            timespan = 10;
            if ~exist('INFO','var'); timedif = input('timedif?  (from header xls file)'); INFO.timedif = timedif; end
            [GPS,GPSerr] = GPS2prh (DN,p,fs,tagon,data,prhfile{jj},50000,timespan,INFO.timedif);
            GPS(1,:) = oi;
            II = ~isnan(GPS(:,1)) & p>10;
            GPS(II,:) = nan;
            GPSerr(II,:) = nan;
%             [fig,ax] = plotMap(GPS,DN,tagon,10);
            surfs = findsurfacings(p,fs,tagon,60);
%             xs = get(ax,'xlim'); ys = get(ax,'ylim');
%             text(xs(1)+diff(xs)/20,ys(2)-diff(ys)/20,{['Number of surfacings: ' num2str(length(surfs))]; ['Number of GPS hits on surfacings: ' num2str(size(GPS(tagon&~isnan(GPS(:,1))),1))]},'parent',ax);
            report = [report; length(surfs) size(GPS(tagon&~isnan(GPS(:,1))),1)];
            save([D{j} '\' prhfile{jj}],'GPS','GPSerr','-append');
            oneup = max(strfind(D{j},'\'));
            if exist([D{j}(1:oneup) 'prh\'],'dir');
                save([D{j}(1:oneup) 'prh\' prhfile{jj}],'GPS','GPSerr','-append');
            end
            whalename = prhfile{jj}(1:strfind(prhfile{jj},' ')-1);
%             set(fig,'units','normalized','outerposition',[0 0 1 1]);
%             if ~exist([D{j} '\QL\'],'dir'); mkdir([D{j} '\QL\']); end
%             savefig(fig,[D{j} '\QL\' whalename ' Map.fig']);
%             saveas(fig,[D{j} '\' whalename ' Map.bmp']);
%             if exist([D{j}(1:oneup) 'Quicklook\'],'dir');
%                 saveas(fig,[D{j}(1:oneup) 'Quicklook\' whalename ' Map.bmp']);
%             end
            lat= GPS(~isnan(GPS(:,1)),1); lat(1) = [];
            long= GPS(~isnan(GPS(:,1)),2); long(1) = [];
            tstart = DN(~isnan(GPS(:,1)))-timespan/2/24/60/60; tstart(1) = [];
            tstop = DN(~isnan(GPS(:,1)))+timespan/2/24/60/60; tstop(1) = [];
            %convoluted system, but it works
            kmlstr = '';
            for ii = 1:length(lat)
                kmlstr = [kmlstr ge_point(long(ii),lat(ii),0,'timeSpanStart',datestr(tstart(ii)-UTC/24,'yyyy-mm-ddTHH:MM:SSZ'),'timeSpanStop',datestr(tstop(ii)-UTC/24,'yyyy-mm-ddTHH:MM:SSZ'))];
            end
            KML = ge_folder('Points',kmlstr);
            try
                ge_output([D{j}(1) ':\' whalename '.kml'],KML);
                movefile([D{j}(1) ':\' whalename '.kml'],[D{j} '\' whalename '.kml']); % some long directory names were giving problems to ge_output
            catch
                outdir = regexp(D{j},'\'); outdir = D{j}(1:outdir(4));
                ge_output([outdir whalename '.kml'],KML);
                movefile([outdir whalename '.kml'],[D{j} '\' whalename '.kml']); % some long directory names were giving problems to ge_output
            end
%             if exist([D{j}(1:oneup) 'Quicklook\'],'dir'); copyfile([D{j} '\' whalename '.kml'],[D{j}(1:oneup) 'Quicklook\' whalename '.kml']); end
        end
    end
    DATA = cell(0);
end
            
        
        
        
    
% end




