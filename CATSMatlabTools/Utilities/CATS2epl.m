% function CATS2epl(geoPtrack, DN, fs, speedJJ, tagon, ID, FS,filedest) %,starttime, endtime 

% uncomment aboe to run this script like a function. This inputs data from a prh
% file and outputs an epl file suitable for plotting animal tracks in
% echoview, including feeding events.

% for n = [1:4 6:length(D)];
% set filedest where your prh file is and filedest2 where you want the
% final epl file
    clearvars -except D n
filedest = 'Y:\Gold Data\Tag Data\CATS\tag_data\bw180827-52 (IOS_Monterey)\';% ['E:\CATS\tag_data\' D{n} '\'];
filedest2 = 'Y:\Gold Data\Echosounders\Monterey\2018\epl files\';
% Istart = 102878;%120000; % to truncate file to a certain segment, set these to the indices of your prh file.
% Iend = 138000;
Tstart = [2018 08 27 09 00 00 ]; % time start and end of prey mapping around tagged animal (in local time)
    Tend = [2018 08 27 14 00 00]; 


D2 = dir(filedest); D2 = {D2.name}; prhfile = D2{~cellfun(@isempty,cellfun(@(x) strfind(x,'prh.mat'),D2,'uniformoutput',false))};
try
    lungefile = D2{~cellfun(@isempty,cellfun(@(x) strfind(x,'lunges.mat'),D2,'uniformoutput',false))};
    if isempty(lungefile); D2 = dir([filedest 'lunges\']); D2 = {D2.name}; lungefile = ['lunges\' D2{~cellfun(@isempty,cellfun(@(x) strfind(x,'lunges.mat'),D2,'uniformoutput',false))}]; end
    load([filedest lungefile]);
catch
    warning('No lunge file detected')
end
load([filedest prhfile]);

[~,Istart] = min(abs(DN-(datenum(Tstart)-30/60/24))); Istart = max(find(tagon,1),Istart); % start the epl file 30 minutes before prey mapping starts (or tag on)
[~,Iend] = min(abs(DN-(datenum(Tend)+30/60/24))); Iende = min(find(tagon,1,'last'),Iend); % start the epl file 30 minutes before prey mapping starts (or tag on)
if ~isinf(Iend) || Istart~=1 %tail = ['-' num2str(Istart) '-' num2str(Iend)];
    T1 = round(datevec(DN(Istart))); T2 = round(datevec(DN(Iend)));
    tail = ['-' sprintf('%02d',T1(4),T1(5),T1(6)) '-' sprintf('%02d',T2(4),T2(5),T2(6))];
    trunc = true;
else tail = ''; trunc = false;
end
FS = 1;
df = round(fs/FS);  if FS>fs || abs(round(fs/FS)-fs/FS)> .001; error('Choose an FS smaller than fs that divides evenly'); end
% heads = {'Ping_date' 'Ping_Time' 'Ping_milliseconds' 'Latitude' 'Longitude' 'Position_status' 'Depth' 'Line_status' 'Ping_status' 'Altitude' 'GPS_UTC_time'};
p2 = decdc(p,df);
DNdec = DN(1:df:end);

c1 = cellstr(datestr(DN(1:df:end)-INFO.UTC/24,'yyyy-mm-dd'));
c2 = cellstr(datestr(DN(1:df:end)-INFO.UTC/24,'HH:MM:SS'));
c3 = cellstr(datestr(DN(1:df:end)-INFO.UTC/24,'fff'));
Gi = find(~isnan(GPS(:,1))); [~,G0] = min(abs(Gi-find(tagon,1))); G1 = GPS(Gi(G0),:);  [x1,y1,z1] = deg2utm(G1(1),G1(2)); [Lats,Longs] = utm2deg(geoPtrack(tagon,1)+x1,geoPtrack(tagon,2)+y1,repmat(z1,sum(tagon),1)); lats = nan(size(tagon)); longs = lats; lats(tagon) = Lats; longs(tagon) = Longs;
c4 = lats(1:df:end);
c5 = longs(1:df:end);
for i = 1:5; eval(['c' num2str(i) ' = c' num2str(i) '(1:length(p2));']); end
c6 = p2;
sp = speed.JJ; sp(p<2) = nan; sp = fixgaps(sp); sp(isnan(sp)) = 0; 
if ~exist('LungeI','var')
    c7 = decdc(sp,df);
else
    c7 = zeros(size(p));
    try llunge = LungeI(LungeC>=2); catch; llunge = LungeI; end
    c7(llunge) = 1;
    speedinc = [false; diff(sp)>0;];
    for ii = 1:length(llunge)
        % to time directly with mouth open, use the next line, else for 20
        % s around lunge, use the line after that
%         c7(llunge(ii):find(speedinc(llunge(ii)+fs/2:end),1)+llunge(ii)+fs/2-1) = 1; % find the first time speed increases greater than a half second after the mouth opening (essentially find mouth closed).
        c7(round(llunge(ii)-10*fs):round(llunge(ii)+10*fs)) = 1;
    end
    c7 = c7(1:df:end);
    if length(c7) == length(c4)+1; c7(end) = []; end
end
c8 = ones(size(c4));
T = table(c1,c2,c3,c4,c5,c6,c7,c8);
tagondec = tagon(1:df:end); tagondec(end) = [];
if trunc
    T = T(ceil(Istart/df):min(floor(Iend/df),size(T,1)),:);
    else T = T(tagondec,:);
end
% DNdec = DNdec(tagondec);
% [~,a] = min(abs(DNdec-starttime));
% [~,b] = min(abs(DNdec-endtime));
% T = T(a:b,:);
%

% csvwrite(,[heads; table2cell(T)]);
% dlmwrite('B:\Dropbox\data\Supergroup\Gold foraging scaling\bwdive.csv',table2cell(T),',');
writetable(T,[filedest INFO.whaleName tail '.epl.txt'],'delimiter',',')
copyfile([filedest INFO.whaleName tail '.epl.txt'],[filedest2 INFO.whaleName tail '.epl']);
disp([filedest2 INFO.whaleName tail '.epl written and copied to ' filedest])
% end