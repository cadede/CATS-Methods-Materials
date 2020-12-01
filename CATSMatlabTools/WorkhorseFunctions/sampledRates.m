function [accHz,gyrHz,magHz,pHz,lHz,GPSHz,UTC,THz,T1Hz] = sampledRates(fileloc,file)

% this looks for the txt file that has information about sample rates
% sometimes txt file doesn't exactly match csv filename with a slightly different time, so find
% the one that matches the best
try
    % first try to automatically import text file
    try
        FILES = dir(fileloc); FILES ={FILES.name};
        FILEST= FILES(~cellfun(@isempty, strfind(FILES,[file(1:end-4) '.txt'])));
        if ~isempty(FILEST); FILES = FILEST; else  FILES= FILES(~cellfun(@isempty, strfind(FILES,[file(1:end-8) '.txt']))); end
        [~,b] = max(cellfun(@(x) sum(x(1:min(length(x),length(file))) == file(1:min(length(x),length(file)))),FILES));
        FILES = char(FILES(b));
        Hzs = importdata([fileloc FILES],'\t'); try Hzs = Hzs.textdata; catch; end %old instead of FILES: fname(1:end-3) 'txt'
    catch
        try
            dirup = regexp(fileloc,'\');
            dirup = fileloc(1:dirup(end-1));
            FILES = dir(dirup); FILES ={FILES.name};
            FILES= FILES(~cellfun(@isempty, strfind(FILES,[file(1:end-4) '.txt'])));
            [~,b] = max(cellfun(@(x) sum(x(1:min(length(x),length(file))) == file(1:min(length(x),length(file)))),FILES));
            FILES = char(FILES(b));
            Hzs = importdata([dirup FILES],'\t'); try Hzs = Hzs.textdata; catch; end %old instead of FILES: fname(1:end-3) 'txt'
        catch
            FILES = dir([fileloc 'raw\']); FILES ={FILES.name};
            FILES= FILES(~cellfun(@isempty, strfind(FILES,[file(1:end-4) '.txt'])));
            [~,b] = max(cellfun(@(x) sum(x(1:min(length(x),length(file))) == file(1:min(length(x),length(file)))),FILES));
            FILES = char(FILES(b));
            Hzs = importdata([fileloc 'raw\' FILES],'\t'); try Hzs = Hzs.textdata; catch; end %old instead of FILES: fname(1:end-3) 'txt'
        end
    end
catch
    [txtfile,txtfileloc] = uigetfile('*.txt','Select meta data txt file with sample rate information');
    Hzs = importdata([txtfileloc txtfile],'\t'); try Hzs = Hzs.textdata; catch; end %old instead of FILES: fname(1:end-3) 'txt'
end
try UTC = Hzs{find(~cellfun(@isempty, strfind(Hzs,'utc')),1,'first')};
    UTC = str2num(UTC(regexp(UTC,'=')+1:end));
catch
    UTC = false;
end
ints = find(~cellfun(@isempty, strfind(Hzs,'interval')));
Hzs = Hzs(ints(1)-3:end); ints = ints-(ints(1))+3+1;
% start reading and downsampling data
acc = find(~cellfun(@isempty,strfind(Hzs,'Accel')),1,'last'); acc = Hzs{ints(find(ints>acc,1,'first'))};
accHz =str2num(acc(regexp(acc,'=')+1:end));
acc = find(~cellfun(@isempty,strfind(Hzs,'Gyro')),1,'last'); acc = Hzs{ints(find(ints>acc,1,'first'))};
gyrHz =str2num(acc(regexp(acc,'=')+1:end));
try acc = find(~cellfun(@isempty,strfind(Hzs,'Compass')),1,'last'); acc = Hzs{ints(find(ints>acc,1,'first'))};
catch; acc = find(~cellfun(@isempty,strfind(Hzs,'Magnet')),1,'last'); acc = Hzs{ints(find(ints>acc,1,'first'))};
end
magHz =str2num(acc(regexp(acc,'=')+1:end));
try % To check again to see if there is a pressure sensor listed, else creating pressure with all 0s
    try acc = find(~cellfun(@isempty,strfind(Hzs,'Pressure')),1,'last'); acc = Hzs{ints(find(ints>acc,1,'first'))};
    catch; acc = find(~cellfun(@isempty,strfind(Hzs,'Depth')),1,'last'); acc = Hzs{ints(find(ints>acc,1,'first'))};
    end
    noPress = false;
catch
    warning('Doesn''t appear to be a pressure sensor.  Inputting 0 for pressure'); noPress = true;
end
try pHz =str2num(acc(regexp(acc,'=')+1:end)); catch; pHz = 10; end
acc = find(~cellfun(@isempty,strfind(Hzs,'Light')),1,'last'); acc = Hzs{ints(find(ints>acc,1,'first'))};
lHz =str2num(acc(regexp(acc,'=')+1:end));
try acc = find(~cellfun(@isempty,strfind(Hzs,'GPS')),1,'last'); acc = Hzs{ints(find(ints>acc,1,'first'))}; GPSHz =str2num(acc(regexp(acc,'=')+1:end)); catch; GPSHz = nan; end
try
acc = find(~cellfun(@isempty,strfind(Hzs,'Temp.')),1,'last'); acc = Hzs{ints(find(ints>acc,1,'first'))};
THz =str2num(acc(regexp(acc,'=')+1:end));
catch; THz = nan;
end
try
acc = find(~cellfun(@isempty,strfind(Hzs,'Temperature')),1,'last'); acc = Hzs{ints(find(ints>acc,1,'first'))};
T1Hz =str2num(acc(regexp(acc,'=')+1:end));
catch; T1Hz = nan;
end
