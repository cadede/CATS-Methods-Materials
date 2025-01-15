function [UTC,cityGPS,citytext,UTCs,TZs] = getUTC(lat,long,DN,UTCloc,cityGPS,citytext,UTCs,TZs)
%input latitude, longitude, datenumber (in local time), output a difference in UTC from the
%closest city with known timezone from UTC location spreadsheet
% set up to save time so that if you want it can input cityGPS and citytext
% so it doesn't have to reread it for every use.
%
if nargin<4 || isempty(UTCloc)
    try
        curdir = pwd;
        oi = strfind(curdir,'MATLAB');
        [D2,F2] = subdir(curdir(1:oi+5));
        oi = find(cellfun(@isempty,cellfun(@(x) strfind(x,'CATSMatlabTools\templates\'),D2,'uniformoutput',false))==0);
        Fnum = find(cellfun(@(y) y==1,cellfun(@(x) strcmp(x,'UTC location spreadsheet.xlsx'),F2{oi(1)},'uniformoutput',false)));
        fileloc = [D2{oi(1)} '\']; filename =  F2{oi(1)}{Fnum};
        [cityGPS,citytext]= xlsread([fileloc filename],'cities');
    catch
        fileloc = 'C:\Users\Dave\Documents\Programs\MATLAB\Tagging\CATS Tools\CATSMatlabTools\templates\';
        filename = 'UTC location spreadsheet.xlsx';
        try    [cityGPS,citytext]= xlsread([fileloc filename],'cities');
        catch
            [filename,fileloc] = uigetfile('*.xlsx','to calculate UTC offset select UTC location spreadsheet (should be in excel files folder)');
            [cityGPS,citytext]= xlsread([fileloc filename],'cities');
        end
    end
    [UTCs, TZs] = xlsread([fileloc filename],'DST');
elseif nargin<5
    [cityGPS,citytext]= xlsread([UTCloc 'UTC location spreadsheet.xlsx'],'cities');
    [UTCs, TZs] = xlsread([UTCloc 'UTC location spreadsheet.xlsx'],'DST');
end
[~,b] = min(arrayfun(@(x,y) distance(lat,long,x,y),cityGPS(:,1),cityGPS(:,2))); %find the closest city
TZ = citytext{b,end}; % timezone of closest city
TZnum = find(cellfun(@(x) strcmp(x,TZ),TZs(:,1)));
UTC = -UTCs(TZnum,4)/60; %reported in wrong direction from convention
DV = datevec(DN);
yy = DV(1);
DSToff = UTCs(TZnum,[8 11]); DSTon = UTCs(TZnum,[16 19]);
if sum(DSToff)~=0 && ~(sum(isnan(DSToff)) == 2);
    dayoff = max(find(weekday(datenum(yy,DSToff(1),1):datenum(yy,DSToff(1),1)+30)==1,DSToff(2)));
    dayon = max(find(weekday(datenum(yy,DSTon(1),1):datenum(yy,DSTon(1),1)+30)==1,DSTon(2)));
    dayoff = datenum(yy,DSToff(1),dayoff); dayon = datenum(yy,DSTon(1),dayon);
    isDST = (DN>dayon&&DN<dayoff)||(DN>dayon&&dayon>dayoff);
    if isDST; UTC = UTC - UTCs(TZnum,14)/60; end
end

