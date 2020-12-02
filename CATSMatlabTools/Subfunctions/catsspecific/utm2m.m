function [xp,UZ] = utm2m(x,y,uz)


oi = cellstr(uz);
[C,~,ic] = unique(oi);
UZ = C{mode(ic)};
xp = nan(length(x),1);
for i = 1:length(x)
%     if strcmp(uz(i,:),UZ); continue; end
    if strcmp(uz(i,1:2),UZ(1:2)); xp(i) = x(i); continue; end %if crosses E-W zones
    umain = str2num(UZ(1:2));
    unew = str2num(uz(i,1:2));
    [Lat,Long] = utm2deg(x(i),y(i),uz(i,:));
    [~,Longmain] = utm2deg(500000,y(i),UZ);
    Longmain = [Longmain-3 + .000000001 Longmain+3-.000000001];
    mainx = [deg2utm(Lat,min(Longmain)) deg2utm(Lat,max(Longmain))];
    if unew>umain; newx = max(mainx); else newx = min(mainx); end
    for k = umain+1:unew
        if k<10; str = ['0' num2str(k)]; else str = num2str(k); end
        [~,Longthis] = utm2deg(500000,y(i),[str uz(i,end-1:end)]);
        Longthis = [Longthis-3 + .000000001 Longthis+3-.000000001];
        thisx = [deg2utm(Lat,min(Longthis)) deg2utm(Lat,max(Longthis))];
        if k~=unew; newx = newx + diff(thisx); continue; end
        newx = newx + (x(i)-min(thisx));
    end
    for k = umain-1:-1:unew
        if k<10; str = ['0' num2str(k)]; else str = num2str(k); end
        [~,Longthis] = utm2deg(500000,y(i),[str uz(i,end-1:end)]);
        Longthis = [Longthis-3 + .000000001 Longthis+3-.000000001];
        thisx = [deg2utm(Lat,min(Longthis)) deg2utm(Lat,max(Longthis))];
        if k~=unew; newx = newx - diff(thisx); continue; end
        newx = newx - (max(thisx) - x(i));
    end
    xp(i) = newx;
end

