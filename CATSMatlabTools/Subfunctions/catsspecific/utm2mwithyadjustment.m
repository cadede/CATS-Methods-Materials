function [xp,yp,UZ] = utm2m(x,y,uz)


oi = cellstr(uz);
[C,~,ic] = unique(oi);
UZ = uz(1,:);%C{mode(ic)};
xp = nan(length(x),1); yp = xp;
for i = 1:length(x)
%     if strcmp(uz(i,:),UZ); continue; end
    if strcmp(uz(i,1:2),UZ(1:2)); xp(i) = x(i); yp(i) = y(i); continue; end %if crosses E-W zones
    umain = str2num(UZ(1:2));
    unew = str2num(uz(i,1:2));
    [Lat,Long] = utm2deg(x(i),y(i),uz(i,:));
    [~,Longmain] = utm2deg(500000,y(i),UZ);
    [~,Longthis] = utm2deg(500000,y(i),uz(i,:)); %gives central meridian
    Longmain = [Longmain-3 + .000000001 Longmain+3-.000000001];
    Longthis = [Longthis-3 + .000000001 Longthis+3-.000000001];
    mainx = nan(1,2); mainy=mainx;
    [mainx(1),mainy(1)] = deg2utm(Lat,min(Longmain));
    [mainx(2),mainy(2)] = deg2utm(Lat,max(Longmain)); % mainy should be the same on both sides of the zone
    thisy = nan(1,2);
    [~,thisy(1)] = deg2utm(Lat,min(Longthis));
    [~,thisy(2)] = deg2utm(Lat,max(Longthis)); % mainy sh
    if diff(mainy)>0.1 || diff(thisy)>.1; error('utm conversion not giving accurate results in this zone'); end
    yp(i) = y(i)-thisy(2)+mainy(2)-2*(y(i)-thisy(2));
  
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

