function [Lats, Longs] = track2latlong(GPS,tagon,geoPtrack)

Gi = find(~isnan(GPS(:,1))); [~,G0] = min(abs(Gi-find(tagon,1))); G1 = GPS(Gi(G0),:);  [x1,y1,z1] = deg2utm(G1(1),G1(2)); Lats = nan(size(tagon)); Longs = Lats; z = str2num(z1(1:2)); L = nan(1,2); L(2) = (z*6)-180; L(1) = L(2) - 6; startx = 0; starty = 0;
for i = find(tagon)'
[Lats(i),Longs(i)] = utm2deg(geoPtrack(i,1)+x1 - startx,geoPtrack(i,2)+y1 - starty,z1); 
if Longs(i)<L(1) || Longs(i) > L(2)
    z = ceil((Longs(i)+180)/6);
    z1 = [num2str(z) z1(end-1:end)];
    [x1,y1,z1] = deg2utm(Lats(i),Longs(i));
    startx = geoPtrack(i,1); starty = geoPtrack(i,2);
    [Lats(i),Longs(i)] = utm2deg(geoPtrack(i,1)+x1 - startx,geoPtrack(i,2)+y1 - starty,z1); 
    L(2) = (z*6)-180; L(1) = L(2) - 6;
end
end