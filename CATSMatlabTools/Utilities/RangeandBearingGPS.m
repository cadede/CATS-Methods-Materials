function [newlat,newlong] = RangeandBearingGPS(lat,long,range,bearing)
%lat is decimal degrees
% long is decimal degrees
% range in m
% bearing in true deg from N

[x,y,utmzone] = deg2utm(lat,long);

offy = -sin((-bearing-90)*pi/180).*range;
offx = -cos((-bearing-90)*pi/180).*range;

[newlat,newlong] = utm2deg(x+offx,y+offy,utmzone);
