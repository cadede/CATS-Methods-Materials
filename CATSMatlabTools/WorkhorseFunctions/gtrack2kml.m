function gtrack2kml (gtrack,tagon,fs,DN,kmlfs,startlat,startlong,UTC,whaleName,fileloc)

oi = false(size(tagon));
oi(fs/kmlfs:fs/kmlfs:end) = true;
[x,y,uzone] = deg2utm(startlat,startlong);
oiDN = DN(tagon&oi);
[lats,longs] = utm2deg(gtrack(tagon&oi,1)+x,gtrack(tagon&oi,2)+y,repmat(uzone,sum(tagon&oi),1));

tstart = oiDN-.5/24/60/60; tstop = oiDN+.5/24/60/60;
kmlstr = '';
for ii = 1:length(lats)
    kmlstr = [kmlstr ge_point(longs(ii),lats(ii),0,'timeSpanStart',datestr(tstart(ii)-UTC/24,'yyyy-mm-ddTHH:MM:SSZ'),'timeSpanStop',datestr(tstop(ii)-UTC/24,'yyyy-mm-ddTHH:MM:SSZ'))];
end
KML = ge_folder('Points',kmlstr);
outfile = [fileloc(1) ':\' whaleName '.kml'];
try ge_output(outfile,KML); catch; d = regexp(fileloc,'\'); outfile = [fileloc(1:d(3)) whaleName '.kml']; ge_output(outfile,KML); end
movefile(outfile,[fileloc '\' whaleName 'geoPtrack.kml'],'f'); % come long directory names were giving problems to ge