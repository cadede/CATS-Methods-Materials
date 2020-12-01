function surfs = findsurfacings(p,fs,tagon,period,thresh)
% returns an index of times when there is a surfacing
% period is chunks of time to look (in minutes) at to 0 out the depth
% tagon is an index of times the tag is on the animal
%
if nargin<5; thresh = 1; end
if nargin<4; period = 60; end
if nargin<3; t1 = find(p>3,1,'first'); t2 = find(p>3,1,'last');
else t1 = find(tagon,1); t2 = find(tagon,1,'last'); end
p = runmean(p,fs);

surfs = [];
for i = t1:period*fs*60:t2
    minp = min(p(i:min(i+period*60*fs-1,t2)));
    [ps, hs]= peakfinder(p(i:min(i+period*fs*60-1,t2)),0.5,minp+thresh,-1);
    oi = find(diff(ps)<10*fs);
    for ii = 1:length(oi)
        [~,b] = max(hs(oi(ii):oi(ii)+1));
        ps(oi(ii)-1+b(1)) = nan;
    end
    ps = ps+i-1;
    surfs = [surfs; ps];
end