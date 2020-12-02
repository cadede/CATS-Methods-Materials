function combos = getcombos(vidDN,vidDurs,viddeploy,maxL,maxgap)

% combine all videos < 5 minutes so long as there is less than a minute
% between them.  Have three seconds of black in between but show the whole
% graph.

% if nargin<3 || isempty(minL); minL = 300; end
if nargin<4 || isempty(maxL); maxL = 30*60; end
if nargin<5; maxgap = 300; end

i = find(~isnan(vidDN)&~isnan(vidDurs),1);
oi = find(isnan(vidDN)); oi(oi<i) = [];
vidDN(oi) = vidDN(oi-1);
vidDurs(oi) = vidDurs(oi-1);

% shorts = find(vidDurs<minL);
v1 = viddeploy(1);
combos = {v1};
n = 1;
L = vidDurs(v1);
for i = 2:length(viddeploy)
    if isnan(vidDN(viddeploy(i)))||isnan(vidDurs(viddeploy(i))); continue; end
    gap = (vidDN(viddeploy(i))-vidDN(viddeploy(i-1)))*24*60*60-vidDurs(viddeploy(i-1));
    if gap>maxgap || L+gap+vidDurs(viddeploy(i))>maxL 
%     if (vidDurs(i)<minL || vidDurs(i-1)<minL) && <maxL
        n = n+1; combos{n,1} = viddeploy(i);
        L = vidDurs(viddeploy(i));
    else
        combos{n} = [combos{n} viddeploy(i)];
        L = L + gap + vidDurs(viddeploy(i));
    end
end