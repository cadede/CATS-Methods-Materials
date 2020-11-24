function filtered = filterCATS(X,k,m,thresh)
% removes data spikes by looking for large values in the difference between
% consecutive values.  Replaces new nans with averages of the values in
% the surrounding column.  k is the max distance (in points) between the
% start of a spike and the end of a spike (versus a jump to a "real" 
% dataset of points). Uses a runmean filter of size m to find the "true
% max" and "true min" of the data.  Thresh is the percentage of the max-min
% dif to set as the threshold for which to look for "large" differences

wasnan = isnan(X);

if nargin < 2
    m = 10;
    k = 5;
    thresh = .5;
elseif nargin < 3
    m = 10;
    thresh = .5;
elseif nargin < 4
    thresh = .5;
end

if size(X,1) == 1; X = X'; T = true; else T = false; end

Xr = X;
for j = 1:size(X,2)
    Xr(:,j) = runmean(fixgaps(X(:,j)),m);
end
M = max(Xr);
mi = min(Xr);
baddiff = (M-mi)*thresh;
for j = 1:length(baddiff)
    poses = find(diff(X(:,j))>baddiff(j))';
    negs = find(diff(X(:,j))<-baddiff(j))';
%     pairs = repmat(poses,2,1);
%     todel = [];
    pairs = nan(2,0);
    for n = 1:length(poses)
        CLOSE = find(abs(negs-poses(n))<=k);
        for nn = 1:length(CLOSE)
            pairs(:,end+1) = [min(poses(n),negs(CLOSE(nn))); max(poses(n),negs(CLOSE(nn)))];
        end
    end
    for n = 1:length(negs)
        CLOSE = find(abs(poses-negs(n))<=k);
        for nn = 1:length(CLOSE)
            pairs(:,end+1) = [min(negs(n),poses(CLOSE(nn)));  max(negs(n),poses(CLOSE(nn)))];
        end
    end

    for n = 1:length(pairs(1,:))
        X(pairs(1,n)+1:pairs(2,n),j) = nan;
    end
    
    X(:,j)= fixgaps(X(:,j));
end

if T; filtered = X'; else filtered = X; end
filtered(wasnan) = nan;

