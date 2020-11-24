function [Aw,Mw,Gw] = applyW(W,startsI,endsI,At,Mt,Gt,skipempty)
% startsI and endsI are the start and end indices of where to apply W
% W is a cell, each of which has a 3x3 rotation matrix
% if skipempty is true, only calculates whale frame for the cells that exist, else it throws an error for missing values

% whale frame variables come out with nans where W is not applied. (e.g.
% off whale)

if ~exist('skipempty','var') || isempty(skipempty); skipempty = false; end
if length(W)~=length(startsI); error('check length W'); end
Aw = nan(size(At)); Mw = Aw; Gw = Aw;
% tagprh = nan(length(W),3);

for i = 1:length(startsI);
    if isempty(W{i}) && skipempty; continue;
    elseif isempty(W{i}); error(['W is empty in calperiod ' num2str(i)]);
    end
    I = startsI(i):endsI(i);
    [Aw(I,:),Mw(I,:),Gw(I,:)] = tagframe2whaleframe(At(I,:),Mt(I,:),Gt(I,:),[],[],[],W{i}); % calperiodI{i}-startsI(i)+1
    if i~=1 && ~isempty(W{i-1}) && endsI(i-1)~=startsI(i)
        I = endsI(i-1):startsI(i);
        [Aw(I,:),Mw(I,:),Gw(I,:)] = applyWduringslip(At(I,:),Mt(I,:),Gt(I,:),W{i-1},W{i});
    end
    if i~=length(W) && ~isempty(W{i+1}) && endsI(i)~=startsI(i+1)
        I = endsI(i):startsI(i+1);
        [Aw(I,:),Mw(I,:),Gw(I,:)] = applyWduringslip(At(I,:),Mt(I,:),Gt(I,:),W{i},W{i+1});
    end
end