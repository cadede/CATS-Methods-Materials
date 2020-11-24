function [W,calperiodI,tagslip,tocal,startsI,endsI,Wchange,Wchangeend] = reconcileSlips(W,calperiodI,tagslip)

% needs W,calperiodI,slips
% if there are fewer slips, keeps the first calperiod that exists from the
% section before and after the slip that was removed.

startsI = tagslip(1:end-1,2);
endsI = tagslip(2:end,1);
Wchange = endsI(1:end-1); Wchangeend = startsI(2:end);
Wtemp = W;
ctemp = calperiodI;
W = cell(size(startsI));
calperiodI = W;
if ~isempty(ctemp)&&length(Wtemp)~=length(startsI) % if there is already a W for every tag slip, no need to do this.  if calperiodI is empty, can't do this
    for i = 1:length(startsI)
        iscal = cellfun(@(x) all(x>startsI(i) & x<endsI(i))&~isempty(x),ctemp);
        if sum(iscal) == 1
            W{i} = Wtemp{iscal};
            calperiodI{i} = ctemp{iscal};
        elseif sum(iscal) == 0;
            W{i} = []; calperiodI{i} = [];
        elseif ~isempty(Wtemp{find(iscal,1)});
            W{i} = Wtemp{find(iscal,1)};
            calperiodI{i} = ctemp{find(iscal,1)};
        else Wdata = ~cellfun(@isempty,Wtemp);
            ii = find(Wdata&iscal,1);
            if ~isempty(ii); W{i} = Wtemp{ii}; calperiodI{i} = ctemp{ii}; end
        end
    end
else W = Wtemp; calperiodI = ctemp;
end
tocal = [startsI endsI]; % periods of deployment left to calibrate
tocal(~cellfun(@isempty, W),:) = []; % if W already exists, do not need to calibrate that period

