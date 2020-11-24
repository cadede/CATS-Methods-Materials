function outs = edgenans(ins)
outs = ins;
if any(sum(isnan(ins))>0)
    for j = find(sum(isnan(ins))>0)
        I = find(isnan(ins(:,j)));
        [s,e] = consec(I);
        m = round(mean([s; e]));
        if s(1) == 1; s(1) = 2; m(1) = 1; end;
        if e(end) == size(ins,1); m(end) = s(end); e(end) = e(end)-1; end
        for i = 1:length(s)
            outs(s(i):m(i),j) = nanmean(outs(max(1,s(i)-5):s(i)-1,j)); %replaces nans with the closest actual value
            outs(m(i):e(i),j) = nanmean(outs(e(i)+1:min(size(ins,1),e(i)+5),j));
        end
    end
end