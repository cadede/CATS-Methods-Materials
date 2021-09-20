function y = runcirc_mean (x,m) %input in radians, not fast

if sum(size(x)>1)<2
X = buffer(x,2*m,2*m-1,x(1)*ones(2*m-1,1));
X = X(:,m:end);
X2 = rot90(buffer(flipud(x(end-(2*m-1):end)),2*m,2*m-1,x(end)*ones(2*m-1,1)),2);
X = [X X2(:,2:m)];
y = circ_mean(X);
if size(x,1)>size(x,2); y = y'; end
else
    y = nan(size(x));
    for i = 1:size(x,2)
        x2 = x(:,i);
        X = buffer(x2,2*m,2*m-1,x2(1)*ones(2*m-1,1));
        X = X(:,m:end);
        X2 = rot90(buffer(flipud(x2(end-(2*m-1):end)),2*m,2*m-1,x2(end)*ones(2*m-1,1)),2);
        X = [X X2(:,2:m)];
        y(:,i) = circ_mean(X)';
    end
end


