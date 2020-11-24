function [y,z] = consec(x)
%CONSEC		Run-length encode consecutive integers as (start,end) runs.
%
% y = consec(x)
%   Given a vector of numbers x that includes sequences of consecutive
%   ascending integers, encode it as (start,end) pairs, where each pair
%   stands for one sequence of consecutive numbers.  Each pair forms a
%   column of the return value y -- that is, y is a 2-by-N array.
%
%   Example:
%          consec([4 5 6 17 18 2 11 12 13]) ==> [4 17 2 11]
%                                               [6 18 2 13]
%
% [starts,stops] = consec(x)
%   As above, but return two 1-by-N vectors instead of one 2-by-N array.

n = length(x);
if (n)
  s = find(x(1:n-1) ~= x(2:n)-1);
  if (nCols(x) == 1)
    y = x([1 s'+1])';			% x is a column vector
    z = x([s' n])';
  else
    y = x([1 s+1]);			% x is a row vector
    z = x([s n]);
  end
else
  y = [];
  z = [];
end

if (nargout < 2)
  y = [y; z];
end
