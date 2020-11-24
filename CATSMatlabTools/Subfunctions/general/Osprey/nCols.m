function n = nCols(array)
%   nCols(array) returns the number of columns in the array.  This is the
%   size of the second dimension, i.e., it is shorthand for size(array,2).  
%   For n-dimensional arrays, this is different from the 'n' value given
%   by "[m,n] = size(array)", as 'size' with two output arguments rolls
%   together all the dimensions after the first one.
%
% See also nRows, size.

n = size(array,2);
