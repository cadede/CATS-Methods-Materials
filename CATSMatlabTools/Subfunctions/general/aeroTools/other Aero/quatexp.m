function qout = quatexp(q)
%  QUATEXP Calculate the exponential of a quaternion.
%   QE = QUATEXP( Q ) calculates the exponential, QE, for a given
%   quaternion, Q.
%
%   Input:
%   Q:      M-by-4 matrix containing M quaternions.
%
%   Output:
%   QE:     M-by-4 matrix containing exponentials of quaternions.
%
%   Examples:
%
%   Determine the exponential of q = [0 0 0.7854 0]:
%       qe = quatexp([0 0 0.7854 0])
%
%   See also QUATCONJ, QUATDIVIDE, QUATINV, QUATMOD, QUATMULTIPLY,
%   QUATNORMALIZE, QUATROTATE, QUATLOG, QUATPOWER, QUATINTERP.

%   Copyright 2015-2018 The MathWorks, Inc.

% Validate input
validateattributes(q,{'numeric'},{'ncols',4,'real','finite','nonnan'})

% Calculate half the rotation
len = size(q,1);
th = arrayfun(@(k) norm(q(k,2:4)),1:len,'UniformOutput',true);
% Calculate exponential
qout = reshape(cell2mat(arrayfun(@(k) exp(q(k,1))*[cos(th(k)) sin(th(k))*q(k,2:4)/th(k)],1:len,'UniformOutput',false)'),len,4);
% Replace NaN's by 0's for singularity
if any(th==0)
    rZero = zeros(len,3);
    qout(th==0,2:4)=rZero(th==0,:);
end