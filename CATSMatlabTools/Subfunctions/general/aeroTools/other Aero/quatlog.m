function qout = quatlog(q)
%  QUATLOG Calculate the natural logarithm of a quaternion.
%   QL = QUATLOG( Q ) calculates the natural logarithm, QL, for a given
%   normalized quaternion, Q. 
%
%   Input: 
%   Q:      M-by-4 matrix containing M quaternions. 
% 
%   Output:
%   QL:     M-by-4 matrix of quaternion natural logarithms.
%
%   Examples:
%
%   Determine the natural logarithm of q = [0.7071 0 0.7071 0]:
%       qlog = quatlog(quatnormalize([0.7071 0 0.7071 0]))
%
%   See also QUATCONJ, QUATDIVIDE, QUATINV, QUATMOD, QUATMULTIPLY, 
%   QUATNORMALIZE, QUATROTATE, QUATEXP, QUATPOWER, QUATINTERP.

%   Copyright 2015 The MathWorks, Inc.

% Validate input
validateattributes(q,{'numeric'},{'ncols',4,'real','finite','nonnan'});

% Calculate modulus
modq = quatmod(q);

% In case the quaternions are not normalized, then normalize them
nflag = logical((modq>1.0+sqrt(eps))+(modq<1.0-sqrt(eps)));
if any(nflag)
    q(nflag,:) = quatnormalize(q(nflag,:));
    warning(message('aero:quatlog:notUnitQuaternion'));
end

% Calculate half the rotation angle
len = size(q,1);
normv = arrayfun(@(k) norm(q(k,2:4)),1:len,'UniformOutput',true)';
th = atan2(normv,q(:,1));

% Initialize outputs
qout = zeros(size(q));

% Calculate logarithm
tmp = arrayfun(@(k) th(k)*q(k,2:4)/normv(k),1:len,'UniformOutput',false);
tmp = reshape(cell2mat(tmp'),length(modq),3);
qout(normv~=0,2:4)=tmp(normv~=0,:);