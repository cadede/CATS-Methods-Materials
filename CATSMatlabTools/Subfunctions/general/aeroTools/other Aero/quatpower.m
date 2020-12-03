function qout = quatpower(q,pow)
%  QUATPOWER Calculate the power of a quaternion.
%   QP = QUATPOWER( Q , POW ) raises Q to the power of POW, for a given
%   normalized quaternion, Q.
%
%   Inputs:
%   Q:      M-by-4 matrix containing M quaternions.
%   POW:    M-by-1 vector of exponents, specified as scalars.
%
%   Output:
%   QP:     M-by-4 matrix of exponentiated quaternions.
%
%   Examples:
%
%   Determine the power of 2 of q = [0.7071 0 0.7071 0]:
%       qp = quatpower(quatnormalize([0.7071 0 0.7071 0]),2)
%
%   See also QUATCONJ, QUATDIVIDE, QUATINV, QUATMOD, QUATMULTIPLY, 
%   QUATNORMALIZE, QUATROTATE, QUATLOG, QUATEXP, QUATINTERP.

%   Copyright 2015 The MathWorks, Inc.

% Validate inputs
validateattributes(q,{'numeric'},{'ncols',4,'real','finite','nonnan'})
validateattributes(pow,{'numeric'},{'real','finite','nonnan'})

% Calculate modulus
modq = quatmod(q);

% Check that the quaternions are normalized, normalize if they are not.
nflag = logical((modq>1.0+sqrt(eps))+(modq<1.0-sqrt(eps)));
if any(nflag)
    q(nflag,:) = quatnormalize(q(nflag,:));
    warning(message('aero:quatlog:notUnitQuaternion'));
end

% Calculate power
len = size(q,1);
qout = cell2mat(arrayfun(@(k) quatexp(pow(k)*quatlog(q(k,:))),1:len,'UniformOutput',false)');