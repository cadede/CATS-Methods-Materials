function qout = quatinterp(p,q,f,varargin)
%  QUATINTERP Calculate the interpolation between two quaternions.
%   QI = QUATINTERP(P,Q,F,METHOD) calculates the interpolation quaternion
%   between two normalized quaternions P and Q by interval fraction F.
%
%   Inputs:
%   P,Q:        M-by-4 matrices containing M quaternions.
%   F:          M-by-1 vector. F varies between 0 and 1.
%   METHOD:     String specifying the method for the interpolation. The
%               options for the methods are: 
%               'SLERP' - (Default) Spherical linear quaternion 
%                         interpolation. 
%               'LERP'  - Linear quaternion interpolation.
%               'NLERP' - Normalized linear interpolation.
%
%   Output:
%   QI:         M-by-4 matrix of interpolated quaternions.
%
%   Examples:
%
%   Determine the quaternion located in the middle between p = [0.7071 0 0.7071 0]
%   and q = [-0.7071 0 0.7071 0] using the 'SLERP' method:
%
%       pn = quatnormalize([0.7071 0 0.7071 0])
%       qn = quatnormalize([-0.7071 0 0.7071 0])
%       qi = quatinterp(pn,qn,0.5,'slerp')
%
%   See also QUATCONJ, QUATDIVIDE, QUATINV, QUATMOD, QUATMULTIPLY, 
%   QUATNORMALIZE, QUATROTATE, QUATLOG, QUATEXP, QUATPOWER.

%   Copyright 2015 The MathWorks, Inc.

% Validate input number
nargoutchk(0,1)

% Validate output number
narginchk(3,4)

% Validate inputs
validateattributes(p,{'single','double'},{'ncols',4,'real','finite','nonnan'});
validateattributes(q,{'single','double'},{'ncols',4,'nrows',size(p,1),'real','finite','nonnan'});
validateattributes(f,{'single','double'},{'real','finite','nonnan','>=',0,'<=',1,'numel',size(p,1)});
if length(varargin)==1
    method = lower(varargin{1});
    method = validatestring(method,{'slerp','nlerp','lerp'});
else
    method = 'slerp';
end
% In case the quaternions are not normalized, then normalize them
modp = quatmod(p);
modqq = quatmod(q);
nflagp = logical((modp>1.0+sqrt(eps)) + (modp<1.0-sqrt(eps)));
nflagq = logical((modqq>1.0+sqrt(eps)) + (modqq<1.0-sqrt(eps)));
warnflag = true;
if any(nflagp)
    p(nflagp,:) = quatnormalize(p(nflagp,:));
    warning(message('aero:quatlog:notUnitQuaternion'));
    warnflag = false;
end
if any(nflagq)
    q(nflagq,:) = quatnormalize(q(nflagq,:));
    if warnflag
        warning(message('aero:quatlog:notUnitQuaternion'));
    end
end
% In case that the interpolation is bigger than 90 degrees, one of the
% quaternions must be changed to its negative equivalent so the
% interpolation is symmetric.
dotpq = dot(p',q')<0;
if any(dotpq)
    q(dotpq,:) = -q(dotpq,:);
end
% Calculate quaternion interpolation depending on the selected method
qout = zeros(size(p));
switch method
    case 'slerp'
        qout = quatmultiply(p,quatpower(quatnormalize(quatmultiply(quatconj(p),q)),f));
    case 'lerp'
        for k=1:size(p,1)
            qout(k,:) = p(k,:)*(1-f(k))+q(k,:)*f(k);
        end
    case 'nlerp'
        for k=1:size(p,1)
            qout(k,:) = p(k,:)*(1-f(k))+q(k,:)*f(k);     
        end
        qout = quatnormalize(qout);
end
