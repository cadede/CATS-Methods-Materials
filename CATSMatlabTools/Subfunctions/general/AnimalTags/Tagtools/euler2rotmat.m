function    Q = euler2rotmat(p,r,h)
%
%      Q = euler2rotmat(p,r,h)
%		 or
%      Q = euler2rotmat(prh)
%      Make a rotation (or direction cosine) matrix out of sets of Euler angles,
%		 pitch, roll, and heading.
%
%		 Inputs:
%		 p is the pitch angle in radians.
%		 r is the roll angle in radians.
%		 h is the heading or yaw angle in radians.
%		 p, r and h must either be the same size vectors (i.e., all sampled at the
%		 same rate) or one or more of p, r or h can be a scalar in which case the 
%		 same angle is used for all 
%	    p, r and h can also be given in a three column matrix: prh=[p,r,h].
%
%		 Returns:
%      Q contains one or more 3x3 rotation matrices. If p, r, and h are all scalars,
%		  Q is a 3x3 matrix, Q = H*P*R where H, P and R are the cannonical rotation 
%		  matrices corresponding to the yaw, pitch and roll rotations, respectively.
%	     To rotate a vector or matrix of triaxial measurements, pre-multiply by Q.
%		  If p, r or h contain multiple values, Q is a 3-dimensional matrix with
%		  size 3x3xn where n is the number of Euler angle triples that are input.
%		  To access the k'th rotation matrix in Q use squeeze(Q(:,:,k)).
%
%		 Example:
%		  Q = euler2rotmat([22 -22 85]*pi/180)
% 	     returns: Q = [0.080809  -0.911425  -0.403453
%   						 0.923656   0.220605  -0.313358
%   						 0.374607  -0.347329   0.859670]
%
%      Valid: Matlab, Octave
%      markjohnson@st-andrews.ac.uk
%      last modified: 15 May 2017


if nargin~=1 && nargin~=3,
   help euler_to_rotmat
   return
end

if nargin==1,
   h = p(:,3) ;
   r = p(:,2) ;
   p = p(:,1) ;
end

n = [length(p) length(r) length(h)] ;
nn = max(n) ;
if n(1)<nn, p = p(1)*ones(nn,1) ; end
if n(2)<nn, r = r(1)*ones(nn,1) ; end
if n(3)<nn, h = h(1)*ones(nn,1) ; end

cp = cos(p) ;
sp = sin(p) ;
cr = cos(r) ;
sr = sin(r) ;
ch = cos(h) ;
sh = sin(h) ;

Q = zeros(3,3,nn) ;
for k=1:nn,
   P = [cp(k) 0 -sp(k);0 1 0;sp(k) 0 cp(k)] ;   % there was a bug here until 14apr12
   R = [1 0 0;0 cr(k) -sr(k);0 sr(k) cr(k)] ;
   H = [ch(k) -sh(k) 0;sh(k) ch(k) 0;0 0 1] ;
   Q(:,:,k) = H*P*R ;
end

Q = squeeze(Q) ;
