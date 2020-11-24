function [Aw,Mw,Gw] = applyWduringslip(At,Mt,Gt,W1,W2)

q1 = dcm2quat(W1);
q2 = dcm2quat(W2);
qd = quatdivide(q2,q1);
ang = 2*acos(qd(1));
n = size(At,1);
% a quaternion is rotation of angle alpha around axes x y z (given a unit magnitude) equiv to
% cos(alpha/2) + sin(alpha/2)(xi + yj + zk).
% See https://www.mathworks.com/help/fusion/ug/rotations-orientation-and-quaternions.html for more.
qdn = [cos(ang/(n-1)/2),sin(ang/(n-1)/2)*qd(2:4)/sin(ang/2)];
Aw = At; Mw = Mt; Gw = Gt;
for j = 1:n;
    if j == 1; qt = q1; else
        qt = quatmultiply(qt,qdn);
    end
    Wt = quat2dcm(qt);
    Aw(j,:) = At(j,:)*Wt;
    Mw(j,:) = Mt(j,:)*Wt;
    Gw(j,:) = Gt(j,:)*Wt;
end
