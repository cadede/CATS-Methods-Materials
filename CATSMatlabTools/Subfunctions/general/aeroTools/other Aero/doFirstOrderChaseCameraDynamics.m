function doFirstOrderChaseCameraDynamics(~,Bodies,h)
%DOFIRSTORDERCHASECAMERADYNAMICS follows the first body as a "chase plane"
% This is the default 'PositionFcn' for Aero.Camera

   % Copyright 2018 The MathWorks, Inc. 

%--- Location of camera is a "chase plane",
%    use first body as target position reference

if ~isempty(Bodies) && isa(Bodies, 'cell') && isa(Bodies{1},'Aero.Body')
    target     = Bodies{1};
    targetPos  = target.Position;
else
    targetPos = h.AimPoint; % don't change anything
end

pos        = h.Offset;
viewExtent = h.ViewExtent;

%--- Extent of view to render

h.xlim = targetPos(1) + viewExtent;
h.ylim = targetPos(2) + viewExtent;
h.zlim = targetPos(3) + viewExtent;

%--- Camera aim point dynamics for [x,y,z] aim point

tauFactPoint = 0.90;
prevPoint = h.AimPoint;
camtgt = prevPoint + tauFactPoint * (targetPos - prevPoint);
h.AimPoint = camtgt;

%--- Camera position dynamics for [x,y,z] position

tauFactPos = 0.99;
prevPos = h.Position;
campos = prevPos + tauFactPos * (targetPos + pos - prevPos);
h.Position = campos;

end

%[EOF] doFirstOrderChaseCameraDynamics.m