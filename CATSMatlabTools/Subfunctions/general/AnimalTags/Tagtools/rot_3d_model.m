function		rot_3d_model(F,prh,speed)

%		rot_3d_model(F,prh)
%		Rotate a 3d model in a gimbal frame for visualizing animal orientation.
%
%		Input:
%		F is a structure containing plot handles from a previous call to plot_3d_model.
%     prh is a vector or matrix of pitch, roll and heading angles in radians, 
%		 prh = [pitch,roll,heading]. Each row of prh defines a separate orientation and 
%		 these are played in sequence. For a single orientation, prh should be a 3-element
%		 row vector.
%     speed is the optional animation speed, very roughly in Hz (although this depends
%		 on how quickly the figure can be updated by the program). The default value is
%		 20 Hz.
%
%		The 3d model must be first created in the current figure using plot_3dmodel.
%
%		Example:
%		 F=plot_3d_model;
%		 rot_3d_model(F,linspace(0,2*pi,100)'*[0 1 0]) ;
% 	    Plots a dolphin in a gimbal and animates it through a corkscrew roll.
%
%     Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     Last modified: 10 May 2017

if nargin<2,
	help rot_3d_model
	return
end
	
if nargin<3,
	t = 0.05 ;
else
	t = 1/speed ;
end
	
for k=1:size(prh,1),
    Q = euler2rotmat(prh(k,:))' ;
    tpts = F.P*Q ;
    set(F.p,'Vertices',[tpts(:,1) -tpts(:,2) tpts(:,3)]) ;
	 lx = F.LX*Q ;
	 set(F.L(1),'xdata',lx(:,1),'ydata',-lx(:,2),'zdata',lx(:,3));
	 set(F.A(1),'xdata',lx(2,1),'ydata',-lx(2,2),'zdata',lx(2,3));
	 ly = F.LX(:,[3 1 2])*Q ;
	 set(F.L(2),'xdata',ly(:,1),'ydata',-ly(:,2),'zdata',ly(:,3));
	 set(F.A(2),'xdata',ly(2,1),'ydata',-ly(2,2),'zdata',ly(2,3));
	 lz = F.LX(:,[2 3 1])*Q ;
	 set(F.L(3),'xdata',lz(:,1),'ydata',-lz(:,2),'zdata',lz(:,3));
	 set(F.A(3),'xdata',lz(2,1),'ydata',-lz(2,2),'zdata',lz(2,3));
	 if ~isempty(F.C),
		cx = F.CX*Q ;
		set(F.C(1),'xdata',cx(:,1),'ydata',-cx(:,2),'zdata',cx(:,3));
		cy = F.CX(:,[3 1 2])*Q ;
		set(F.C(2),'xdata',cy(:,1),'ydata',-cy(:,2),'zdata',cy(:,3));
		cz = F.CX(:,[2 3 1])*Q ;
		set(F.C(3),'xdata',cz(:,1),'ydata',-cz(:,2),'zdata',cz(:,3));
	 end
    drawnow ;
    pause(t);
end
