function		F = plot_3d_model(fname)

%		F = plot_3d_model(fname)
%		Plot a 3d model in a gimbal frame for visualizing animal orientation.
%
%		Input:
%		fname is the optional name of a pair of files containing the points and
%		 connections of a wire frame. These files should have suffix .pts and .knx,
%		 respectively. If fname is not given, files for a 3d plot of a dolphin will
%		 be loaded.
%
%		Returns:
%		F is a structure containing handles to the components of the display which
%		 can then be manipulated using rot_3dmodel.
%
%		The 3d model is plotted in the current figure. This function clears any plot
%		already in the figure to ensure that the plot will appear.
%
%		Example:
%		 F=plot_3d_model;
%		 rot_3d_model(F,linspace(0,2*pi,100)'*[0 1 0]) ;
% 	    Plots a dolphin in a gimbal and animates it through a corkscrew roll.
%
%     Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     Last modified: 10 May 2017

if nargin<1 || isempty(fname)
	fname = 'dolphin' ;
end
	
SZ = 5 ;			% size of the gimbal
CSZ = 1.07 ;	% multiplier for size of the compass wheel
PLOT_RINGS = 0 ;

% load wire frame
warning off
P=load('-ascii',[fname '.pts']);
K=load('-ascii',[fname '.knx']);
warning on
clf

% convert from
P = [-P(:,1) P(:,2:3)] ;
p = patch('faces',K(:,1:3)+1, 'vertices',[P(:,1) -P(:,2) P(:,3)]);
set(p, 'facealpha',1)
colormap(flipud(bone))
brighten(-0.3)
set(p, 'FaceVertexCData',P(:,3));
set(p, 'FaceColor', 'flat');
shading flat
axis square
hold on

c = exp(j*2*pi*(0:500-1)'/500);
CX = SZ*[real(c) imag(c) zeros(500,1)];

if PLOT_RINGS
	C(1)=plot3(CX(:,1),-CX(:,2),CX(:,3),'k');
	C(2)=plot3(CX(:,3),-CX(:,1),CX(:,2),'k');
	C(3)=plot3(CX(:,2),-CX(:,3),CX(:,1),'k');
	set(C,'color',0.6*[1 1 1])
else
	C = [] ;
end
	
COMP=plot3(CSZ*CX(:,1),CSZ*CX(:,2),CX(:,3),'k');
CTICK=plot3(CSZ*SZ*[1 0 -1 0;CSZ 0 -CSZ 0],CSZ*SZ*[0 1 0 -1;0 CSZ 0 -CSZ],zeros(2,4),'k') ;
set(COMP,'LineWidth',1)
set(CTICK,'LineWidth',1)
text(CSZ^3*SZ,0,0,'N','Color','k')
LX = SZ*[-1 0 0;1 0 0] ;
L(1)=plot3(LX(:,1),-LX(:,2),LX(:,3),'r');
L(2)=plot3(LX(:,3),-LX(:,1),LX(:,2),'g');
L(3)=plot3(LX(:,2),-LX(:,3),LX(:,1),'b');
set(L,'LineWidth',1.5)
AX = SZ*[1 0 0] ;
A(1)=plot3(LX(2,1),-LX(2,2),LX(2,3),'r.');
A(2)=plot3(LX(2,3),-LX(2,1),LX(2,2),'g.');
A(3)=plot3(LX(2,2),-LX(2,3),LX(2,1),'b.');
set(A,'MarkerSize',14);
axis([-1 1 -1 1 -1 1]*SZ*sqrt(2))
view([30 20])
axis off
F.p = p ;
F.P = P ;
F.LX = LX ;
F.L = L ;
F.CX = CX ;
F.C = C ;
F.A = A ;
