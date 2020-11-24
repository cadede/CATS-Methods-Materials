function    [Y,G] = spherical_cal_T(X,n,T,method)

%     [Y,G] = spherical_cal(X)
%		or
%     [Y,G] = spherical_cal(X,n)
%		or
%     [Y,G] = spherical_cal(X,n,method)
%		Deduce the calibration constants for a triaxial field sensor,
%		such as an accelerometer or magnetometer, based on movement data.
%		This can be used to do a 'bench' calibration of a sensor.
%
%		Inputs:
%     X is the segment of triaxial sensor data to calibrate. It must
%		 be a 3-column matrix. X can come from any triaxial field sensor
%		 and can be in any unit and any frame.
%     n is the target field magnitude e.g., 9.81 for accelerometer 
%      data using m/s2 as the unit.
%     method is an optional string selecting the type of calibration.
%		 The default is to calibrate for offset and scaling only. Other
%		 options are:
%      'gain' adjust gain of axes 2 and 3 relative to 1.
%      'cross' adjust gain and remove cross-axis correlations
%
%		Results:
%     Y is the matrix of converted sensor values. These will have the same
%		 units as for input argument n. The size of Y is the same as the size
%		 of X and it has the same frame and sampling rate.
%     G is the calibration structure containing fields:
%		G.poly is a matrix of polynomials. The first column of G.poly is the 
%		 three scale factors applied to the columns of X. The second column 
%		 is the offset added to each column of X after scaling.
%		G.cross is a 3x3 matrix of cross-factors. If there are no cross-terms, this
%		 is the identity matrix. Off-axis terms correct for cross-axis sensitivity.
%
%		The function reports the residual and the axial balance of the data.
%		A low residual e.g., <5% indicates that the data can be calibrated
%		well and there is not much noise. The axial balance indicates whether
%		the movement in X is suitable for data-driven calibration. If the
%		movement covers all directions fairly equally, the axial balance will
%		be high. A balance <20% may lead to unreliable calibration. For bench
%		calibrations, a high axial balance is achieved by rotating the sensor
%		through the full 3-dimensions.
%		Sampling rate and frame of Y are the same as the input data so Y
%		has the same size as X. The units of Y are the same as the units used for n. 
%		If n is not specified, the units of Y are the same as for the input data.
%		It is a good idea to low-pass filter and/or remove outliers from
%		the sensor data before using this function to reduce errors from 
%		specific acceleration and sensor noise. 
%		Notes: this function uses a Simplex search for optimal calibration parameters
%		and so can be slow if the data size is large. For this reason it is most
%		suitable for bench calibrations rather than field data.
%		This function is only usable for field sensors. It will not work
%		for gyroscope data.

%		Example:
%		 C = spherical_cal()
% 	    returns: C=?.
%
%     Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     Last modified: 10 May 2017

G = []; Y=[];
if nargin<1,
	help spherical_cal
	return
end

if nargin<2,
	n = [] ;
end

if nargin<3 || isempty(T);
    T = zeros(size(X,1),1);
end


if nargin<4,
	method = [] ;
end
		
% remove any rows in X with NaNs
k = find(all(~isnan(X),2)) ;
X = X(k,:) ;

nv1 = 3 ;		% number of variables for offset
nv2 = 5 ;		% number of variables for gain and offset
nv3 = 8 ;		% number of variables for gain, offset and cross
nv4 = 14;       % number of variables for all + Temp cal

% start by estimating offsets using linear least squares. This ensures
% that the iterative search starts fairly close to a solution.
bsq = sum(X.^2,2) ;
XX = [2*X ones(size(X,1),1)];
R = XX'*XX ;
P = sum(repmat(bsq,1,4).*XX) ;
H = -inv(R)*P' ;
offs = H(1:3) ;
X = X+repmat(offs',size(X,1),1) ;

% now try up to three calibration scenarios using simplex search
C = zeros(nv4,4) ;
C(1:nv1,1) = fminsearch(@(c) ccostT(c,X,T),zeros(nv1,1)) ;		% offset only cal

if strcmp(method,'gain') || strcmp(method,'cross'),
	C(1:nv2,2) = fminsearch(@(c) ccostT(c,X,T),C(1:nv2,1)) ;	% offset and gain cal
end

if strcmp(method,'cross'),
	C(1:nv3,3) = fminsearch(@(c) ccostT(c,X,T),C(1:nv3,2)) ;		% offset, gain and cross cal
end

 options = optimset('MaxFunEvals',10^4,'MaxIter',10^4);
[C(1:nv4,4),~,~,output] = fminsearch(@(c) ccostT(c,X,T),C(:,3),options) ;		% offset, gain and cross cal


[m,k] = min(ccostT(C,X,T)) ;  		% pick the best performer
C = C(:,k) ;
[Y,C] = appcalT(X,C,T) ;		% apply the calibration
nn = norm2(Y) ;
fprintf('Residual: %2.1f%%\n',100*std(nn)/mean(nn)) ;
R = Y'*Y ;
fprintf('Axial balance: %2.1f%%\n',100/cond(R)) ;

if ~isempty(n),
   sf = n/mean(nn) ;
   Y = Y*sf ;
else
	sf = 1 ;
end

G.off0 = offs;
% G.poly = [(1+C(:,2))*sf (offs.*(1+C(:,2))+C(:,1))*sf] ;
G.poly = [(1+C(:,2))*sf C(:,1)*sf];
G.Temp = [C(:,5) C(:,4)]*sf;
G.cross = 0.5*[2 C(1,3) C(3,3);C(1,3) 2 C(2,3);C(3,3) C(2,3) 2] ;
return


function		p = ccost(C,X)
for k=1:size(C,2),
	n = sqrt(sum(appcal(X,C(:,k)).^2,2)) ;
	p(k) = std(n)/mean(n) ;
end
return

function		p = ccostT(C,X,T)
for k=1:size(C,2),
	n = sqrt(sum(appcalT(X,C(:,k),T).^2,2)) ;
	p(k) = std(n)/mean(n) ;
end
return

function    [Y,C] = appcal(X,C)
% C is a vector of up to 8 parameters
% Only the first of these may be provided - the remainder are 0.
C(length(C)+1:8) = 0 ;
C = [C(1:3);0;C(4:8)] ;		% add the col1 fixed gain of 0
C = reshape(C,3,[]) ;
%	At this point:
%	C(:,1) are the offsets for each column of X
%	C(:,2) are the gain adjustments for each column of X (column 1 is always 0)
%	C(:,3) are the cross terms
Y = X*diag(1+C(:,2))+repmat(C(:,1)',size(X,1),1) ;
xcm = 0.5*[2 C(1,3) C(3,3);C(1,3) 2 C(2,3);C(3,3) C(2,3) 2] ;
Y = Y*xcm ;
return

function    [Y,C] = appcalT(X,C,T)
% C is a vector of up to 14 parameters
% Only the first of these may be provided - the remainder are 0.
C(length(C)+1:14) = 0 ;
C = [C(1:3);0;C(4:14)] ;		% add the col1 fixed gain of 0
C = reshape(C,3,[]) ;
%	At this point:
%	C(:,1) are the offsets for each column of X
%	C(:,2) are the gain adjustments for each column of X (column 1 is always 0)
%	C(:,3) are the cross terms
%   C(:,4) are the Temp. offsets
%   C(:,5) are the Temp. gains
gains = 1+repmat(C(:,2)',size(X,1),1);
offs = repmat(C(:,1)',size(X,1),1);
gains = gains + repmat(C(1,5)*T+C(2,5)*T.^2+C(3,5)*T.^3,1,3);
offs = offs + repmat(C(1,4)*T+C(2,4)*T.^2+C(3,4)*T.^3,1,3);
% Y = X*diag(1+C(:,2))+repmat(C(:,1)',size(X,1),1) ;
Y = X.*gains + offs;
xcm = 0.5*[2 C(1,3) C(3,3);C(1,3) 2 C(2,3);C(3,3) C(2,3) 2] ;
Y = Y*xcm ;
return


