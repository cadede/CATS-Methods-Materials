function    [Y,G,output] = spherical_calwk(X,n,method,bound)

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
	
if nargin<3 || isempty(method),
	method = [] ;
end

if nargin<4
    bound = false;
end

		
% remove any rows in X with NaNs
k = find(all(~isnan(X),2)) ;
X = X(k,:) ;

nv1 = 3 ;		% number of variables for offset
nv2 = 5 ;		% number of variables for gain and offset
nv3 = 8 ;		% number of variables for gain, offset and cross

% start by estimating offsets using linear least squares. This ensures
% that the iterative search starts fairly close to a solution.
bsq = sum(X.^2,2) ;
XX = [2*X ones(size(X,1),1)];
R = XX'*XX ;
P = sum(repmat(bsq,1,4).*XX) ;
H = -inv(R)*P' ;
offs = H(1:3) ;
X = X+repmat(offs',size(X,1),1) ;

options = optimset('MaxFunEvals',10000,'MaxIter',10000,'TolX',0.00001,'TolFun',0.00001);
% now try up to three calibration scenarios using simplex search

C = zeros(nv3,3) ;
if ~bound
    [C(1:nv1,1),~,~,output] = fminsearch(@(c) ccost(c,X),zeros(nv1,1),options) ;		% offset only cal
%     disp(output);
    if strcmp(method,'gain') | strcmp(method,'cross'),
        [C(1:nv2,2),~,~,output] = fminsearch(@(c) ccost(c,X),C(1:nv2,1),options) ;	% offset and gain cal
    end
%     disp(output);
    if strcmp(method,'cross'),
        [C(:,3),~,~,output] = fminsearch(@(c) ccost(c,X),C(:,2),options) ;		% offset, gain and cross cal
    end
%     disp(output);
else
    LB = [min(min(X))*.1*ones(1,3) 0.1 0.1 -10 -10 -10];
    UB = [max(max(X))*.1*ones(1,3) 10 10 10 10 10];%, options) ;
    [C(1:nv1,1),~,~,output] = fminsearchbnd(@(c) ccost(c,X),zeros(nv1,1),LB(1:nv1),UB(1:nv1),options) ;		% offset only cal
%     disp(output);
    if strcmp(method,'gain') | strcmp(method,'cross'),
        [C(1:nv2,2),~,~,output] = fminsearchbnd(@(c) ccost(c,X),C(1:nv2,1),LB(1:nv2),UB(1:nv2),options) ;	% offset and gain cal
    end
%     disp(output);
    if strcmp(method,'cross'),
        [C(:,3),~,~,output] = fminsearchbnd(@(c) ccost(c,X),C(:,2),LB,UB,options) ;		% offset, gain and cross cal
    end
%     disp(output);
%    [C(1:nv2,2),~,~,output] = fminsearchbnd(@(c) ccost(c,X),zeros(5,1),[min(min(X))*2*ones(1,3) 0.1 0.1],[max(max(X))*2*ones(1,3) 10 10], options) ;
 
end


% if ~forcegain
[m,k] = min(ccost(C,X)) ;  		% pick the best performer
% else
%     [m,k] = min(ccost(C(:,2:end),X)) ;
%     k = k+1;
% end
C = C(:,k) ;
[Y,C] = appcal(X,C) ;		% apply the calibration
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

G.poly = [(1+C(:,2))*sf (offs.*(1+C(:,2))+C(:,1))*sf] ;
G.cross = 0.5*[2 C(1,3) C(3,3);C(1,3) 2 C(2,3);C(3,3) C(2,3) 2] ;
return


function		p = ccost(C,X)
for k=1:size(C,2),
	n = sqrt(sum(appcal(X,C(:,k)).^2,2)) ;
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

