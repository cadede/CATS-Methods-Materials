function    [p,pc] = fix_pressure(p,t,fs)

%     [p,pc] = fix_pressure(p,t,maxp)		% p and t are sensor structures
%		or
%     [p,pc] = fix_pressure(p,t,fs,maxp)		% p and t are vectors
%
%     Correct a depth or altitude profile for offsets caused by
%     mis-calibration and temperature. This function finds minima
%     in the dive/altitude profile that are consistent with
%     surfacing/landing. It uses the depth/height at these points
%     to fit a temperature regression.
%
%     Inputs:
%     p is a sensor structure or vector of depth/altitude in meters.
%     t is a sensor structure or vector of temperature in degrees Celsius.
%     fs is the sampling rate of p and t in Hz. This is only needed if
%		 p and t are not sensor strucures. The depth and temperature
%      must both have the same sampling rate (use decdc.m or resample.m 
%		 if needed to achieve this).
%		maxp is the maximum depth or altitude reading in the pressure data
%		 for which the animal could actually be at the surface. This is a
%		 rough measurement of the potential error in the pressure data. The
%		 unit is meters. Start with a small value, e.g., 2 m and re-run fix_depth
%		 with a larger value if there are still obvious temperature-related
%		 errors in the resulting depth/altitude profile.
%
%     Results:
%     p is a sensor structure or vector of corrected depth/altitude 
%		 measurements at the same sampling rate as the input data. If the
%		 input is a sensor structure, the output will also be. 
%     pc is a structure containing the pressure offset and temperature 
%		 correction coefficients. It has fields:
%		 pc.tref is the temperature reference in degrees Celsius (this will
%		  always be 20.
%		 pc.tcomp is the temperature compensation polynomial. This is used
%		  within the function to correct pressure as follows:
%       p = p+polyval(pc.tcomp,t-pc.tref) ;
%
%		This function makes a number of assumptions about the depth/altitude
%     data and about the behaviour of animals:
%     - the depth data should have few incorrect outlier (negative) values
%       that fall well beyond the surface. These can be reduced using
%       median_filter.m before calling fix_depth.
%     - the animal is assumed to be near the surface at least 2% of the
%       time. If the animal is less frequently at the surface, you may need
%       to change the value of PRCTSURF near the start of the function.
%     - potential surfacings are detected by looking for zero-crossings in
%       the vertical speed and this requires defining a threshold in vertical
%       speed that must be crossed by each zero crossing. The value used is
%       0.05 m/s but this may be too high for animals that move very slowly
%       near the surface. In which case, change MAXSPEED near the start of
%       the function.
%
%		Example:
%		 loadnc('mn12_186a_raw')
%		 [PP,pc] = fix_depth(P,T,10);
%		 plott(P,PP)
%		 % lower plot shows the compensated pressure which is closer to 0 when
%		 % the animal is at the surface
%
%     Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     Last modified: 3 June 2017

MAXSPEED = 0.05 ;    % maximum speed in metres/second of points at the surface to accept
ASYMM = 0.2 ;        % maximum assymmetry between positive and negative residuals
TREF = 20 ;          % standard temperature reference to use
PRCTSURF = 2 ;       % minimum percent of time animal is near surface

pc = [] ;
if nargin<2,
   help fix_pressure ;
	return
end

if isstruct(p),
	[pp,tt,fs] = sens2var(p,t,'regular') ;
	if isempty(pp), return, end
	if isfield(p,'cal_tcomp'),		% if there is already some temperature compensation, undo it.
		pp = pp - polyval([p.cal_tcomp,0],tt - p.cal_tref) ;
	end
	if nargin>3,
		maxp = fs ;
	else
		maxp = [] ;
	end
else
	if nargin<3,
		fprintf(' Sampling rate is a required input when p and t are not sensor structures\n') ;
		return
	end
	pp = p(:) ;
	tt = t(:) ;
	if length(p)~=length(t),
		fprintf(' Error: p and t must have the same length\n') ;
		return
	end
	if nargin<4,
		maxp = [] ;
	end
end

if isempty(maxp),
	maxp = 5 ;
end
	

if fs>5,
   df = round(fs/5) ;
   pp = decdc(pp,df) ;   % decimate depth and temperature to around 5Hz
   tt = decdc(tt,df) ;
   fs = fs/df ;
end

v = depth_rate(pp,fs) ;     % compute vertical velocity with 5s timeconstant
if isempty(v),
	return
end

% do initial offset correction - just using the 2%ile of the depth. This
% assumes animals spend at least 2% of their time at or close to the
% surface. Lower the percentile if this is not the case.
p0 = -prctile(pp,PRCTSURF) ;
pp = pp+p0 ;
k = find(pp<maxp) ;
pp = pp(k) ;
tt = tt(k) ;
v = v(k) ;

[K,s,KK] = zero_crossings(v,MAXSPEED) ; % find zero crossings of vertical velocity
KK=KK(s>0,:) ;     % pick just the positive zero crossings
% these are when the animal goes from descending to ascending if flying
% or from ascending to descending is swimming

% select depth samples around each zero crossing
k=zeros(length(v),1);
for kk=1:size(KK,1),
	k(KK(kk,1):KK(kk,2))=1 ;
end
k=find(k);
ps=pp(k);      % pick just the 'surface' samples of pressure
ts=tt(k)-TREF;      % and temperature

% do several iterations of regression followed by removal of the largest
% positive residuals - these are non-surface points that have survived the
% previous data selection steps. This approach relies on there being a
% relatively small proportion of non-surface samples by this point,
% certainly less than 50%. If this is not true, then the previous data
% selection must be improved.
for k=1:10,    % put an upper limit on the number of iterations
   [b,bi,r,ri,stats]=regress(ps,[ts.^2 ts ones(length(ps),1)]);
   rr = sqrt([sum(r(r>0).^2) sum(r(r<0).^2)]) ;  % RMS of +ve and -ve residuals
   if abs(diff(rr))/mean(rr)<ASYMM, break, end  % if ratio of +ve and -ve is similar, break
   kk=find(r<prctile(r,90));  % otherwise eliminate the top 10% depth samples
   ps=ps(kk);
   ts=ts(kk);
end

% compute the correction terms
pc.tref = TREF ;
if stats(3)>0.05,			% if the regression didn't help, just keep the offset adjustment
	pc.tcomp = [0 0] ;
	pc.poly = [1 p0] ;
else
	pc.tcomp = -b(1:2)' ;
	pc.poly = [1 p0-b(3)] ;
end

% correct the original pressure
if ~isstruct(p),
	p = p+polyval([pc.tcomp,pc.poly(2)],t-TREF) ;
	return
end

p.data = p.data + polyval([pc.tcomp,pc.poly(2)],t.data - TREF) ;
p.cal_tref = TREF ;
p.cal_tcomp = pc.tcomp ;

if isfield(p,'cal_poly'),
	p.cal_poly = p.cal_poly*pc.poly(1) + [0 pc.poly(2)] ;
else
	p.cal_poly = pc.poly ;
end

if ~isfield(p,'history') || isempty(p.history),
	p.history = 'fix_depth' ;
else
	p.history = [p.history ',fix_depth'] ;
end
