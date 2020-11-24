function    [K,s,KK] = zero_crossings(x,TH,Tmax)
%
%    	 [K,s] = zero_crossings(x,TH,Tmax)
%      Find zero-crossings in a vector using a hysteretic detector. This is 
%		 useful, e.g., to locate cyclic postural changes due to propulsion.
%	
%		 Inputs:
%		 x is a vector of data. This can be from any sensor and with any sampling rate.
%    	 TH is the magnitude threshold for detecting a zero-crossing. A zero-crossing
%		  is only detected when values in x pass from -TH to +TH or vice versa.
%      Tmax is the (optional) maximum duration in samples between threshold 
%       crossings. To be accepted as a zero-crossing, the signal must pass from below 
%       -TH to above TH, or vice versa, in no more than Tmax samples. This is useful
%		  to eliminate slow transitions. If Tmax is not given, there is no limit on the
%		  number of samples between threshold crossings.
%  
%      Result:
%		 K is a vector of cues (in samples) to zero-crossings in x.
%		 s is a vector containing the sign of each zero-crossing (1 = positive-going, 
%		  -1 = negative-going). s is the same size as K. If no zero-crossings are found
%		  K and s will be empty
%
%		Example:
%		 [K,s] = zero_crossings(sin(2*pi*0.033*(1:100)'),0.3)
% 	    returns: K=[15.143,30.286,45.429,60.628,75.771,90.914]'
%					 s=[-1,1,-1,1,-1,1]'
%
%     Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     Last modified: 10 May 2017

K = [] ; s = [] ;

if nargin<2,
   help zero_crossings
   return
end

% find all positive and negative threshold crossings
xtp = diff(x>TH) ;
xtn = diff(x<-TH) ;
kpl = find(xtp>0)+1 ;  % leading edges of positive threshold crossings
kpt = find(xtp<0) ;  % trailing edges of positive threshold crossings
knl = find(xtn>0)+1 ;  % leading edges of negative threshold crossings
knt = find(xtn<0) ;  % trailing edges of negative threshold crossings

K = zeros(length(kpl)+length(knl),3) ; % prepare space for the results
cnt = 0 ;
if min(kpl)<min(knl),      % find which direction zero-crossing comes first
   SIGN = 1 ;
else
   SIGN = -1 ;
end

while 1,
   if SIGN==1,
      if isempty(kpl), break, end
      kk = max(find(knt<=kpl(1))) ;
      if ~isempty(kk),
         cnt = cnt+1 ;
         K(cnt,:) = [knt(kk),kpl(1),SIGN] ;
         knt = knt(kk+1:end) ;
         knl = knl(knl>kpl(1)) ;
         kpl = kpl(2:end) ;
      end
      SIGN = -1 ;

   else
      if isempty(knl), break, end
      kk = max(find(kpt<=knl(1))) ;
      if ~isempty(kk),
         cnt = cnt+1 ;
         K(cnt,:) = [kpt(kk),knl(1),SIGN] ;
         kpt = kpt(kk+1:end) ;
         kpl = kpl(kpl>knl(1)) ;
         knl = knl(2:end) ;
      end
      SIGN = 1 ;
   end
end
K = K(1:cnt,:) ;

if nargin==3,
   k = find(K(:,2)-K(:,1)<=Tmax) ;
   K = K(k,:) ;
end

s = K(:,3) ;
KK = K(:,1:2) ;
X = [x(K(:,1)) x(K(:,2))] ;
K = (X(:,2).*K(:,1)-X(:,1).*K(:,2))./(X(:,2)-X(:,1)) ;
