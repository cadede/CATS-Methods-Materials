function    [TOL,fc]=spec2tol(SL,f)

%    	[TOL,fc]=spec2tol(SL,f)
%	  	Estimate third octave levels from FFT power spectra.
%
%	 	Inputs:
%    	SL is a vector or matrix of power spectra in dB re U^2/Hz
%      where U is any appropriate unit. SL can be produced by
%		 wavSL. If SL is a vector, it is treated as a single spectrum.
%		 If SL is a matrix, each column is treated as a separate spectrum.
%     f is the centre frequency of each row in SL.
%
%	   Returns:
%	   TOL is a matrix of third octave levels in dB re U^2 RMS
%	   fc is a vector with the centre frequencies of the third octaves.
%		 Only the third octaves that can be estimated from SL are
%		 returned. These are determined by the frequency resolution and
%      upper frequency limit of SL.
%
%		Example:
%		 TBD
%
%	   Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     Last modified: 18 Feb 2017

if nargin<2,
	help spec2tol
	return
end
	
if size(SL,1)==1,
	SL = SL(:) ;
end
	
Fc = 1000*((2^(1/3)).^(-16:1:30));        % Exact center freq. 	
f1 = Fc/(2^(1/6)); 
f2 = Fc*(2^(1/6)); 
fres = f(2)-f(1) ;
bw = f2-f1 ;
kf = find(fres<bw & f2<=max(f)) ; 
P = 10.^(SL/10) ;
top = NaN*ones(length(kf),size(P,2)) ;
for k=1:length(kf),
   kk = find(f>=f1(kf(k)) & f<f2(kf(k))) ;
   if length(kk)==1,
      top(k,:) = P(kk,:) ;
   else
      top(k,:) = mean(P(kk,:)) ;
   end
end

TOL = 10*log10(top)+repmat(10*log10(bw(kf)'),1,size(top,2)) ;
fc = Fc(kf) ;
