function    [SL,f] = spectrum_level(x,nfft,fs,w,nov)

%    	[SL,f] = spectrum_level(x,nfft,fs)
%		or
%    	[SL,f] = spectrum_level(x,nfft,fs,w)
%		or
%    	[SL,f] = spectrum_level(x,nfft,fs,w,nov)
%
%     Spectrum level of a signal x, i.e., the amount of power per 1Hz
%		band. This replaces Matlab's deprecated psd function. The input
%		signal is divided into overlapping pieces equal in length to the
%		required Fast Fourier Transform (FFT) length. Each piece is
%		windowed and the FFT computed. The spectral power is then estimated
%		from the mean of the spectral magnitudes squared. Power is scaled
%		to account for the scale factor of the FFT and the window. The power
%		is also scaled by 10log10 of the bin width in Hz (i.e., the sampling
%		rate divided by the FFT length) to convert the per-bin powers into
%		approximate per-Hz powers. This scaling method is suitable for
%		wideband signals (i.e., with bandwidth wider than the bin width) but
%		is NOT suitable for narrow band and tonal signals.
%
%		Inputs:
%     x is a vector or matrix containing the signal(s) to be processed. 
%		 For signals with multiple channels, each channel should be in a 
%		 column of x.
%     nfft is the length of the FFT to use. Choose a power of two for
%      fastest operation.
%     fs is the sampling rate of the signals in x in Hz.
%     w is the optional window length. The default value is nfft. If w<nfft,
%      each segment of w samples is zero-padded to nfft. Use w=[] to use
%		 the default value if you want to specify nov.
%     nov is the number of samples to overlap each segment. The default
%      value is half of the window length.
%
%     Returns:
%     SL is the spectrum level at each frequency in dB RMS re root-Hz.
%      The spectrum is single-sided and extends to fs/2.
%      The reference level is 1.0 (i.e., white noise with unit variance
%      will have a spectrum level of 3-10*log10(fs). The 3dB is because
%      both the negative and positive spectra are added together so that
%      the total power in the signal is the same as the total power in 
%      the spectrum.
%     f is the vector of frequencies in Hz at which SL is calculated.
%
%		Example:
%		 TBD
%
%     Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     Last modified: 23 July 2017

if nargin<3
	help spectrum_level
	return
end

if nargin<4 || isempty(w),
   w = nfft ;
end

if nargin<5,
   if length(w)==1,
      nov = round(w/2) ;
   else
      nov = round(length(w)/2) ;
   end
end

if length(w)==1,
   w = hanning(w) ;
end

P = zeros(nfft/2,size(x,2)) ;
for k=1:size(x,2),
   [X,z] = buffer(x(:,k),length(w),nov,'nodelay') ;
   X = detrend(X).*repmat(w,1,size(X,2)) ;
   F = abs(fft(X,nfft)).^2 ;
   P(:,k) = sum(F(1:nfft/2,:),2) ;
end

ndt = size(X,2) ;

% these two lines give correct output for randn input
% SL of randn should be -10*log10(fs/2)

slc = 3-10*log10(fs/nfft)-10*log10(sum(w.^2)/nfft) ;

% 3 is to go from a double-sided spectrum to a single-sided (positive frequency) spectrum.
% fs/nfft is to go from power per bin to power per Hz
% sum(w.^2)/nfft corrects for the window

SL = 10*log10(P)-10*log10(ndt)-20*log10(nfft)+slc ;

% 10*log10(ndt) corrects for the number of spectra summed in P (i.e., turns the sum into a mean)
% 20*log10(nfft) corrects the nfft scaling in matlab's fft

f = (0:nfft/2-1)/nfft*fs ;
