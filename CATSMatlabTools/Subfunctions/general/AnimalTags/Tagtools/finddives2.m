function    T = finddives2(p,fs,th,surface,findall)
% finddives2 counts the first dive (treats findall as findfirst, replace
% first p with 0 to ensure this is caught;


%    T = finddives(p,fs,[th,surface,findall])
%    Find time cues for the edges of dives.
%    p is the depth time series in meters, sampled at fs Hz.
%    th is the threshold in m at which to recognize a dive - dives
%    more shallow than th will be ignored. The default value for th is 10m.
%    surface is the depth in meters at which it is considered that the
%    animal has reached the surface. Default value is 1.
%    findall = 1 forces the algorithm to include incomplete dives at the
%    start and end of the record. Default is 0
%    T is the matrix of cues with columns:
%    [start_cue end_cue max_depth cue_at_max_depth mean_depth mean_compression]
%
%    If there are n dives deeper than th in p, then T will be an nx6 matrix. Partial
%    dives at the beginning or end of the recording will be ignored - only dives that
%    start and end at the surface will appear in T. 
%
%    mark johnson, WHOI
%    mjohnson@whoi.edu
%    last modified: 25 October 2005

if p(1)~=0; p(1) = 0; end

if nargin<2,
   help('finddives') ;
   return
end

if nargin<3 | isempty(th),
   th = 10 ;
end

if nargin<4 | isempty(surface),
   surface = 1 ;        % maximum p value for a surfacing (was 2)
end

if nargin<5,
   findall = 1 ;
end

if fs>1000,
   fprintf('Suspicious fs of %d Hz - check\n', round(fs)) ;
   return
end

searchlen = 20 ;        % how far to look in seconds to find actual surfacing
dpthresh = 0.25 ;        % vertical velocity threshold for surfacing
dp_lp = 0.5 ;           % low-pass filter frequency for vertical velocity

% first remove any NaN at the start of p
% (these are used to mask bad data points and only occur in a few data sets)
kgood = find(~isnan(p)) ;
p = p(kgood) ;
tgood = (min(kgood)-1)/fs ;

% find threshold crossings and surface times
tth = find(diff(p>th)>0) ;
tsurf = find(p<surface) ;
ton = 0*tth ;
toff = ton ;
k = 0 ;

% sort through threshold crossings to find valid dive start and end points
for kth=1:length(tth) ;
   if all(tth(kth)>toff),
      ks0 = find(tsurf<tth(kth)) ;
      ks1 = find(tsurf>tth(kth)) ;
      if (findall & ~isempty(ks1)) | (~isempty(ks0) & ~isempty(ks1)) % modified to only go for the first dive, not interrupted last dives
         k = k+1 ;
         if isempty(ks0),
            ton(k) = 1 ;
         else
            ton(k) = max(tsurf(ks0)) ;
         end
         if isempty(ks1),
            toff(k) = length(p) ;
         else
            toff(k) = min(tsurf(ks1)) ;
         end
      end
   end
end

% truncate dive list to only dives with starts and stops in the record
ton = ton(1:k) ;
toff = toff(1:k) ;

% filter vertical velocity to find actual surfacing moments
try [b a] = butter(4,dp_lp/(fs/2)) ;
catch; [b a] = butter(4,.99) ;
end
dp = filtfilt(b,a,[0;diff(p)]*fs) ;

% for each ton, look back to find last time whale was at the surface
% for each toff, look forward to find next time whale is at the surface
dmax = zeros(length(ton),2) ;
for k=1:length(ton),
   ind = ton(k)+(-round(searchlen*fs):0) ;
   ind = ind(find(ind>0)) ;
   ki = max(find(dp(ind)<dpthresh)) ;
   if isempty(ki),
      ki=1 ;
   end
   ton(k) = ind(ki) ;
   ind = toff(k)+(0:round(searchlen*fs)) ;
   ind = ind(find(ind<=length(p))) ;
   ki = min(find(dp(ind)>-dpthresh)) ;
   if isempty(ki),
      ki=1 ;
   end
   toff(k) = ind(ki) ;
   [dm km] = max(p(ton(k):toff(k))) ;
   dmax(k,:) = [dm (ton(k)+km-1)/fs+tgood] ;
end

% measure dive statistics
pmean = 0*ton ;
pcomp = pmean ;
for k=1:length(ton),
   pdive = p(ton(k):toff(k)) ;
   pmean(k) = mean(pdive) ;
   pcomp(k) = mean((1+0.1*pdive).^(-1)) ;
end
 
% assemble output
T = [[ton toff]/fs+tgood dmax pmean pcomp] ;

