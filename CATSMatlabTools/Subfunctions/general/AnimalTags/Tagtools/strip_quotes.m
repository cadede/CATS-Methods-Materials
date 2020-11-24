function    s = strip_quotes(s)

%     s = strip_quotes(s)
%     Remove bracketing double quotes from string or
%     cell array containing strings.
%
%     Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     last modified: 12 July 2017

if nargin<1,
   help strip_quotes
end

if nargin<2,
   nowhite = 0 ;
end

if ischar(s),
   k = find(isspace(s)==0) ;
   s = s(k(1):k(end)) ;
   if length(s)>=2 && all(s([1 end])=='"'),
      s = s(2:end-1) ;
   end
   return
end

[m,n] = size(s) ;
for k=1:m*n,
   ss = s{k} ;
   if ~isempty(ss) && ischar(ss),
      kk = find(isspace(ss)==0) ;
      ss = ss(kk(1):kk(end)) ;
      if length(ss)>=2 && all(ss([1 end])=='"'),
         ss = ss(2:end-1) ;
      end
      s{k} = ss ;
   end
end
return
