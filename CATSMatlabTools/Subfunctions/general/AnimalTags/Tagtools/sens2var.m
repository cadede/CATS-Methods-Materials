function    [X,Y,fs] = sens2var(Sx,Sy,r)

%    [X,fs] = sens2var(Sx)     % regularly sampled data
%    or
%    [X,T] = sens2var(Sx)      % irregularly sampled data
%    Extract data from a sensor structure.
%    Can also be called with two variables, in which case, they
%    are checked for compatibility (i.e., same length and sampling):
%    [X,Y,fs] = sens2var(Sx,Sy); % regularly sampled data
%    or
%    [X,Y,T] = sens2var(Sx,Sy);  % irregularly sampled data
%    Can also be called with a trailing string 'regular' to check if
%    the sensor structures are regularly sampled. If not, X will be
%    returned empty.
%
%    Inputs:
%    Sx, Sy must be sensor structures. If not X will be returned empty.
%
%    Returns:
%    X, Y are vectors or matriced of sensor data. These are in the unit and
%     frame as stored in the input sensor structures.
%    fs is the sampling rate of the sensor data in Hz (samples per second).
%    T is the time in seconds of each measurement in data for irregularly
%     sampled data. The time reference (i.e., the 0 time) is with
%     respect to the start time of the data in the sensor structure.
%
%    Example:
%     loadnc('testset3')
%     [pca_dur,pca_time] = sens2var(PCA)
%
%    Valid: Matlab, Octave
%    markjohnson@st-andrews.ac.uk
%    Last modified: 22 July 2017

X = [] ; Y = [] ; fs = [] ;
if nargin<1,
   help sens2var
   return
end

if ~(isstruct(Sx) && isfield(Sx,'data') && isfield(Sx,'sampling')),
   fprintf(' Error: input argument must be a sensor structure\n') ;
   return
end

if nargin==1,
   r = [] ; Sy = [] ;
elseif nargin==2
   if ischar(Sy),
      r = Sy ;
      Sy = [] ;
   else
      r = [] ;
   end
end

if nargin>1 && ~isempty(Sy),
   if ~isstruct(Sy) || ~isfield(Sy,'data') ~isfield(Sy,'sampling')
      fprintf(' Error: input argument must be a sensor structure\n') ;
      return
   end
end

R = [1 1]*strcmp(Sx.sampling,'regular') ;
if ~isempty(Sy),
   R(2) = strcmp(Sy.sampling,'regular') ;
end
if ~isempty(r) && strcmpi(r,'regular') && sum(R)<2,
   fprintf(' Error: input argument must be regularly sampled\n') ;
   return
end 

if sum(R)==1,
   fprintf(' Error: input arguments must both be sampled in the same way\n') ;
   return
end 
   
if R(1),
   X = Sx.data ;
   fs = Sx.sampling_rate ;
else   
   fs = Sx.data(:,1) ;
   if size(Sx.data,2)>1,
      X = Sx.data(:,2:end) ;
   else
      X = ones(length(fs),1) ;
   end
end

if isempty(Sy),
   Y = fs ;
   return
end
   
% here on for two input variables
if R(2),
   if fs ~= Sy.sampling_rate,
      fprintf(' Error: input arguments must both have the same sampling rate\n') ;
      X = [] ; Y = [] ; fs = [] ;
      return
   end
   Y = Sy.data ;
      
else
   if size(Sy.data,2)>1,
      Y = Sy.data(:,2:end) ;
   else
      Y = ones(size(Sy.data,1),1) ;
   end
end

if size(X,1)~=size(Y,1),
   fprintf(' Error: input arguments must both have the same number of samples\n') ;
   X = [] ; Y = [] ; fs = [] ;
   return
end
   
