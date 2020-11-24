function    [ax,h]=plott(varargin)

%     [ax,h]=plott(X)			% X is a sensor structure
%	   or
%     [ax,h]=plott(X,r)			% X is a sensor structure
%	   or
%     [ax,h]=plott(X,fsx)		% X is a vector or matrix of sensor data
%	   or
%     [ax,h]=plott(X,fsx,r)	% X is a vector or matrix of sensor data
%	   or
%     [ax,h]=plott(X,Y,...)	% X, Y etc are sensor structures
%	   or
%     [ax,h]=plott(X,fsx,Y,fsy,...)	% X, Y etc are vectors or matrices of sensor data
%     Plot sensor time series against time in a single or multi-paneled figure with linked
%	   x-axes. This is useful for comparing measurements across different sensors. The
%		time axis is automatically displayed in seconds, minutes, hours, or days according
%		to the span of the data.
%
%	   Inputs:
%		X, Y, etc, are sensor structures or vectors/matrices of time series data.
%		fsx, fsy, etc, are the sampling rates in Hz for each data object. Sampling
%		 rates are not needed when the data object is a sensor structure. 
%		r is an optional argument which has several uses:
%      If r='r', the direction of the y-axis is reversed for the data
%      object being plotted. This is useful for plotting dive profiles which 
%      match the physical situation i.e., with greater depths lower in the
%	    display. Note that 'r' is the default for a sensor structure if the
%      axes field has a value of 'D' for down.
%      If r='i', the preceding data is taken as irregularly sampled. The
%      data must have at least 2 columns with the first one being the
%      sampling times. In this case, data are plotted as single points
%      rather than a continuous line.
%      If r is a number, it specifies the number of seconds time offset 
%		 for the preceding data object. A positive value means that these 
%      data were collected later than the other objects and so should be 
%      plotted more to the right-hand side. Note that the correct time
%      offset is automatically used if the input is a sensor structure with
%      a start_offset field.
%
%		Returns:
%   	ax is a vector of handles to the axes created.
%		h is a cell array of vectors of handles to the lines plotted. There is a cell of
%		 handles for each axis.
%
%		This is a flexible plotting tool which can be used to display and explore sensor
%		data with different sampling rates on a uniform time grid. Zooming any of the
%	   panels should cause all of the panels to zoom in or out to match the x-axis.
%
%		Example:
%		 loadnc('testset3');
%		 plott(P,A,PCA)				% plot depth, acceleration and prey-capture-attempts
%
%     Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     Last modified: 02 August 2017 by RJS

ax=[]; h=[];
if nargin<1,		% must have at least one input argument
   help plott
   return
end

brk = [0,2e3,2e4,5e5] ;		% break points for plots in seconds, mins, hours, days
div = [1 60 3600 24*3600] ;	% corresponding time multipliers
L = {'s','min','hr','day'} ;	% and xlabels

% each data object can have one or two qualifying arguments. Scan through varargin
% to find the objects and their qualifiers.
X = {} ;
T = {} ;
% fsrt is a variable to collect the sampling rate, reverse flag and time
% offset of each panel
fsrt = zeros(length(varargin),3) ;
for k=1:nargin,
	x = varargin{k} ;
	if isstruct(x),		% this input is a sensor structure
		if isfield(x,'sampling') && isfield(x,'data'),
         if strcmp(x.sampling,'regular'),
            X{end+1} = x.data ;
            T{length(X)} = [] ;
   			fsrt(length(X),1) = x.sampling_rate ;
         else
            if size(x.data,2)>1,
               X{end+1} = x.data(:,2:end) ;
            else
               X{end+1} = ones(length(x.data),1) ;
            end
            T{length(X)} = x.data(:,1) ;
            fsrt(length(X),1) = -1 ;
         end
         if isfield(x,'start_offset'),
            fsrt(length(X),3) = x.start_offset ;
         end
         if isfield(x,'axes') && length(x.axes)==1,
            fsrt(length(X),2) = upper(x.axes)=='D' ;
         end
		else
			fprintf('Error: sensor structure must have data and sampling fields\n') ;
			return
		end
	elseif size(x,1)>1 || size(x,2)>1,	% this input is a vector or matrix
		X{end+1} = x ;
		T{length(X)} = [] ; 
      
	else		% this input is a qualifier
		if ischar(x),
         switch x(1),
            case 'r'
               fsrt(length(X),2) = 1 ;
            case 'i'
               fsrt(length(X),1) = -1 ;
               T{length(X)} = X{end}(:,1) ;
               if size(X{end},2)>1,
                  X{end} = X{end}(:,2:end) ;
               else
                  X{end} = ones(length(X{end}),1) ;
               end
            otherwise
               fprintf('Unknown option to plott "%c", skipping\n',x(1)) ;
         end
		else
			if fsrt(length(X),1)==0,
				fsrt(length(X),1) = x ;
			else
				fsrt(length(X),3) = x ;
			end
		end
	end
end

fsrt = fsrt(1:length(X),:) ;
if any(fsrt(:,1)==0),
	fprintf('Error: sampling rate undefined for data object %d\n',find(fsrt(:,1)==0,1)) ;
	return
end
	
ax = zeros(length(X),1) ;
nst = 0 ;
ned = 0 ;
for k=1:length(X),
   ax(k) = subplot(length(X),1,k) ;
   if fsrt(k,1)>0,
      nst = min(nst,fsrt(k,3)) ;
   	ned = max(ned,size(X{k},1)/fsrt(k,1)+fsrt(k,3)) ;
   else
      nst = min(nst,fsrt(k,3)+min(T{k})) ;
      ned = max(ned,max(T{k})+fsrt(k,3)) ;
   end
end

spann = ned-nst ;
for divk=length(brk):-1:1,
   if spann>=brk(divk), break, end
end

ddiv = div(divk) ;
xlims = [nst ned]/ddiv ;
h = {} ;
for k=1:length(X),		% now we are ready to plot
	axes(ax(k)) ;
   if fsrt(k,1)>0,
   	h{k}=plot(((1:size(X{k},1))/fsrt(k,1)+fsrt(k,3))*(1/ddiv),X{k}); grid
   else
   	h{k}=plot((T{k}+fsrt(k,3))*(1/ddiv),X{k},'.'); grid
   end
	set(ax(k),'XLim',xlims)
	if fsrt(k,2)==1,
		set(ax(k),'ydir','reverse') ;
	end
end

xlab = sprintf('Time (%s)',L{divk}) ;
xlabel(ax(end),xlab) ;

if (numel (ax) > 1)==1 %RJS updated 2017-08-02
  linkaxes(ax,'x')
end

if nargout<1 %RJS updated 2017-08-02
   clear ax
end 
