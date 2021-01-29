function [speed,speedstats] = applySpeed(JigRMS,label1,flownoise,label2,tagon,p,pitch,roll,fs,speedstatsfilepath)
% JigRMS,flownoise and speedstats are required (flownoise or JigRMS can be
% all nans)
% the third input can either be the speedstats variable or a file path from
% which to load speedstats variable (i.e. from another tag)

try load(speedstatsfilepath,'speedstats'); catch; speedstats = speedstatsfilepath; end
if isempty(flownoise) || length(flownoise) == sum(isnan(flownoise)); numRMS = 1; else numRMS = 2; end

try % if JigRMS is a table, else if it's just a vector of RMS data (like flownoise), use it like it is.
JX = JigRMS.X;
JY = JigRMS.Y;
JZ = JigRMS.Z;
J = JigRMS.Mag;
catch
    J = JigRMS;
end

if speedstats.sections_end_index(end) == find(tagon,1,'last');
    disp(['Applying ' num2str(length(speedstats.sections_end_index)) ' speed calibration periods']);
    speedper = [[find(tagon,1); speedstats.sections_end_index(1:end-1)] speedstats.sections_end_index];
    plotnewspeed = false;
else
    disp('Applying whole calibration from speedstats to this deployment in a single speedperiod (from speedperiod 1 of the example file)');
    speedper = [find(tagon,1) find(tagon,1,'last')];
    speedstats.sections_end_index = find(tagon,1,'last');
    plotnewspeed = true;
end

dd = runmean(p,round(speedstats.Thresh.filterSize/2*fs)); %smooth depth to the same filterSize size as your RMS data
dd = [diff(dd); 0];
spitch = runcirc_mean(pitch,round(speedstats.Thresh.filterSize*2*fs)); %smoothed pitch
spitch = [circ_mean([spitch(1:end-1) spitch(2:end)],[],2); nan]; % average each point with the one after it (since each depth deviation is happening between two pitches)
speedSP = -dd./sin(spitch)*fs;
minPitch =speedstats.Thresh.minPitch; if length(minPitch) == 1; minPitch = minPitch*ones(size(speedper,1),1); end
minDepth =speedstats.Thresh.minDepth; if length(minDepth) == 1; minDepth = minDepth*ones(size(speedper,1),1); end
try maxDepth =speedstats.Thresh.maxDepth; catch; maxDepth = max(p); end
    if length(maxDepth) == 1; maxDepth = maxDepth*ones(size(speedper,1),1); end
maxRR =speedstats.Thresh.maxRollRate; if length(maxRR) == 1; maxRR = maxRR*ones(size(speedper,1),1); end
rollrate = abs(runcirc_mean(wrapToPi([diff(roll); 0]),round(fs/2))*180/pi)*fs;


fun = @(a,X)(a(1)*exp(a(2)*(a(3)*X(:,1)+a(4)*X(:,2)+a(5)*X(:,3))));%*(1-(a(3)<0|a(3)>1))*(1-(a(4)<0|a(4)>1))*(1-(a(5)<0|a(5)>1)))*(sum(a(3:5))>.9&sum(a(3:5))<1.1);
speedJJ = nan(size(tagon)); JRMS = speedJJ;
% apply other statistics 
P95 = nan(length(tagon),2);
P68 = P95;
C95 = P95;
P952 = P95; P682 = P68; C952 = C95;
R2 = nan(length(tagon),2);
section = nan(size(tagon));
speed = table(p);


for i = 1:length(speedper(:,1))
    I = speedper(i,1):speedper(i,2);
    oi = speedSP(I);
    oi(abs(spitch(I)*180/pi)<minPitch(i) | p(I)<minDepth(i) | p(I)>maxDepth(i) | rollrate(I) > maxRR(i)) = nan; 
    speedSP(I) = oi;
    try
        m = speedstats.multiModels{i}.Estimate;
        JRMS(I) = log(fun([1 1 m(3:end)'],[JX(I) JY(I) JZ(I)])); %use the weighted contributions of each axis
        %         JRMS(I) = m(1)*exp(m(2)*(m(3)*JX(I) + m(4)*JY(I) + m(5)*JZ(I))); % speed using contributions from each axis
    catch
        disp(['No multi-axis model in section ' num2str(i) ]);
        JRMS(I) = J(I);
    end
    m = [speedstats.Models{i,1}.a speedstats.Models{i,1}.b];
    speedJJ(I) = m(1)*exp(m(2)*JRMS(I));
    P95(I,:) = predint(speedstats.Models{i,1},JRMS(I,1),0.95,'observation','off');
    P68(I,:) = predint(speedstats.Models{i,1},JRMS(I,1),0.68,'observation','off');
    C = confint(speedstats.Models{i,1});
    C95(I,:) = [C(1,1)*exp(C(1,2)*JRMS(I,1)) C(2,1)*exp(C(2,2)*JRMS(I,1))];
    section(I) = i*ones(size(section(I)));
    speed.sectionUsed(I,1) = cellstr(repmat(speedstats.r2used.sectionUsed{i},size(size(section(I)))));
    R2(I,1) = speedstats.r2used.([label1 'r2'])(i)*ones(size(section(I)));
end
    
sectionUsed = speed.sectionUsed;

speedFN = nan(size(speedJJ)); FNsection = zeros(size(speedFN));
useFN = false;
if numRMS == 2
    try speedper2 = [[find(tagon,1); speedstats.(label2).sections_end_index(1:end-1)] speedstats.(label2).sections_end_index];
    catch; speedper2 = speedper;
    end
    for i = 1:length(speedper2(:,1))
        I = speedper2(i,1):speedper2(i,2);
        m = [speedstats.(label2).Models{i,1}.a speedstats.(label2).Models{i,1}.b];
        speedFN(I) = m(1)*exp(m(2)*flownoise(I));
        P952(I,:) = predint(speedstats.(label2).Models{i,1},flownoise(I),0.95,'observation','off');
        P682(I,:) = predint(speedstats.(label2).Models{i,1},flownoise(I),0.68,'observation','off');
        C2 = confint(speedstats.(label2).Models{i,1});
        C952(I,:) = [C2(1,1)*exp(C2(1,2)*flownoise(I)) C2(2,1)*exp(C2(2,2)*flownoise(I))];
        FNsection(I) = i*ones(size(section(I)));
        try ssFN = speedstats.(label2).r2used; useFN = true; catch; ssFN = speedstats.r2used; end
        speed.FNsectionUsed(I,1) = cellstr(repmat(ssFN.sectionUsed{i},size(size(section(I)))));
        R2(I,2) = ssFN.([label2 'r2'])(i)*ones(size(section(I)));
    end
    FNsectionUsed = speed.FNsectionUsed;
end
        
if numRMS == 2
    speed = table(speedSP, speedJJ, JRMS(:,1),P68, P95, C95, R2(:,1), speedFN,flownoise,P682, P952, C952, R2(:,2),section,sectionUsed,...
        'VariableNames',{'SP',label1,'JRMS',[label1 'P68'],[label1 'P95'],[label1 '95'],[label1 'r2'],label2,[label2 'RMS'],[label2 'P68'],[label2 'P95'],[label2 'C95'],[label2 'r2'],'section','sectionUsed'});
else
    speed = table(speedSP, speedJJ(:,1), JRMS(:,1), P68, P95, C95, R2(:,1),nan(size(speedSP)),section,sectionUsed, ...
        'VariableNames',{'SP',label1,[label1 'RMS'],[label1 'P68'],[label1 'P95'],[label1 '95'],[label1 'r2'],'FN','section','sectionUsed'});
end

if useFN; speed.FNsection = FNsection; speed.FNsectionUsed = FNsectionUsed; end % only if there was a separate FN section used

if plotnewspeed
    speedstats.Models = speedstats.Models(1);
    speedstats.ModelFits = speedstats.ModelFits(1,:);
    speedstats.r2used = speedstats.r2used(1,:); for ii = 1:size(speedstats.r2used,2)-1; speedstats.r2used{1,ii} = nan; end
    speedstats.info = ['Used speed calibration from ' speedstatsfilepath];
    speedstats.multiModels = speedstats.multiModels{1};
    I = 1:length(p);
    figure(30); clf; ax = plotyy(I,p,I,speed.JJ);
    set(ax(1),'ydir','rev'); ylabel('Depth');
    ylabel('SpeedJJ','parent',ax(2));
end

end
        

%public subfunctions
% circ_mean.m
% (c) Phillipp Berens, 2009
% https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-
function [mu ul ll] = circ_mean(alpha, w, dim)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit 
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 3
  dim = 1;
end

if nargin < 2 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
  t = circ_confmean(alpha,0.05,w,[],dim);
  ul = mu + t;
  ll = mu - t;
end
end

%
% runcirc_mean.m
% (c) David Cade, 2016
function y = runcirc_mean (x,m) %input in radians, not fast

if sum(size(x)>1)<2
X = buffer(x,2*m,2*m-1,x(1)*ones(2*m-1,1));
X = X(:,m:end);
X2 = rot90(buffer(flipud(x(end-(2*m-1):end)),2*m,2*m-1,x(end)*ones(2*m-1,1)),2);
X = [X X2(:,2:m)];
y = circ_mean(X);
if size(x,1)>size(x,2); y = y'; end
else
    y = nan(size(x));
    for i = 1:size(x,2)
        x2 = x(:,i);
        X = buffer(x2,2*m,2*m-1,x2(1)*ones(2*m-1,1));
        X = X(:,m:end);
        X2 = rot90(buffer(flipud(x2(end-(2*m-1):end)),2*m,2*m-1,x2(end)*ones(2*m-1,1)),2);
        X = [X X2(:,2:m)];
        y(:,i) = circ_mean(X)';
    end
end
end

%
% runmean.m
% (c) Jos van der Geest 2006
% http://www.mathworks.com/matlabcentral/fileexchange/10113-runmean
function Y = runmean(X, m, dim, modestr) 
% RUNMEAN - Very fast running mean (aka moving average) filter
%   For vectors, Y = RUNMEAN(X,M) computes a running mean (also known as
%   moving average) on the elements of the vector X. It uses a window of
%   2*M+1 datapoints. M an positive integer defining (half) the size of the
%   window. In pseudo code: 
%     Y(i) = sum(X(j)) / (2*M+1), for j = (i-M):(i+M), and i=1:length(X) 
%
%   For matrices, Y = RUNMEAN(X,M) or RUNMEAN(X,M,[]) operates on the first   
%   non-singleton dimension of X. RUNMEAN(X,M,DIM) computes the running
%   mean along the dimension DIM.
%
%   If the total window size (2*M+1) is larger than the size in dimension
%   DIM, the overall average along dimension DIM is computed.
%
%   As always with filtering, the values of Y can be inaccurate at the
%   edges. RUNMEAN(..., MODESTR) determines how the edges are treated. MODESTR can be
%   one of the following strings:
%     'edge'    : X is padded with first and last values along dimension
%                 DIM (default)
%     'zero'    : X is padded with zeros
%     'mean'    : X is padded with the mean along dimension DIM 
%
%   X should not contains NaNs, yielding an all NaN result. NaNs can be
%   replaced by using, e.g., "inpaint_nans" created by John D'Errico.
%
%   Examples
%     runmean([1:5],1) 
%       % ->  1.33  2  3  4 4.67
%     runmean([1:5],1,'mean') 
%       % ->  2 2 3 4 4
%     runmean([2:2:10],1,1) % dimension 1 is larger than 2*(M=1)+1 ...
%       % -> 2 4 6 8 10
%     runmean(ones(10,7),3,2,'zero') ; % along columns, using mode 'zero'
%     runmean(repmat([1 2 4 8 NaN 5 6],5,1),2,2) ; 
%       % -> all NaN result
%     A = rand(10,10) ; A(2,7) = NaN ;
%     runmean(A,3,2) ; 
%       % -> column 7 is all NaN
%     runmean(1:2:10,100) % mean
%       % -> 5 5 5 5 5
%
%   This is an incredibly fast implementation of a running mean, since
%   execution time does not depend on the size of the window.
%
%   See also MEAN, FILTER

% for Matlab R13
% version 3.0 (sep 2006)
% Jos van der Geest
% email: jos@jasen.nl

% History:
%   1.0 (2003) created, after a snippet from Peter Acklam (?)
%   1.1 (feb 2006) made suitable for the File Exchange (extended help and
%       documentation)
%   1.2 (feb 2006) added a warning when the window size is too big
%   1.3 (feb 2006) improved help section
%   2.0 (sep 2006) working across a dimension of a matrix. 
%   3.0 (sep 2006) several treatments of the edges. 

% Acknowledgements: (sep 2006) Thanks to Markus Hahn for the idea of
% working in multi-dimensions and the way to treat edges. 

error(nargchk(2,4,nargin)) ;

if ~isnumeric(m) || (numel(m) ~= 1) || (m < 0) || fix(m) ~= m,
    error('The window size (M) should be a positive integer') ;
end

if nargin == 2,
    dim = [] ;
    modestr = 'edge' ;
elseif nargin==3,
    if ischar(dim),
        % no dimension given
        modestr = dim ;
        dim = [] ;
    else 
        modestr = 'edge' ;
    end
end

modestr = lower(modestr) ;

% check mode specifier
if ~ismember(modestr,{'edge','zero','mean'}),
    error('Unknown mode') ;
end

szX = size(X) ;
if isempty(dim),
    dim = min(find(szX>1)) ;
end

if m == 0 || dim > ndims(X),
    % easy
    Y = X ;
else
    mm = 2*m+1 ;
    if mm >= szX(dim),
        % if the window is larger than X, average all
        sz2 = ones(size(szX)) ; 
        sz2(dim) = szX(dim) ;
        Y = repmat(mean(X,dim),sz2) ;
    else
        % here starts the real stuff
        % shift dimensions so that the desired dimensions comes first
        [X, nshifts] = shiftdim(X, dim-1); 
        szX = size(X) ;
        % make the rest of the dimensions columns, so we have a 2D matrix
        % (suggested of Markus Hahn)
        X = reshape(X,szX(1),[]) ; 
        % select how to pad the matrix
        switch (modestr),
            case 'edge'
                % pad with first and last elements
                Xfirst = repmat(X(1,:),m,1) ;
                Xlast = repmat(X(end,:),m,1) ;
            case 'zero'
                % pad with zeros
                Xfirst = zeros(m,1) ;
                Xlast= zeros(m,1) ;
            case 'mean',
                % pad with the average
                Xfirst = repmat(mean(X,1),m,1) ;
                Xlast = Xfirst ;
        end        
        % pad the array
        Y = [zeros(1,size(X,2)) ; Xfirst ; X ; Xlast] ;       
        % the cumsum trick (by Peter Acklam ?)
        Y = cumsum(Y,1) ;
        Y = (Y(mm+1:end,:)-Y(1:end-mm,:)) ./ mm ;
        
        % reshape into original size
        Y = reshape(Y,szX)   ;
        % and re-shift the dimensions
        Y = shiftdim(Y,ndims(Y)-nshifts) ;
    end
end

% =====================
%  CODE OF VERSION 1.3 
% =====================

% function Y = runmean(X,m) ;
% % RUNMEAN - Very fast running mean filter for vectors
% %   Y = RUNMEAN(X,M) computes a running mean on vector X using a window of
% %   2*M+1 datapoints. X is a vector, and M an positive integer defining
% %   (half) the size of the window. In pseudo code:
% %     Y(i) = sum(X(j)) / (2*M+1), for j = (i-M):(i+M), and i=1:length(X)
% %
% %   If the total window size (2M+1) is larger than the length of the vector, the overall
% %   average is returned.
% %
% %   Example:
% %     runmean(1:10,1) % ->
% %     [1.3333 2 3 4 5 6 7 8 9 9.6667]
% %
% %   This is an incredibly fast implementation of a running average, since
% %   execution time does not depend on the size of the window.
% %
% %   X should not contains NaNs (a NaN will result in a all NaN result)
% %   At both ends the values of Y can be inaccurate, as the first and last
% %   values of X are used multiple times. 
% %
% %   See also MEAN
% 
% % for Matlab R13
% % version 1.3 (feb 2006)
% % Jos van der Geest
% % email: jos@jasen.nl
% 
% % History:
% % 1.0 (2003) created, after a snippet from Peter Acklam (?)
% % 1.1 (feb 2006) made suitable for the File Exchange (extended help and
% % documentation)
% % 1.2 (feb 2006) added a warning when the window size is too big
% % 1.3 (feb 2006) improved help section
% 
% error(nargchk(2,2,nargin)) ;
% 
% sz = size(X) ;
% 
% if numel(sz) ~= 2 || (min(sz) ~= 1),
%     error('X should be a vector') ;
% end
% 
% if any(isnan(X)),
%     error('NaNs cannot be dealt with') ;
% end
% 
% if ~isnumeric(m) || (numel(m) ~= 1) || (m < 0) || fix(m) ~= m,
%     error('The window size (M) should be a positive integer') ;
% elseif m == 0,
%     Y = X ;
%     return ;
% end
% 
% mm = 2*m+1 ;
% 
% if mm >= prod(sz),
%     % if the window is larger than X, average all
%     warning('Window size is larger than the length of the vector.')
%     Y = repmat(mean(X),sz) ;
% else
%     % the cumsum trick ...
%     Y = [repmat(X(1),m,1) ; X(:) ; repmat(X(end),m,1)] ;
%     Y = [0 ; cumsum(Y)] ;
%     Y = (Y(mm+1:end)-Y(1:end-mm)) / mm ;
%     Y = reshape(Y,sz) ;
% end
end
