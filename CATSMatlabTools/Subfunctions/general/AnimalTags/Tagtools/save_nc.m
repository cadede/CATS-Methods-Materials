function    save_nc(fname,X,varargin)
%
%     save_nc(fname,X,...)
%     Save one or more variables to a NetCDF archive file.
%     Warning, this will overwrite any previous NetCDF file with the same name.
%		it is created. The file is assumed to be in the current working directory 
%		unless a pathname is added to the beginning of fname.
%
%		Inputs:
%		fname is the name of the metadata file. If the name does not include a .nc
%		 suffix, this will be added automatically.
%		X (and any subsequent inputs) is a sensor or metadata structure, or a set of sensor
%		 structures. Only these kind of variables can be saved in a NetCDF file because the 
%		 supporting information in these structures is needed to describe the contents of 
%		 the file. For non-archive and non-portable storage of variables, consider using the 
%		 usual 'save' function in Matlab and Octave.
%
%		Example:
%		 save_nc('dog17_124a',A,M,P,info)
% 	    generates a file dog17_124a.nc and adds variables A, M and P, and a metadata
%		 structure.
%
%     Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     last modified: 12 July 2017


if nargin<2,
   help save_nc
   return
end

% append .nc suffix to file name if needed
if length(fname)<3 || ~all(fname(end+(-2:0))=='.nc'),
   fname(end+(1:3))='.nc';
end

if exist(fname,'file'),
	delete(which(fname)) ;
end
	
if ~isstruct(X),
	fprintf('savenc can only save sensor or metadata structures\n') ;
	return
end
	
if isfield(X,'info'),		% X is a set of sensor structures
	f = fieldnames(X) ;
	for k=1:length(f),
		add_nc(fname,X.(f{k})) ;
	end
else
	add_nc(fname,X) ;
end

% save the remaining variables to the file
for k=1:nargin-2, 
   add_nc(fname,varargin{k}) ;
end
