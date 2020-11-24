function 	S = csv2struct(dirpath,fname)

%		S = csv2struct(fname)		% file is on the search path
%		or
%		S = csv2struct(dirpath,fname)		% specify where the file is
%		Read a CSV metadata file and convert it into a metadata structure.
%		A metadata file is a text file containing a line for each metadata
%		entry. The first comma-separated field in each line is the name of the
%		entry which is a dot-separated list of structure elements, e.g.,
%		'animal.dbase.url'. The last field in each line contains the value to
%		be assigned to this metadata entry. The value can be a string or number
%		but is always saved as a string in the structure - it is up to downstream
%		users of the metadata to parse/decode the entries.
%
%		Inputs:
%		dirpath is a string containing the path to the file. If the file is on
%		 the search path, skip this argument.
%		fname is the name of the metadata file. If the name does not include a .csv
%		 suffix, this will be added automatically.
%
%		Returns:
%		S is a metadata structure populated from the file.
%
%		Example:
%		 S = csv2struct('testset1')
% 	    returns: S with fields including S.depid='md13_134a'.
%
%     Valid: Matlab, Octave
%		Rene Swift (rjs@st-andrews.ac.uk)
%     and markjohnson@st-andrews.ac.uk
%     last modified: 12 July 2017

S = struct ;
if nargin<1
   help csv2struct
   return
end

if nargin==1,
	fname = dirpath ;
	dirpath = [] ;
end

if length(fname)<4 || ~all(fname(end+(-3:0))=='.csv'),
   fname(end+(1:4))='.csv';
end

if ~isempty(dirpath),
	if ismember(dirpath(end),'\/'),
		dirpath = dirpath(1:end-1) ;
	end
	fname=[dirpath '\' fname];	
end

% Check to see if there is a header field
fid = fopen(fname, 'r');
if fid<=0,
	fprintf('Error: cannot find file %s\n',fname) ;
	return
end
	
fmtspec='%s';
C = textscan(fid,fmtspec,1,'Delimiter','\r');
fclose(fid);
x=strfind(C{1},'depid'); 	% Check for header
fid = fopen(fname, 'r');
if ~isempty(x)
   C = textscan(fid,fmtspec,'Delimiter','\r');
else
   C = textscan(fid,fmtspec,'Delimiter','\r','HeaderLines',1);
end
fclose(fid);

T=C{1} ;
for k=1:length(T),		% for each line of the file...
	t = T(k) ;
	t = t{1} ;
   f = t(1:find(t==',',1)-1) ;   % Field name is up to the first comma
	f(find(f=='.')) = '_' ;
	if t(end)=='"',
		v = t(find(t(1:end-1)=='"',1,'last')+1:end-1) ;
	else
		v = t(find(t==',',1,'last')+1:end) ;
	end
	S.(f) = v ;
end
