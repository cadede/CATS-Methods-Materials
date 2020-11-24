function       [X,headers] = read_csv(fname,hdr,rr,delim)
%
%     [X,fields] = read_csv(fname,hdr,rr,delim)
%     Read data from a CSV file. The file is assumed to have a list of column
%     names in the first line followed by lines of data. Each line must have the
%     same number of elements as the first line. Empty cells are allowed.
%     read_csv returns a structure array of strings with the field names taken
%     from the headers. For a rectangular cell array of strings (nentries x nfields)
%     use struct2cell(X)
%     To convert a column containing only numbers to a vector of numbers,
%     use: str2num(strvcat(X(:).field))  where field is the name of a field
%     in structure X.
%     To read only the header of the file, put hdr=1.
%     If the CSV file does not have a header, put hdr=-1. In this case,
%     the function will return a cell array with the contents of the file.
%
%		Example:
%		 S=read_csv('sensor_names.csv')
% 	    Returns: S is a structure with fields matching the names on the
%               first line of sensor_names.csv.
%
%     Valid: Matlab, Octave
%     markjohnson@st-andrews.ac.uk
%     last modified: 12 July 2017

if nargin<1,
   help read_csv ;
   return
end

if nargin<2 | isempty(hdr),
   hdr = 0 ;
end

if nargin<3,
   rr = [] ;      % default is to read all rows
end

if nargin<4 || isempty(delim),
   delim = [44,9] ;           % ',' and '\t'
else
   delim = abs(delim) ;
end

X = {} ;
headers = {} ;
f = fopen(fname,'rt') ;
if f<0,
   fprintf('Cannot open file %s\n',fname) ;
   return
end

lk = 1 ;                                 % line count
if hdr>0,
   hh = fgetl(f) ;
   lk = lk+1 ;
   headers = parseline(hh,delim) ;     % break header into field names
end

nfields = length(headers) ;

% skip lines if rr(1)>0,
if length(rr)>1,
   while lk<rr(1),
      fgetl(f) ;
      lk = lk+1 ;
   end
   rr = rr(2) ;
end

% read remainder of the file into a cell of strings
S = {} ;
while 1,
   if ~isempty(rr) && lk>rr, break, end
   ss = fgetl(f) ;
   lk = lk+1 ;
   if isempty(ss) || ss(1)<0, break, end
   x = parseline(ss,delim) ;
   if isempty(S) && nfields==0,
      nfields = length(x) ;
   end
   if length(x)==nfields,
      [S{1:nfields,end+1}] = deal(x{:}) ;
      continue
   end

   if length(x)>nfields,
      fprintf('Too many columns on line %d of file %s\n',lk,fname) ;
      fclose(f) ;
      return
   elseif ~isempty(x) && length(x)<nfields,
      fprintf('Too few columns on line %d of file %s\n',lk,fname) ;
      fclose(f) ;
      return
   end
end

fclose(f) ;
if isempty(headers), X = S' ; return, end

% make structure with field names taken from the header
headers = strip_quotes(headers) ;
% check field names to make sure they have no illegal characters
for k=1:length(headers),
   h = headers{k} ;
   h(ismember(h,'-/\().:, ')) = '_' ;
   headers{k} = h ;
end
if ~isempty(S),
   X = cell2struct(S,headers,1) ;
end
return


function    x = parseline(s,delim)
%
%
x = {} ;
s = deblank(s) ;
if ismember(abs(s(end)),delim),
   endfield = 1 ;
else
   endfield = 0 ;
end

while ~isempty(s),
   if s(1)=='"',
      k = strfind(s,'",') ;
      if isempty(k),
         if s(end)=='"',
            x{end+1} = s ;
         end
         break ;
      else
         x{end+1} = s(1:k(1)) ;
         s = s(k(1)+2:end) ;
      end
   else
      k = find(ismember(abs(s),delim)) ;
      if isempty(k),
         x{end+1} = s(1:end) ;
         break ;
      else
         x{end+1} = s(1:k(1)-1) ;
         s = s(k(1)+1:end) ;
      end
   end
end

if endfield==1,
   x{end+1} = '' ;
end
return
