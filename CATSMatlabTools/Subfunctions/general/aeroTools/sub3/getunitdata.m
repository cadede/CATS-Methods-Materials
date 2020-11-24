function ConvertStruct = getunitdata( mtype )
% GETUNITDATA internal function retrieving data structure for unit
% conversion. 
  
%   Copyright 2000-2006 The MathWorks, Inc.

error(nargoutchk(0,1, nargout,'struct'));

% load unit conversion data structure
ConvertStruct = aeroconvertdata(mtype);
