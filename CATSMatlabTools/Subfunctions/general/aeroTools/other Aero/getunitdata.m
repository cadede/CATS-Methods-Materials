function ConvertStruct = getunitdata( mtype,varargin )
% GETUNITDATA internal function retrieving data structure for unit
% conversion. 
  
%   Copyright 2000-2016 The MathWorks, Inc.

% load unit conversion data structure
if isempty(varargin)
    ConvertStruct = aeroconvertdata(mtype);
elseif length(varargin)==1
    ConvertStruct = aeroconvertdata(mtype,varargin{1});
else
    ConvertStruct = aeroconvertdata(mtype,varargin{1},varargin{2});
end
