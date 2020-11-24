function [ varargout ] = unitconversion( mtype, iunit, ounit, useport )
% UNITCONVERSION internal function containing common algorithm for unit
% conversion. 
  
%   Copyright 2000-2011 The MathWorks, Inc.

error(nargoutchk(0,3, nargout,'struct'));

if ~ischar( mtype )
    error(message('aero:unitconversion:notChar1'));
end

if ~ischar( iunit )
    error(message('aero:unitconversion:notChar2'));
end

if ~ischar( ounit )
    error(message('aero:unitconversion:notChar3'));
end

if ~isnumeric( useport )
    error(message('aero:convlength:notNumeric4'));
end

% load unit conversion data structure
ConvertStruct = getunitdata(mtype);

units = {ConvertStruct.tdata.unit};

% find index for input unit conversion data
inidx = find(strcmpi(iunit,units),1);

if isempty( inidx )
    error(message('aero:unitconversion:unknownInputUnit', mtype, iunit));
end

% find index for output unit conversion data
outidx = find(strcmpi(ounit,units),1);

if isempty( outidx )
    error(message('aero:unitconversion:unknownOutputUnit', mtype, ounit));
end

slope_in = ConvertStruct.tdata(inidx).slope;
slope_out = ConvertStruct.tdata(outidx).slope;
varargout{1} = slope_in/slope_out;

if strcmpi(mtype,'temperature conversion')
    bias_in = ConvertStruct.tdata(inidx).bias;
    bias_out = ConvertStruct.tdata(outidx).bias;
    varargout{2} = ( bias_in - bias_out )/slope_out;
end

if useport
    ports(1).txt = ConvertStruct.tdata(inidx).unit;
    ports(2).txt = ConvertStruct.tdata(outidx).unit;
    if strcmpi(ConvertStruct.tdata(inidx).unit,'naut mi')
        ports(1).txt = 'n.mi';
    end
    if strcmpi(ConvertStruct.tdata(outidx).unit,'naut mi')
        ports(2).txt = 'n.mi';
    end
    if strcmpi(mtype,'temperature conversion')
        varargout{3} = ports;
    else
        varargout{2} = ports;
    end
end
