function x = eml_scalar_eg(varargin)
%MATLAB Code Generation Private Function
%
%   Returns a scalar "example" of the additive combined type of the input
%   arguments.  Specifically, if a = varargin{1}, b = varargin{2}, and so
%   forth, the class of the return value x is
%
%   class(a),                                                 nargin == 1
%   class(cast(0,class(a))+cast(0,class(b))),                 nargin == 2
%   class(cast(0,class(a))+cast(0,class(b))+cast(0,class(c)), nargin == 3
%   etc.
%
%   The result x is complex if any of the arguments are complex.
%
%   When all inputs are float, integer, logical, or char, the output x
%   is guaranteed to satisfy x == 0.
%
%   Opaque, struct, and enumeration inputs are supported, but the first
%   input argument must be nonempty, and all inputs after the first
%   argument are ignored.
%
%   Some consequential behavior to note:
%
%       eml_scalar_eg(false) is logical.
%       eml_scalar_eg(false,false) is double.
%       eml_scalar_eg('a') is char.
%       eml_scalar_eg('a','a') is double.

%   Copyright 2005-2013 The MathWorks, Inc.
%#codegen

eml_allow_enum_inputs;
x = coder.internal.scalarEg(varargin{:});
