function [varargout] = pvalue(varargin)
%PVALUE Unwrap a PhasorArray to its underlying 3D array; pass anything else
%   through unchanged.
%
%   This is the seam between the legacy array-3D functions and the PhasorArray
%   class (which arrived later). A PhasorArray returns its .value; any raw
%   array (double, sym, sdpvar, ndsdpvar) is returned as-is so the polymorphic
%   element types are preserved. Accepts multiple inputs (returns a cell).
if nargin==1
    if isa(varargin{1},'PhasorArray')
        varargout{1}=value(varargin{1});
    else
        varargout{1}=varargin{1};
    end
else
    varargout=cellfun(@pvalue,varargin,'UniformOutput',0);
end

