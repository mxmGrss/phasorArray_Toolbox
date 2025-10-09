function [varargout] = pvalue(varargin)
%pavlue output array value of phasor array, or directly the array,
%depending on the input
if nargin==1
    if isa(varargin{1},'PhasorArray')
        varargout{1}=value(varargin{1});
    elseif isa(varargin{1},'double')
        varargout{1}=varargin{1};
    else
        error('non handled class')
    end
else
    varargout=cellfun(@pvalue,varargin,'UniformOutput',0);
end

