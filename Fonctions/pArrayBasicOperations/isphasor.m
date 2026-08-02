function [out] = isphasor(pA1)
    %ISPHASOR Checks if each element of input is an instance of PhasorArray.
    %   out = isphasor(pA1) returns a logical cell array where each element 
    %   corresponds to whether the element in pA1 is a PhasorArray. If there is
    %   only one element in pA1, the output is a single logical value.
    %
    %   Input:
    %       - pA1: A cell array of objects to check for PhasorArray type.
    %
    %   Output:
    %       - out: Logical cell array or a single logical value indicating 
    %              whether each element of pA1 is a PhasorArray.
    %       If pA1 is a single element, out is a logical value.
    %       
    %   See also: isa
arguments (Repeating)
    pA1
end
out=cell(size(pA1));
for ii=1:numel(pA1)
    out{ii}=logical(isa(pA1{ii},"PhasorArray"));
end
if numel(pA1)==1
    out=out{:};
end
end

