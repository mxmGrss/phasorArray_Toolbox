function [out] = isphasor(o1)
    %ISPHASOR Checks if each element of input is an instance of PhasorArray.
    %   out = isphasor(o1) returns a logical cell array where each element 
    %   corresponds to whether the element in o1 is a PhasorArray. If there is
    %   only one element in o1, the output is a single logical value.
    %
    %   Input:
    %       - o1: A cell array of objects to check for PhasorArray type.
    %
    %   Output:
    %       - out: Logical cell array or a single logical value indicating 
    %              whether each element of o1 is a PhasorArray.
    %       If o1 is a single element, out is a logical value.
    %       
    %   See also: isa
arguments (Repeating)
    o1
end
out=cell(size(o1));
for ii=1:numel(o1)
    out{ii}=logical(isa(o1{ii},"PhasorArray"));
end
if numel(o1)==1
    out=out{:};
end
end

