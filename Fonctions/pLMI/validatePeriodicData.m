function [A,Q] = validatePeriodicData(A,Q)
%VALIDATEPERIODICDATA Validate square real numeric A and symmetric Q.
arguments
    A
    Q
end
A = validatePeriodicReal(A); Q = validatePeriodicReal(Q);
n = size(A,1);
if size(A,2) ~= n || size(Q,1) ~= n || size(Q,2) ~= n
    error('pLMI:systemDimensions','A and Q must be square matrices of the same size.');
end
plmiValidateSymbol(Q,true);
end
