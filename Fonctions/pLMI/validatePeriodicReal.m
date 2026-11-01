function A = validatePeriodicReal(A)
%VALIDATEPERIODICREAL Require finite numeric coefficients, real-valued in phase.
arguments
    A
end
A = PhasorArray(A); C = value(A);
if ~isnumeric(C) || any(~isfinite(C(:)))
    error('pLMI:numericSystem','System data must have finite numeric coefficients.');
end
defect = C-conj(flip(C,3));
if norm(defect(:)) > 1e-12*max(1,norm(C(:)))
    error('pLMI:realSystem','System data must be real-valued in time.');
end
end
