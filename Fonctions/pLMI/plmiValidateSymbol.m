function plmiValidateSymbol(M,requireReal)
%PLMIVALIDATESYMBOL Check an affine Hermitian polynomial before compilation.
%   plmiValidateSymbol(M,true) additionally requires real-valued time data.
%   Defects are measured on the YALMIP basis coefficients, with a relative
%   tolerance of 1e-12.
arguments
    M (:,:) PhasorArray
    requireReal (1,1) logical = false
end
if size(M,1) ~= size(M,2)
    error('pLMI:notSquare','A matrix inequality requires a square symbol.');
end
C = value(M);
vectorC = reshape(C,[],1);
if isnumeric(C)
    basisC = vectorC;
else
    if ~all(is(vectorC,'linear'))
        error('pLMI:notAffine','Coefficients must be affine in the decision variables.');
    end
    basisC = getbase(vectorC);
end
if ~all(isfinite(nonzeros(basisC)))
    error('pLMI:nonfinite','Finite coefficients are required.');
end
tolerance = 1e-12*max(1,norm(basisC,'fro'));
n = size(M,1); count = size(C,3);
% Flatten first: ndsdpvar has no conj on 3-D arrays.
for k = 1:count
    current = reshape(vectorC((k-1)*n*n+(1:n*n)),n,n);
    opposite = reshape(vectorC((count-k)*n*n+(1:n*n)),n,n);
    if basisNorm(current-opposite') > tolerance
        error('pLMI:notHermitian','The symbol must be Hermitian in time.');
    end
    if requireReal && basisNorm(current-conj(opposite)) > tolerance
        error('pLMI:notReal','The real Gram basis requires a real-valued symbol.');
    end
end
end

function magnitude = basisNorm(coefficients)
vector = reshape(coefficients,[],1);
if isnumeric(vector), basis = vector; else, basis = getbase(vector); end
magnitude = norm(basis,'fro');
end
