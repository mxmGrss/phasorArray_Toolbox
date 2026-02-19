function PhDet = detLeibnizHmc(A)
%DETLEIBNIZHMC Compute the determinant of a PhasorArray purely in the harmonic domain.
%
%   PhDet = DETLEIBNIZHMC(A) computes the determinant of the periodic matrix A(t)
%   using the Leibniz formula, operating entirely in the harmonic (phasor) domain.
%
%   Unlike the standard DET method (which uses IFFT -> pointwise det -> FFT),
%   this function computes det(A(t)) via the classical Leibniz expansion:
%
%       det(A) = sum_{sigma in S_n} sgn(sigma) * prod_{i=1}^{n} A_{i, sigma(i)}
%
%   Each product is a chain of Cauchy products (phasor convolutions), and the
%   sum is a phasor addition. This yields an EXACT result in the harmonic domain
%   with no FFT truncation artifacts.
%
%   OUTPUT HARMONIC ORDER:
%       If A has harmonic order h, the result has harmonic order up to n*h,
%       where n = size(A,1). This is exact: the Leibniz product of n scalar
%       phasors each of order h yields a scalar phasor of order n*h.
%
%   COMPLEXITY:
%       The computation involves n! permutations, each requiring (n-1) Cauchy
%       products. This makes the function practical only for small matrices:
%           n=2 :  2 permutations,  2 products   -> trivial
%           n=3 :  6 permutations, 12 products   -> fast
%           n=4 : 24 permutations, 72 products   -> feasible
%           n>=5: n! grows factorially            -> use det() instead
%
%   KEY ADVANTAGE:
%       This function works with SDPVAR (ndsdpvar) and SYMBOLIC PhasorArrays,
%       where the FFT-based det() cannot be used. This enables formulating
%       determinant constraints (e.g., det(A(t)) > 0) in SDP/LMI problems.
%
%   INPUT:
%       A  - A square PhasorArray object (n x n x (2h+1)).
%
%   OUTPUT:
%       PhDet - A scalar PhasorArray (1 x 1) representing det(A(t)).
%
%   EXAMPLES:
%       % Compute determinant of a 2x2 PhasorArray
%       A = PhasorArray.random(2, 2, 5);
%       d = detLeibnizHmc(A);
%
%       % Compare with FFT-based det
%       d_fft = det(A);
%       fprintf('Harmonic order: Leibniz=%d, FFT=%d\n', d.h, d_fft.h);
%
%       % Works with symbolic PhasorArrays
%       Asym = PhasorArray.sym(2, 2, 3, "M");
%       dsym = detLeibnizHmc(Asym);
%
%       % Works with YALMIP sdpvar PhasorArrays
%       Asdp = PhasorArray.ndsdpvar(2, 2, 3);
%       dsdp = detLeibnizHmc(Asdp);
%
%   ALGORITHM:
%       1. Enumerate all n! permutations of {1, ..., n}.
%       2. For each permutation sigma, compute the sign sgn(sigma).
%       3. Compute the product: sgn(sigma) * A{1,sigma(1)} * ... * A{n,sigma(n)}
%          using PhasorArray scalar multiplication (Cauchy product).
%       4. Sum all n! signed products.
%
%   See also: DET, PHASORDET, INV, TRACE.

arguments
    A
end

% Validate input
n = size(A, 1);
assert(n == size(A, 2), ...
    'detLeibnizHmc:notSquare', ...
    'Input must be a square PhasorArray, got %d x %d.', n, size(A, 2));

% Warn for large matrices
if n >= 5
    warning('detLeibnizHmc:largeDimension', ...
        'Matrix dimension n=%d results in %d! = %d permutations. Consider using det() instead.', ...
        n, n, factorial(n));
end

% Special case: scalar
if n == 1
    PhDet = A{1, 1};
    return
end

% Special case: 2x2 (avoid perms overhead)
if n == 2
    PhDet = A{1,1} * A{2,2} - A{1,2} * A{2,1};
    return
end

% General case: Leibniz expansion
P = perms(1:n);          % n! x n matrix of all permutations
nPerms = size(P, 1);

PhDet = PhasorArray(0);  % accumulator (scalar, h=0)

for p = 1:nPerms
    sigma = P(p, :);
    s = permutationSign(sigma);
    
    % Build the product: A{1,sigma(1)} * A{2,sigma(2)} * ... * A{n,sigma(n)}
    term = A{1, sigma(1)};
    for i = 2:n
        term = term * A{i, sigma(i)};
    end
    
    PhDet = PhDet + s * term;
end

end


function s = permutationSign(sigma)
%PERMUTATIONSIGN Compute the sign (parity) of a permutation.
%   s = PERMUTATIONSIGN(sigma) returns +1 if sigma is an even permutation,
%   -1 if sigma is an odd permutation. The sign is determined by counting
%   the number of transpositions (cycles of even length).


n = length(sigma);
visited = false(1, n);
nTranspositions = 0;

for i = 1:n
    if ~visited(i)
        % Follow the cycle starting at i
        cycleLen = 0;
        j = i;
        while ~visited(j)
            visited(j) = true;
            j = sigma(j);
            cycleLen = cycleLen + 1;
        end
        % A cycle of length k contributes (k-1) transpositions
        nTranspositions = nTranspositions + (cycleLen - 1);
    end
end

if mod(nTranspositions, 2) == 0
    s = 1;
else
    s = -1;
end

end
