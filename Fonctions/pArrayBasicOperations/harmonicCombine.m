function Mt = harmonicCombine(M, E)
%HARMONICCOMBINE  Combine the harmonic dimension of a phasor array against a basis.
%
%   Mt = HARMONICCOMBINE(M, E), with M of size [n x m x p] and E of size
%   [p x q], returns [n x m x q]:
%
%       Mt(:,:,k) = sum_j M(:,:,j) * E(j,k)
%
%   E is the harmonic basis on a time grid: exp(1i*k*omega*t) (exponential
%   form) or [sin; cos] (trigonometric form).
%
%   If size(M,3) ~= size(E,1), the larger is centre-truncated to the
%   smaller's count -- equal pages/rows dropped from each end, keeping the
%   DC term (the middle entry) aligned. An odd excess drops from the
%   high-harmonic end. n and m are unchanged.
%
%   Example
%       h = 5; t = linspace(0, 2*pi, 128);
%       E  = exp(1i*(-h:h)'*t);
%       Mt = harmonicCombine(randn(3,3,2*h+1), E);   % 3 x 3 x 128
%
%   See also PhasorArray2time, evalTime, evalp.

arguments
    M {mustBeNumeric}
    E {mustBeNumeric}
end

[n, m, p] = size(M);
q = size(E, 1);

if p ~= q
    excess  = abs(p - q);
    dropLow = floor(excess / 2);
    dropHigh = excess - dropLow;
    if p > q
        M = M(:, :, dropLow+1 : end-dropHigh);
    else
        E = E(dropLow+1 : end-dropHigh, :);
    end
end

% reshape+matmul, not tensorprod(M,E,3,1): toolbox supports R2021b, tensorprod is R2022a+.
Mt = reshape(reshape(M, n*m, []) * E, n, m, []);
end
