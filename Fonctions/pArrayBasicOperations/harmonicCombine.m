function Mt = harmonicCombine(M, E)
%HARMONICCOMBINE  Combine the harmonic dimension of a phasor array against a basis.
%
%   Mt = HARMONICCOMBINE(M, E) with M of size [n x m x p] and E of size [p x q]
%   returns [n x m x q], where Mt(:,:,k) is the combination of the p slices of
%   M weighted by column k of E:
%
%       Mt(:,:,k) = sum_j M(:,:,j) * E(j,k)
%
%   E holds the harmonic basis sampled on a grid — exp(1i*k*omega*t) for the
%   exponential form, [sin; cos] for the trigonometric one — so this is the
%   evaluation step of every time-domain reconstruction.
%
%   A mode-3 contraction is a plain matrix product on the unfolding, so this
%   needs no version gate, unlike tensorprod(M, E, 3, 1) which is R2022a while
%   the toolbox supports R2021b. Same result, same run time.
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

[n, m, ~] = size(M);
Mt = reshape(reshape(M, n*m, []) * E, n, m, []);
end
