function [k, f, ushifted] = find2piAntecedant(th, u, firstRotIsNaN)
%FIND2PIANTECEDANT Find 2π antecedent indices with fractional interpolation.
%
%   [K, F] = FIND2PIANTECEDANT(TH) finds for each TH(i) the largest index
%   K(i) such that TH(K(i)) <= TH(i) - 2π, and the fractional part F(i):
%       TH(i) - 2π  ≈  TH(K(i)) + F(i) · (TH(K(i)+1) - TH(K(i)))
%
%   For indices where no such antecedent exists (transient, first revolution),
%   K(i) = 1, F(i) = 0  (unless firstRotIsNaN = true).
%
%   [K, F, USHIFTED] = FIND2PIANTECEDANT(TH, U) also computes the signal
%   values at phase TH - 2π via linear interpolation:
%       USHIFTED(i) = U(K(i)) + F(i) · (U(K(i)+1) - U(K(i)))
%
%   Non-decreasing phases are supported, including an initial plateau.
%   Decreasing or non-finite phases are rejected: the forward scan assumes
%   that the phase never reverses.
%
%   Algorithm: O(n) forward scan using a persistent last_k pointer.
%
%   Inputs:
%     TH           : numeric vector of phase values (non-decreasing after unwrap)
%     U            : (optional) signal vector, same length as TH
%     firstRotIsNaN: (optional, default false) return NaN instead of 1/0
%                    for the transient region
%
%   Outputs:
%     K        : integer index vector (row), same size as TH
%     F        : fractional factor vector (row), same size as TH, ∈ [0, 1)
%     USHIFTED : signal at TH-2π (row), or [] if U not provided
%
%   See also: angularsft, computeAngularPhasors

arguments
    th {mustBeNumeric, mustBeVector, mustBeReal, mustBeFinite}
    u = []
    firstRotIsNaN = false
end

th = th(:)';
n  = length(th);
if any(diff(th)<0)
    error('find2piAntecedant:NotIncreasing','Phase must be non-decreasing.');
end
if ~isempty(u) && (~isvector(u) || numel(u)~=n)
    error('find2piAntecedant:SignalSize','U must be a vector with one value per phase.');
end

k = zeros(size(th));
f = zeros(size(th));

if n <= 1
    if n == 1
        k(1) = 1;
        f(1) = 0;
    end
    ushifted = u;
    if firstRotIsNaN && n==1
        k(:)=NaN; f(:)=NaN;
        if ~isempty(u), ushifted(:)=NaN; end
    end
    return;
end

last_k = 1;

for ii = 1:n
    target = th(ii) - 2*pi;

    while last_k < n && th(last_k) <= target
        last_k = last_k + 1;
    end
    last_k = max(last_k - 1, 1);

    if th(last_k) <= target
        k(ii) = last_k;
        if last_k < n
            f(ii) = (target - th(last_k)) / (th(last_k+1) - th(last_k));
        end
    else
        if firstRotIsNaN
            k(ii) = NaN;
            f(ii) = NaN;
        else
            k(ii) = 1;
            f(ii) = 0;
        end
    end
end

if ~isempty(u)
    ushifted = zeros(size(u));
    for ii = 1:n
        if isnan(k(ii))
            ushifted(ii) = NaN;
        else
            ushifted(ii) = u(k(ii)) + f(ii) * (u(k(ii)+1) - u(k(ii)));
        end
    end
else
    ushifted = [];
end

end
