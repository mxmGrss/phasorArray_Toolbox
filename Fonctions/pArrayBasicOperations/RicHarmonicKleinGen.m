function [K, S, info] = RicHarmonicKleinGen(A, B, Q, R, E, K0, T, varargin)
%RICHARMONICKLEINGEN  Descriptor periodic Riccati — compatibility shim.
%
%   [K, S, info] = RicHarmonicKleinGen(A, B, Q, R, E, K0, T, ...)
%
%   The descriptor and base solvers are merged; prefer in new code:
%       RicHarmonicKlein(A, B, Q, R, K0, T, 'E', E, ...)   or   hare(...)
%
%   E = [] means identity here (historical behaviour: the descriptor path is
%   still taken). Call RicHarmonicKlein without 'E' for the faster path.
%   Name-value options are forwarded verbatim; see its help.
%
%   See also: RicHarmonicKlein, hare, KalHarmonicKleinGen, PhasorArray/lyapG

arguments
    A
    B
    Q
    R
    E                    = []
    K0                   = []
    T   (1,1) double     = 2*pi
end
arguments (Repeating)
    varargin
end

% [] historically meant "identity", and the descriptor path was still used.
% Materialise it so this entry point stays numerically identical to before;
% RicHarmonicKlein called without 'E' is the way to opt into the fast path.
if isempty(E), E = PhasorArray(eye(size(A, 1))); end

[K, S, info] = RicHarmonicKlein(A, B, Q, R, K0, T, 'E', E, varargin{:});
end
