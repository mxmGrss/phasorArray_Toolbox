function [K, X, info] = hare(A, B, Q, R, nvp)
%HARE  Harmonic Algebraic Riccati Equation solver ((G)HARE) for a periodic system.
%
%   [K, X, info] = HARE(A, B, Q, R) solves the periodic (LTP) control PDRE for
%   the system  x' = A(t) x + B(t) u  as a HARE via harmonic Kleinman iteration:
%
%       dX/dt + A'X + XA - X B R^-1 B' X + Q = 0 ,   K = R^-1 B' X .
%
%   With a descriptor mass matrix E (name-value 'E'), it solves the
%   generalized PDRE (GPDRE) for  E(t) x' = A(t) x + B(t) u  as a GHARE,
%   without inverting E(t) or forming dE/dt:
%
%       d/dt(E'XE) + A'XE + E'XA - E'XB R^-1 B'XE + Q = 0 ,  K = R^-1 B'XE .
%
%   MATLAB-style: A is the object; B, Q, R are positional; everything else is
%   name-value. Mirrors care/icare in spirit.
%
%   Name-value nvp:
%     E               Descriptor mass matrix; [] → base HARE (default []).
%     direction       'backward' (control, default) | 'forward' (filter /
%                     covariance integration, the dual HARE).
%     derivativeForm  'product' (default) | 'sandwich' (where dE/dt lives).
%     K0              Initial stabilising gain; [] → grow-h LQR-Toeplitz
%                     fallback inside the solver (default []).
%     T               Period (default 2*pi).
%     h               Fixed harmonic truncation; [] → auto (default []).
%     maxh            Hard upper bound on h; [] → h0*20 (default []).
%     autoUpdateh     Adaptive h-refinement (default true).
%     maxIter         Max Kleinman iterations (default 100).
%     thresholdResidual  Convergence threshold (default 1e-8).
%     skipValidate    Skip the Floquet check of a provided K0 (default true).
%                     Set false to re-check it and fall back if it fails.
%     verbose         Iteration print level (default 0).
%
%   Deeper solver tuning — updateMethod, stagnationWindow, warmStartFraction,
%   reduceThreshold — stays on RicHarmonicKlein, which this forwards to.
%
%   Outputs:
%     K     PhasorArray optimal gain.
%     X     PhasorArray Riccati solution.
%     info  Solver diagnostics struct.
%
%   All cases route to RicHarmonicKlein, which takes the mass matrix as
%   name-value 'E' and branches internally where the formulations differ.
%
%   See also: lyap, RicHarmonicKlein, RicHarmonicKleinGen, KalHarmonicKleinGen.

arguments
    A
    B
    Q
    R
    nvp.E                                   = []
    nvp.direction      {mustBeMember(nvp.direction,     {'backward','forward'})} = 'backward'
    nvp.derivativeForm {mustBeMember(nvp.derivativeForm,{'product','sandwich'})} = 'product'
    nvp.K0                                  = []
    nvp.T              (1,1) double         = 2*pi
    nvp.h                                   = []
    nvp.maxh                                = []
    nvp.autoUpdateh    (1,1) logical        = true
    nvp.maxIter        (1,1) {mustBeInteger, mustBePositive} = 100
    nvp.thresholdResidual (1,1) double      = 1e-8
    nvp.skipValidate   (1,1) logical        = true
    nvp.verbose        (1,1) {mustBeInteger, mustBeNonnegative} = 0
end

common = {'maxIter', nvp.maxIter, 'h', nvp.h, 'maxh', nvp.maxh, ...
          'skipValidate', nvp.skipValidate, ...
          'autoUpdateh', nvp.autoUpdateh, ...
          'thresholdResidual', nvp.thresholdResidual, ...
          'verbose', nvp.verbose};

% Forward 'E' only if given: an absent E keeps the fast path (~1.6x cheaper
% than solving with E = I). direction/derivativeForm work on both paths.
if isempty(nvp.E)
    mass = {};
else
    mass = {'E', nvp.E};
end

[K, X, info] = RicHarmonicKlein(A, B, Q, R, nvp.K0, nvp.T, common{:}, ...
    'direction', nvp.direction, 'derivativeForm', nvp.derivativeForm, mass{:});
end
