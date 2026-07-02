function [K, X, info] = hare(A, B, Q, R, options)
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
%   Name-value options:
%     E               Descriptor mass matrix; [] → base HARE (default []).
%     direction       'backward' (control, default) | 'forward' (filter /
%                     covariance integration, the dual HARE).
%     derivativeForm  'product' (default) | 'sandwich' (where dE/dt lives).
%     K0              Initial stabilising gain; [] → grow-h LQR-Toeplitz
%                     fallback inside the solver (default []).
%     T               Period (default 2*pi).
%     h               Fixed harmonic truncation; [] → auto (default []).
%     autoUpdateh     Adaptive h-refinement (default true).
%     maxIter         Max Kleinman iterations (default 100).
%     thresholdResidual  Convergence threshold (default 1e-8).
%     verbose         Iteration print level (default 0).
%
%   Outputs:
%     K     PhasorArray optimal gain.
%     X     PhasorArray Riccati solution.
%     info  Solver diagnostics struct.
%
%   The standard (E = [], backward, product) case routes to RicHarmonicKlein;
%   any descriptor / forward / sandwich request routes to RicHarmonicKleinGen
%   (with E = I when E is omitted).
%
%   See also: lyap, RicHarmonicKlein, RicHarmonicKleinGen, KalHarmonicKleinGen.

arguments
    A
    B
    Q
    R
    options.E                                   = []
    options.direction      {mustBeMember(options.direction,     {'backward','forward'})} = 'backward'
    options.derivativeForm {mustBeMember(options.derivativeForm,{'product','sandwich'})} = 'product'
    options.K0                                  = []
    options.T              (1,1) double         = 2*pi
    options.h                                   = []
    options.autoUpdateh    (1,1) logical        = true
    options.maxIter        (1,1) {mustBeInteger, mustBePositive} = 100
    options.thresholdResidual (1,1) double      = 1e-8
    options.verbose        (1,1) {mustBeInteger, mustBeNonnegative} = 0
end

common = {'maxIter', options.maxIter, 'h', options.h, ...
          'autoUpdateh', options.autoUpdateh, ...
          'thresholdResidual', options.thresholdResidual, ...
          'verbose', options.verbose};

useGen = ~isempty(options.E) ...
      || strcmp(options.direction, 'forward') ...
      || strcmp(options.derivativeForm, 'sandwich');

if useGen
    E = options.E;
    if isempty(E), E = PhasorArray(eye(size(A,1))); end   % E = I → reduces to base HARE
    [K, X, info] = RicHarmonicKleinGen(A, B, Q, R, E, options.K0, options.T, ...
        common{:}, 'direction', options.direction, 'derivativeForm', options.derivativeForm);
else
    [K, X, info] = RicHarmonicKlein(A, B, Q, R, options.K0, options.T, common{:});
end
end
