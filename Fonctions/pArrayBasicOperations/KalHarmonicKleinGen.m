function [L, Y, info] = KalHarmonicKleinGen(A, C, W, V, E, L0, T, nvp)
%KALHARMONICKLEINGEN  Periodic descriptor Kalman filter via dual GPDRE-LTP.
%
%   Computes the steady-state (T-periodic) Kalman filter for the linear
%   time-periodic (LTP) descriptor system
%
%       E(t) x' = A(t) x + w,        cov(w) = W(t) >= 0
%            y  = C(t) x + v,        cov(v) = V(t) >  0
%
%   WITHOUT inverting E(t) and WITHOUT forming dE/dt. The filter is realised
%   in descriptor (observer) form, also free of inv(E) and dE/dt:
%
%       E(t) xhat' = A(t) xhat + L(t) * ( y - C(t) xhat )
%
%   NB: the mass is OUTSIDE the derivative, mirroring the plant — NOT
%   d/dt(E xhat), which would equal E xhat' + dE/dt xhat and inject a
%   spurious dE/dt term (biased estimator). Proof by the error e = x - xhat:
%
%       E e' = E x' - E xhat' = (A - L C) e + w - L v
%
%   a clean descriptor error loop whose covariance equation is exactly the
%   closed-loop form of (FILT) below. Equivalently, in z = E x the optimal
%   filter reads z' = (A + dE/dt)E^-1 z + L(y - C E^-1 z); substituting
%   zhat = E xhat cancels the dE/dt terms and leaves the realisation above.
%   So: plant and observer are mass-outside; it is the ADJOINT and the
%   covariance equation that are of sandwich (mass-inside) type.
%
%   ========================================================================
%   THEORY — derivation of the generalized filter GPDRE (LTP)
%   ========================================================================
%
%   (1) Explicit filter. With Atil = E^-1 A and Wtil = E^-1 W E^-T, the
%       standard (forward) periodic filter Riccati equation for the
%       state-error covariance Y(t) = cov(x - xhat) reads
%
%         Y' = Atil Y + Y Atil.' - Y C.' V^-1 C Y + Wtil.
%
%   (2) Polynomial (inverse-free) form. Multiply on the left by E and on
%       the right by E.'. Every inverse cancels and one obtains the
%       GENERALIZED FILTER GPDRE:
%
%         E Y' E.' = A Y E.' + E Y A.' - E Y C.' V^-1 C Y E.' + W      (FILT)
%
%       *** The derivative is SANDWICHED: E*(dY/dt)*E.', NOT d/dt(E Y E.').
%       The two differ by the dE/dt terms and only coincide when E is
%       constant. One can also see it through z = E x: the covariance of z
%       is Pz = E Y E.' and its filter equation (for the explicit z-system
%       z' = (A + dE/dt) E^-1 z + w) generates dE/dt terms on BOTH sides of
%       d/dt(E Y E.') that cancel exactly, leaving E Y' E.'. Forgetting this
%       and using the product derivative silently solves the Riccati of a
%       DIFFERENT system — the residual of (FILT) does not vanish (see
%       test_phase3, Test 4, which caught precisely this).
%
%       Filter gain (innovation injection in the descriptor observer above):
%
%         L(t) = E(t) Y(t) C(t).' V(t)^-1
%
%       (Explicit-form gain, if one insists on xhat' = ...: Lx = Y C.' V^-1.)
%
%   (3) Why the sandwich is structural: adjoint duality. The adjoint of the
%       mass-OUTSIDE descriptor E x' = A x is the mass-INSIDE system
%       d/dt(E.' lambda) = -A.' lambda. Control GPDREs for the mass-outside
%       form carry the product derivative d/dt(E.'XE) (RicHarmonicKleinGen,
%       'product'); GPDREs for the mass-inside form carry the sandwiched
%       derivative E.'(dX/dt)E ('sandwich'). Filtering being the control of
%       the adjoint, the filter GPDRE is necessarily of sandwich type.
%
%   (4) Duality mapping = transposition + FORWARD direction. (FILT) is a
%       forward equation; plain transposition alone would hand the control
%       solver an ANTICAUSAL problem — a symmetric, positive, plausible and
%       wrong answer (test_phase3, Test 2 demonstrates a residual ratio
%       ~1e11). The causality is restored by flipping the sign of the
%       derivative term: with direction='forward' the control solver solves
%
%         -Ec.'(dX/dt)Ec + Ac.' X Ec + Ec.' X Ac
%                        - Ec.' X Bc Rc^-1 Bc.' X Ec + Qc = 0    (CTRL-fwd-sdw)
%
%       and the plain transposed mapping
%
%         Ac = A.'   Bc = C.'   Ec = E.'   Qc = W   Rc = V
%
%       turns (CTRL-fwd-sdw) exactly into (FILT) with Y = X — no time
%       reversal, no flip back. (Equivalent legacy route: transposition +
%       TIME REVERSAL trretro(.) = (.)(-t).' with the backward solver and
%       recovery Y = retro(X); the two routes build the same harmonic
%       operator and agree EXACTLY — kept as a cross-validation test.)
%
%   (5) HSS realisation. The dual problem is solved by the existing chain
%
%         RicHarmonicKleinGen(derivativeForm='sandwich', direction='forward')
%           -> lyapG(...) -> SylvHarmonicGen:
%              [ -T(Eb.'(x)Ea)*N + T(M_S) ] vec(X) = -vec(C)
%
%       where N = I (x) diag(j*k*omega) and T(.) is the Toeplitz-block
%       operator. Note the operator ordering T(M_E)*N (sandwich) versus
%       N*T(M_E) (product): that single transposition of composition order
%       is the entire difference between primal control and dual filtering
%       in the HSS algebra; the SIGN of that term is the causality
%       (backward = control, forward = filtering).
%
%   ========================================================================
%
%   Usage:
%     [L, Y, info] = KalHarmonicKleinGen(A, C, W, V)
%     [L, Y, info] = KalHarmonicKleinGen(A, C, W, V, E)
%     [L, Y, info] = KalHarmonicKleinGen(A, C, W, V, E, L0, T, 'maxIter', 50, ...)
%
%   Inputs:
%     A     PhasorArray, n x n   — dynamics
%     C     PhasorArray, p x n   — measurement matrix
%     W     State noise covariance, n x n >= 0 (double or PhasorArray)
%     V     Measurement noise covariance, p x p > 0 (double or PhasorArray)
%     E     Mass matrix, n x n, PhasorArray or numeric (possibly time-
%           varying). [] -> identity (standard periodic Kalman filter).
%     L0    Initial filter gain (n x p); [] -> DC generalized fallback via
%           the dual icare (default: []). Must stabilise the descriptor
%           error loop  E e' = (A - L0*C) e.
%     T     Period (default: 2*pi)
%
%   Options (name-value) — forwarded to RicHarmonicKleinGen:
%     maxIter, h, autoUpdateh, maxh, systemType, updateMethod,
%     thresholdResidual, relChangeThreshold, reduceThreshold,
%     warmStartFraction, stagnationWindow, stagnationRatio,
%     verbose, storeIterates, skipValidate
%     (derivativeForm='sandwich' and direction='forward' are forced
%     internally — see theory.)
%
%   Outputs:
%     L     Filter gain (PhasorArray, n x p), L = E*Y*C.'*V^-1 — computed
%           directly from Y so that (L, Y) are exactly consistent.
%     Y     Periodic state-error covariance (PhasorArray, n x n), solution
%           of (FILT); cov(E*x - E*xhat) = E*Y*E.'
%     info  Diagnostics from the dual control solve (RicHarmonicKleinGen
%           fields), augmented with:
%       .resFiltAbs       ||residual of (FILT)||_F  (recomputed on the
%                         filter side with the original, non-reversed data —
%                         independent check of the dual mapping)
%       .resFiltnorm      .resFiltAbs / ||W||_F
%       .gainConsistency  ||Kc.' - L|| / ||L||: distance between the last
%                         Kleinman-iterate gain (dual-mapped) and the gain
%                         induced by the returned covariance. Of the order
%                         of the final Newton increment, NOT eps.
%
%   See also: RicHarmonicKleinGen, PhasorArray/lyapG, SylvHarmonicGen

arguments
    A    PhasorArray
    C    PhasorArray
    W                                                                       % double or PhasorArray
    V                                                                       % double or PhasorArray
    E                                                    = []               % mass matrix; [] → identity
    L0                                                   = []               % initial gain; [] → DC fallback
    T                            (1,1) double            = 2*pi             % period
    nvp.maxIter              (1,1) {mustBeInteger, mustBePositive} = 100
    nvp.h                                            = []               % [] → auto
    nvp.autoUpdateh          (1,1) logical           = false
    nvp.maxh                                         = []               % [] → h0*20
    nvp.systemType           {mustBeMember(nvp.systemType,    {'rectangle','square'})}    = 'rectangle'
    nvp.updateMethod         {mustBeMember(nvp.updateMethod,  {'adaptive','incremental'})} = 'adaptive'
    nvp.thresholdResidual    (1,1) double            = 1e-6
    nvp.relChangeThreshold   (1,1) double            = 1e-14
    nvp.reduceThreshold      (1,1) double            = 0               % 0 → τ_ric/100
    nvp.warmStartFraction    (1,1) double            = 0.7             % ∈ (0,1]
    nvp.stagnationWindow     (1,1) {mustBeInteger, mustBePositive} = 5
    nvp.stagnationRatio      (1,1) double            = 0.05
    nvp.verbose              (1,1) {mustBeInteger, mustBeNonnegative} = 0
    nvp.storeIterates        (1,1) logical           = false
    nvp.skipValidate         (1,1) logical           = true
end

W = PhasorArray(W);
V = PhasorArray(V);
if isempty(E), E = PhasorArray(eye(size(A, 1))); end
if ~isa(E, 'PhasorArray'), E = PhasorArray(E); end

if size(C, 2) ~= size(A, 1)
    error('PhasorArray:KalHarmonicKleinGen:dimensionMismatch', ...
        'cols(C)=%d must equal n=%d.', size(C, 2), size(A, 1))
end

%% --- Dual (transposed) control problem — see theory §4 ---
% Causality is carried by direction='forward' (sign of the derivative
% operator), so plain transposition suffices: no time reversal, no flip back.
Ac = A.';
Bc = C.';
Ec = E.';
if isempty(L0)
    K0c = [];
else
    if ~isa(L0, 'PhasorArray'), L0 = PhasorArray(L0); end
    K0c = L0.';
end

%% --- Solve the dual control GPDRE (forward + sandwich — see theory §3-§4) ---
opts = namedargs2cell(nvp);
[Kc, Y, info] = RicHarmonicKleinGen(Ac, Bc, W, V, Ec, K0c, T, ...
    opts{:}, 'derivativeForm', 'sandwich', 'direction', 'forward');

%% --- Map back to filter quantities ---
% Y = X directly (no retro);  L = E*Y*C.'*V^-1, exactly consistent with Y.
Vinv = inv(V);
L = E * Y * C.' * Vinv;

%% --- Independent filter-side diagnostics ---
% Residual of (FILT), recomputed with the ORIGINAL (non-reversed) data:
% catches any convention slip in the dual mapping (incl. the sandwich vs
% product derivative distinction).
RES  = E*d(Y, T)*E.' - A*Y*E.' - E*Y*A.' + E*Y*C.'*Vinv*C*Y*E.' - W;
info.resFiltAbs  = norm(RES.value, 'fro');
info.resFiltnorm = info.resFiltAbs / (norm(W.value, 'fro') + eps);

% Distance between the dual-mapped last-iterate gain and L: ~ final Newton
% increment (diagnostic only, not an error indicator).
info.gainConsistency = norm(value(Kc.' - L), 'fro') / (norm(value(L), 'fro') + eps);

end % KalHarmonicKleinGen
