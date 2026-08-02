function [res, info] = lyapG(pA1, pA2, pA3, pA4, pA5, nvp)
%LYAPG  Periodic generalized (descriptor) Lyapunov / Sylvester solver.
%
%   Solves descriptor-form periodic equations WITHOUT inverting E(t) and
%   WITHOUT forming dE/dt (handled implicitly by the Toeplitz algebra).
%
%   Generalized Lyapunov mode:  [P, info] = lyapG(A, Q, E, ...)
%     Solves:  d/dt(E.'(t) P E(t)) + A.'(t) P E(t) + E.'(t) P A(t) + Q = 0
%     for the descriptor system E(t) x' = A(t) x, with V = x.' (E.' P E) x.
%     A, Q, E must be same-size square PhasorArray (or convertible).
%     E = [] defaults to identity (→ standard lyap equation).
%
%   Generalized Sylvester mode: [X, info] = lyapG(A, B, C, Ea, Eb, ...)
%     Solves:  d/dt(Ea(t) X Eb(t)) + A(t) X Eb(t) + Ea(t) X B(t) + C = 0
%     A, Ea square nA x nA; B, Eb square nB x nB; C is nA x nB.
%     Ea = [] or Eb = [] default to identity. If BOTH are empty, use lyap.
%
%   Mode detection: Lyapunov when pA4 and pA5 are omitted/empty (3rd arg = E),
%   Sylvester otherwise (3rd arg = C).
%
%   Mass matrices Ea/Eb may be time-varying (PhasorArray): the solver never
%   needs dE/dt explicitly.
%
%   ======================================================================
%   CHOOSING derivativeForm — 'product' vs 'sandwich'
%   ======================================================================
%   The two forms differ by WHERE the dE/dt terms live. They coincide when
%   the mass matrices are constant (an LTI test cannot discriminate them).
%
%   'product' (default) — derivative of the full congruence:
%       d/dt( E.'(t) P E(t) ) + A.'PE + E.'PA + Q = 0
%   This is the PRIMAL Lyapunov equation of the mass-outside descriptor
%       E(t)*dx/dt = A(t)*x.
%   It comes from the Lyapunov function
%       V(x,t) = x.' ( E.'(t) P(t) E(t) ) x
%   whose derivative along trajectories differentiates the whole product:
%   the terms dE/dt.'*P*E + E.'*dP/dt*E + E.'*P*dE/dt are PHYSICALLY part
%   of the equation (handled implicitly by the Toeplitz algebra). Then
%       dV/dt = -x.'*Q*x,  and Q > 0 with P > 0 certifies stability;
%   x0.'*(E.'PE)(t0)*x0 is the cost-to-go int( x.'Qx )dt.
%   Use for: stability certificates, LQR/GPDRE descriptor (control side).
%
%   'sandwich' — derivative of P alone, mass matrices outside:
%       E.'(t) * dP/dt * E(t) + A.'PE + E.'PA + Q = 0
%   No dE/dt term appears. This is the form arising on the DUAL/adjoint
%   side: the dual of the mass-outside descriptor E*dx/dt = A*x is the
%   mass-INSIDE adjoint d/dt(E.'*lambda) = -A.'*lambda, and pushing the
%   mass through the derivative cancels the dE/dt contributions. Concretely
%   it is the filtering GPDRE structure: in the covariance substitution
%   Pz = E*Y*E.' the dE/dt terms cancel on both sides, leaving the
%   sandwiched derivative E*dY/dt*E.' (NOT d/dt(E*Y*E.')).
%   Use for: Kalman / error-covariance equations of the descriptor filter.
%
%   In the harmonic domain the difference is only the operator ordering:
%       'product' :  N * T(M_E)    (derivative applied AFTER the product)
%       'sandwich':  T(M_E) * N    (derivative applied BEFORE the product)
%
%   The 'direction' option is ORTHOGONAL to this choice: derivativeForm
%   selects what the derivative acts on, direction selects its sign
%   (backward = cost-to-go type, forward = covariance type). All four
%   combinations are valid; e.g. the descriptor Kalman covariance is
%   typically (forward, sandwich).
%   ======================================================================
%
%   Options (name-value) — identical to PhasorArray/lyap:
%     T                   Period of the periodic system (default: 2*pi)
%     h                   Fixed harmonic truncation order; [] = auto (default: [])
%     thresholdResidual   Convergence threshold on relative residual (default: 1e-6)
%     autoUpdateh         Enable adaptive h-refinement loop (default: false)
%     maxh                Hard upper bound on h during adaptive loop; [] = h0*20 (default: [])
%     stagnationWindow    Look-back window for stagnation detection (default: 5)
%     stagnationRatio     Min relative improvement to avoid stagnation (default: 0.05)
%     verbose             Print iteration table to console (default: false)
%     storeResidualPhasor Store full residual PhasorArray in info (default: false)
%     systemType          'rectangle' (overdetermined, default) or 'square' (Galerkin)
%     updateMethod        'adaptive' (regime-based jump, default) or 'incremental' (h+1)
%     derivativeForm      'product' (default) or 'sandwich' — what the
%                         derivative acts on. See the dedicated section
%                         "CHOOSING derivativeForm" above.
%     direction           'backward' (default) or 'forward' — sign of the
%                         derivative term. backward = cost-to-go type (as
%                         written above), forward = covariance type, e.g.
%                         d/dt(E.'PE) = A.'PE + E.'PA + Q. Orthogonal to
%                         derivativeForm (see section above).
%
%   Returns:
%     res   PhasorArray solution
%     info  Solver diagnostics struct — same fields as PhasorArray/lyap:
%       .status .statusMsg .resrelnorm .resnorm .h .h_history .resrel_history
%       .res_history .time_history .regime_history .s_alg_history .s_exp_history
%       .resPsym .residualPhasor
%
%   EXAMPLES:
%     % Descriptor Lyapunov, E(t) periodic, adaptive truncation
%     A = PhasorArray.random(2,2,10) - 1*eye(2);
%     E = PhasorArray.random(2,2,10) + 1*eye(2);
%     Q = PhasorArray.eye(2);
%     [P, info] = lyapG(A, Q, E, autoUpdateh=true, verbose=true);
%     RES = d(E.'*P*E, 2*pi) + A.'*P*E + E.'*P*A + Q;  % continuous residual
%
%   See also: SylvHarmonicGen, PhasorArray/lyap, RicHarmonicKlein

arguments
    pA1                                                        % A (both modes)
    pA2                                                        % Q (Lyapunov) or B (Sylvester)
    pA3                               = []                     % E (Lyapunov) or C (Sylvester)
    pA4                               = []                     % Ea (Sylvester only)
    pA5                               = []                     % Eb (Sylvester only)
    nvp.T                   (1,1) double  = 2*pi
    nvp.h                                 = []            % [] triggers autoUpdateh
    nvp.thresholdResidual   (1,1) double  = 1e-6
    nvp.autoUpdateh         (1,1) logical = true
    nvp.maxh                              = []            % default: h0 * 20
    nvp.stagnationWindow    (1,1) {mustBeInteger, mustBePositive} = 5
    nvp.stagnationRatio     (1,1) double  = 0.05
    nvp.verbose             (1,1) logical = false
    nvp.storeResidualPhasor (1,1) logical = false
    nvp.systemType          {mustBeMember(nvp.systemType,  {'rectangle','square'})}   = 'rectangle'
    nvp.updateMethod        {mustBeMember(nvp.updateMethod,{'adaptive','incremental'})} = 'adaptive'
    nvp.derivativeForm      {mustBeMember(nvp.derivativeForm,{'product','sandwich'})}  = 'product'
    nvp.direction           {mustBeMember(nvp.direction,     {'backward','forward'})}  = 'backward'
end

T                 = nvp.T;
omega             = 2*pi / T;
thresholdResidual = nvp.thresholdResidual;
autoUpdateh       = nvp.autoUpdateh;

%% --- Mode dispatch, validation, and h initialisation ---

isLyapunov = isempty(pA4) && isempty(pA5);

if isLyapunov
    % Generalized Lyapunov mode: d/dt(E.'PE) + A.'PE + E.'PA + Q = 0
    if ~isa(pA1,'PhasorArray'), pA1 = PhasorArray(pA1); end
    if ~isa(pA2,'PhasorArray'), pA2 = PhasorArray(pA2); end
    if size(pA1,1) ~= size(pA1,2)
        error('PhasorArray:lyapG:nonSquareA', 'A must be square (got %dx%d).', size(pA1,1), size(pA1,2))
    end
    if size(pA1,1)~=size(pA2,1) || size(pA1,2)~=size(pA2,2)
        error('PhasorArray:lyapG:dimensionMismatch', 'A (%dx%d) and Q (%dx%d) must be the same size.', size(pA1,1), size(pA1,2), size(pA2,1), size(pA2,2))
    end
    if isempty(pA3), pA3 = PhasorArray(eye(size(pA1,1))); end    % E = [] → identity
    if ~isa(pA3,'PhasorArray'), pA3 = PhasorArray(pA3); end
    if size(pA3,1)~=size(pA1,1) || size(pA3,2)~=size(pA1,2)
        error('PhasorArray:lyapG:dimensionMismatch', 'E (%dx%d) must be the same size as A (%dx%d).', size(pA3,1), size(pA3,2), size(pA1,1), size(pA1,2))
    end
    h = nvp.h;
    if isempty(h)
        h = max([pA1.h, pA2.h, pA3.h]);
        autoUpdateh = true;
    end
    % Remap to generalized Sylvester:
    %   d/dt(Ea X Eb) + A_s X Eb + Ea X B_s + C = 0
    %   with  Ea = E.',  Eb = E,  A_s = A.',  B_s = A,  C = Q
    E  = pA3;
    pA3 = pA2;        % C  = Q
    pA2 = pA1;        % B  = A
    pA1 = pA1.';      % A_s = A.'
    pA4 = E.';       % Ea = E.'
    pA5 = E;         % Eb = E

else
    % Generalized Sylvester mode: d/dt(Ea X Eb) + A X Eb + Ea X B + C = 0
    if ~isa(pA1,'PhasorArray'), pA1 = PhasorArray(pA1); end
    if ~isa(pA2,'PhasorArray'), pA2 = PhasorArray(pA2); end
    if ~isa(pA3,'PhasorArray'), pA3 = PhasorArray(pA3); end
    if size(pA1,1) ~= size(pA1,2)
        error('PhasorArray:lyapG:nonSquareA', 'Sylvester: A must be square (got %dx%d).', size(pA1,1), size(pA1,2))
    end
    if size(pA2,1) ~= size(pA2,2)
        error('PhasorArray:lyapG:nonSquareB', 'Sylvester: B must be square (got %dx%d).', size(pA2,1), size(pA2,2))
    end
    if size(pA1,1) ~= size(pA3,1)
        error('PhasorArray:lyapG:dimensionMismatch', 'Sylvester: rows(A)=%d must equal rows(C)=%d.', size(pA1,1), size(pA3,1))
    end
    if size(pA2,2) ~= size(pA3,2)
        error('PhasorArray:lyapG:dimensionMismatch', 'Sylvester: cols(B)=%d must equal cols(C)=%d.', size(pA2,2), size(pA3,2))
    end
    if isempty(pA4), pA4 = PhasorArray(eye(size(pA1,1))); end    % Ea = [] → identity
    if isempty(pA5), pA5 = PhasorArray(eye(size(pA2,1))); end    % Eb = [] → identity
    if ~isa(pA4,'PhasorArray'), pA4 = PhasorArray(pA4); end
    if ~isa(pA5,'PhasorArray'), pA5 = PhasorArray(pA5); end
    if size(pA4,1)~=size(pA1,1) || size(pA4,2)~=size(pA1,2)
        error('PhasorArray:lyapG:dimensionMismatch', 'Ea (%dx%d) must be the same size as A (%dx%d).', size(pA4,1), size(pA4,2), size(pA1,1), size(pA1,2))
    end
    if size(pA5,1)~=size(pA2,1) || size(pA5,2)~=size(pA2,2)
        error('PhasorArray:lyapG:dimensionMismatch', 'Eb (%dx%d) must be the same size as B (%dx%d).', size(pA5,1), size(pA5,2), size(pA2,1), size(pA2,2))
    end
    h = nvp.h;
    if isempty(h)
        h = max([pA1.h, pA2.h, pA3.h, pA4.h, pA5.h]);
        autoUpdateh = true;
    end
end

% Singularity guard on the mass matrices, evaluated on a time grid.
% If E(t) approaches singularity at some instant, X(t) develops a quasi-pole:
% its harmonics then decay only algebraically and the truncation residual
% stalls regardless of h (regular index-0 descriptor assumption violated).
tGrid  = linspace(0, T, 129);  tGrid(end) = [];
detMinA = minAbsDetOnGrid(pA4, T, tGrid);
detMinB = minAbsDetOnGrid(pA5, T, tGrid);
if detMinA < 1e-3 || detMinB < 1e-3
    warning('PhasorArray:lyapG:nearSingularMass', ...
        ['A mass matrix is (close to) singular at some time instant ' ...
         '(min|det Ea(t)|=%.2e, min|det Eb(t)|=%.2e on a %d-point grid). ' ...
         'lyapG assumes a regular (index-0) descriptor: expect slow ' ...
         '(algebraic) harmonic convergence and a stalling residual.'], ...
        detMinA, detMinB, numel(tGrid))
end

%% --- Single-order solve callback ---

% Closes over the operands so adaptiveHSolve only has to vary h.
solveAtH = @(hh) solveOnce(pA1, pA2, pA3, pA4, pA5, hh, omega, T, nvp);

%% --- Fixed-h: early return ---

if ~autoUpdateh
    [res, resnorm, resrelnorm, resPhasor] = solveAtH(h);
    info = packInfo(3, sprintf('Fixed h=%d.', h), ...
        resrelnorm, resnorm, h, [], [], [], [], {}, [], [], res, resPhasor, nvp);
    info.hForTargetResidual = NaN;   % no refinement, nothing to extrapolate
    info.targetResidual     = NaN;
    return
end

%% --- Adaptive-h refinement ---

n_sys    = size(pA1, 1) * size(pA2, 1);   % nA*nB: size of the vectorised solution
% Spectral width of the harmonic operator (fixed): widest of
% T(Eb.'⊗A), T(B.'⊗Ea), the derivative product T(Eb.'⊗Ea), and C.
h_op     = max([pA1.h + pA5.h, pA2.h + pA4.h, pA4.h + pA5.h, pA3.h]);
isSquare = strcmp(nvp.systemType, 'square');

cfg = struct( ...
    'thresholdResidual', thresholdResidual, ...
    'maxh',              nvp.maxh, ...
    'stagnationWindow',  nvp.stagnationWindow, ...
    'stagnationRatio',   nvp.stagnationRatio, ...
    'updateMethod',      nvp.updateMethod, ...
    'verbose',           nvp.verbose, ...
    'hOp',               h_op, ...
    'hOutFcn',           @(hh) hh*isSquare + (hh + h_op)*~isSquare, ...
    'preamble',          buildPreamble(pA3, n_sys, h_op, isLyapunov, isSquare), ...
    'label',             'h');

[best, trace] = adaptiveHSolve(solveAtH, h, cfg);

res  = best.sol;
info = packInfo(trace.status, trace.statusMsg, best.resrelnorm, best.resnorm, best.h, ...
    trace.h_history, trace.resrel_history, trace.res_history, trace.time_history, ...
    trace.regime_history, trace.s_alg_history, trace.s_exp_history, ...
    best.sol, best.resPhasor, nvp);
% The order a near-zero residual would need, extrapolated at the exit.
info.hForTargetResidual = trace.hForTargetResidual;
info.targetResidual     = trace.targetResidual;

end % lyapG

%% =========================================================================
function [res, resnorm, resrelnorm, resPhasor] = solveOnce(pA1, pA2, pA3, pA4, pA5, h, omega, T, nvp)
%SOLVEONCE  Solve the generalized Sylvester equation at a single harmonic order.
res = PhasorArray(SylvHarmonicGen(pA1, pA2, pA3, pA4, pA5, h, omega, ...
    nvp.systemType, nvp.derivativeForm, nvp.direction));
[resnorm, resrelnorm, resPhasor] = computeResidual(res, pA1, pA2, pA3, pA4, pA5, T, ...
    nvp.derivativeForm, nvp.direction);
end

%% =========================================================================
function s = buildPreamble(pA3, n_sys, h_op, isLyapunov, isSquare)
%BUILDPREAMBLE  Header text printed above the verbose refinement table.
n1 = size(pA3, 1);  n2 = size(pA3, 2);
s = sprintf('\nPhasorArray.lyapG — adaptive h refinement\n');
if isLyapunov
    s = [s sprintf('  Equation:  d/dt(E''·P·E) + A''·P·E + E''·P·A + Q = 0   [gen. Lyapunov, %dx%d]\n', n1, n2)];
    s = [s sprintf('  Each step solves:  [N·T(E''⊗E'') + T(M_A)] · vec(P) = -vec(Q)\n')];
else
    s = [s sprintf('  Equation:  d/dt(Ea·X·Eb) + A·X·Eb + Ea·X·B + C = 0   [gen. Sylvester, %dx%d]\n', n1, n2)];
    s = [s sprintf('  Each step solves:  [N·T(Eb''⊗Ea) + T(M_S)] · vec(X) = -vec(C)\n')];
end
if isSquare
    s = [s sprintf('  where  n = n1·n2 = %d·%d = %d,  h_op = %d,  hOut = h  [square mode]\n', n1, n2, n_sys, h_op)];
    s = [s sprintf('    M  :  [n·(2·h+1)] × [n·(2·h+1)]   (square — border harmonics discarded)\n')];
else
    s = [s sprintf('  where  n = n1·n2 = %d·%d = %d,  h_op = %d,  hOut = h + h_op\n', n1, n2, n_sys, h_op)];
    s = [s sprintf('    M  :  [n·(2·hOut+1)] × [n·(2·h+1)]   (rectangle — full residual L2 coverage)\n')];
end
s = [s sprintf('    X  :  solution truncated to h harmonics\n\n')];
end

%% =========================================================================
function detMin = minAbsDetOnGrid(E, T, tGrid)
%MINABSDETONGRID  Minimum of |det(E(t))| sampled on a time grid.
%   evalTime signature is (obj, T, t): the period must be passed explicitly,
%   otherwise tGrid would be consumed as the period.
Et = evalTime(E, T, tGrid);
detMin = Inf;
for k = 1:numel(tGrid)
    detMin = min(detMin, abs(det(Et(:, :, k))));
end
end

%% =========================================================================
function [resnorm, resrelnorm, resPhasor] = computeResidual(res, pA1, pA2, pA3, pA4, pA5, T, derivativeForm, direction)
%COMPUTERESIDUAL  Evaluate the generalized Sylvester residual and its norms.
%   'product' : RES = ±d/dt(Ea*X*Eb) + A*X*Eb + Ea*X*B + C
%               (derivative on the FULL product — product rule handled by d)
%   'sandwich': RES = ±Ea*(dX/dt)*Eb + A*X*Eb + Ea*X*B + C
%   Sign: + for 'backward', − for 'forward'.
if strcmp(derivativeForm, 'sandwich')
    derivTerm = pA4 * d(res, T) * pA5;
else
    derivTerm = d(pA4*res*pA5, T);
end
if strcmp(direction, 'forward'), derivTerm = -derivTerm; end
resPhasor  = derivTerm + pA1*res*pA5 + pA4*res*pA2 + pA3;
resnorm    = norm(resPhasor.value, 'fro');
Cnorm      = norm(pA3.value, 'fro');
resrelnorm = resnorm / (Cnorm + eps);
end

%% =========================================================================
function info = packInfo(status, statusMsg, resrelnorm, resnorm, h, ...
        h_history, resrel_history, res_history, time_history, ...
        regime_history, s_alg_history, s_exp_history, res, resPhasor, nvp)
%PACKINFO  Build the info struct with all fields always present.
info.status         = status;
info.statusMsg      = statusMsg;
info.resrelnorm     = resrelnorm;
info.resnorm        = resnorm;
info.h              = h;
info.h_history      = h_history;      % [] when autoUpdateh=false
info.resrel_history = resrel_history; % [] when autoUpdateh=false
info.res_history    = res_history;    % [] when autoUpdateh=false
info.time_history   = time_history;   % [] when autoUpdateh=false
info.regime_history = regime_history; % {} when autoUpdateh=false
info.s_alg_history  = s_alg_history;  % [] when incremental or no algebraic detected
info.s_exp_history  = s_exp_history;  % [] when incremental or no exponential detected

if size(res,1) == size(res,2)
    info.resPsym = norm(value(res - res'), 'fro');
else
    info.resPsym = NaN;
end

if nvp.storeResidualPhasor
    info.residualPhasor = resPhasor;
else
    info.residualPhasor = [];
end
end
