function [K, S, info] = RicHarmonicKleinGen(A, B, Q, R, E, K0, T, options)
%RICHARMONICKLEINGEN  Periodic generalized (descriptor) Riccati via Kleinman.
%
%   Solves the GPDRE-LTP for the descriptor system  E(t) x' = A(t) x + B(t) u
%   WITHOUT inverting E(t) and WITHOUT forming dE/dt:
%
%     d/dt(E.'*X*E) + A.'*X*E + E.'*X*A - E.'*X*B*R^-1*B.'*X*E + Q = 0
%
%   The optimal feedback is  u = -K(t) x  with  K = R^-1 * B.' * X * E.
%   With E = I this reduces exactly to RicHarmonicKlein.
%
%   Descriptor Kalman filter via duality (forward + sandwich):
%   The filter GPDRE  E*dY/dt*E.' = A*Y*E.' + E*Y*A.' - E*Y*C.'*V^-1*C*Y*E.' + W
%   of the plant  E(t)*dx/dt = A(t)*x + w,  y = C(t)*x + v  is solved by
%       [Kd, Y] = RicHarmonicKleinGen(A.', C.', W, V, E.', [], T, ...
%                     direction='forward', derivativeForm='sandwich')
%   with observer gain  L = Kd.' = E*Y*C.'*V^-1  (no trretro flip needed:
%   the sign flip of the derivative operator plays the time reversal).
%   This is what KalHarmonicKleinGen wraps.
%
%   Algorithm (Newton-Kleinman, see theory_generalized_lyapunov_sylvester.md §3):
%
%     given K0 stabilising the descriptor loop  E x' = (A - B*K0) x
%     for k = 1, 2, ...
%       Ak   = A - B*Kk                          % closed loop
%       Yk   = Kk.'*R*Kk                          % quadratic cost term
%       solve generalized Lyapunov (lyapG):
%           d/dt(E.'*S*E) + Ak.'*S*E + E.'*S*Ak + Yk + Q = 0
%       Kk+1 = R^-1 * B.' * Sk * E                % gain update
%       outer checks on the FULL GHARE residual (see above)
%     end
%
%   At convergence the Kleinman linearisation recovers the GHARE since
%   B.'*X*E = R*K  →  Ak.'*X*E + E.'*X*Ak + K.'*R*K = A.'*X*E + E.'*X*A - K.'*R*K.
%
%   Usage:
%     [K, S, info] = RicHarmonicKleinGen(A, B, Q, R, E)
%     [K, S, info] = RicHarmonicKleinGen(A, B, Q, R, E, K0, T, 'maxIter', 50, ...)
%
%   Inputs:
%     A, B          PhasorArray system and input matrices
%     Q, R          State and control weighting (double or PhasorArray)
%     E             Mass matrix, PhasorArray or numeric (possibly time-varying).
%                   [] → identity (→ standard RicHarmonicKlein problem).
%     K0            Initial feedback gain; [] → DC generalized LQR fallback
%                   via icare(A0, B0, Q0, R0, [], E0) (default: [])
%     T             Period (default: 2*pi)
%
%   Options (name-value) — identical to RicHarmonicKlein, plus:
%     maxIter, h, autoUpdateh, maxh, systemType, updateMethod,
%     thresholdResidual, relChangeThreshold, reduceThreshold,
%     warmStartFraction, stagnationWindow, stagnationRatio,
%     verbose, storeIterates, skipValidate
%     derivativeForm   'product' (default): GPDRE derivative term is
%                      d/dt(E.'XE) — primal descriptor E·x' = A·x + B·u.
%                      'sandwich': derivative term is E.'·(dX/dt)·E — GPDRE
%                      of the mass-inside form d/dt(E·x) = A·x + B·u (the
%                      adjoint structure; used by KalHarmonicKleinGen for
%                      the dual filtering problem). Coincide for constant E.
%     direction        'backward' (default) or 'forward' — sign of the GPDRE
%                      derivative term, propagated to lyapG, to the outer
%                      residual, and to the K0 Floquet check (in forward mode
%                      the physical loop is the TRANSPOSED descriptor
%                      (Acl.', E.') — Floquet exponents differ in LTP).
%                      Orthogonal to derivativeForm; the descriptor Kalman
%                      dual is (forward, sandwich).
%
%   Outputs:
%     K     Final feedback gain (PhasorArray), K = R^-1*B.'*S*E
%     S     Final GPDRE solution X (PhasorArray)
%     info  Diagnostics struct — same fields as RicHarmonicKlein
%
%   See also: RicHarmonicKlein, PhasorArray/lyapG, SylvHarmonicGen

arguments
    A    PhasorArray
    B    PhasorArray
    Q                                                                       % double or PhasorArray
    R                                                                       % double or PhasorArray
    E                                                    = []               % mass matrix; [] → identity
    K0                                                   = []               % initial gain; [] → DC fallback
    T                            (1,1) double            = 2*pi             % period
    options.maxIter              (1,1) {mustBeInteger, mustBePositive} = 100
    options.h                                            = []               % [] → auto
    options.autoUpdateh          (1,1) logical           = false
    options.maxh                                         = []               % [] → h0*20
    options.systemType           {mustBeMember(options.systemType,    {'rectangle','square'})}    = 'rectangle'
    options.updateMethod         {mustBeMember(options.updateMethod,  {'adaptive','incremental'})} = 'adaptive'
    options.thresholdResidual    (1,1) double            = 1e-6
    options.relChangeThreshold   (1,1) double            = 1e-14
    options.reduceThreshold      (1,1) double            = 0               % 0 → τ_ric/100
    options.warmStartFraction    (1,1) double            = 0.7             % ∈ (0,1]
    options.stagnationWindow     (1,1) {mustBeInteger, mustBePositive} = 5
    options.stagnationRatio      (1,1) double            = 0.05
    options.verbose              (1,1) {mustBeInteger, mustBeNonnegative} = 0
    options.storeIterates        (1,1) logical           = false
    options.skipValidate         (1,1) logical           = true
    options.derivativeForm       {mustBeMember(options.derivativeForm, {'product','sandwich'})} = 'product'
    options.direction            {mustBeMember(options.direction,      {'backward','forward'})} = 'backward'
end

maxIter           = options.maxIter;
thresholdResidual = options.thresholdResidual;
reduceThresh      = options.reduceThreshold;
if reduceThresh == 0, reduceThresh = thresholdResidual / 100; end

Q = PhasorArray(Q);
R = PhasorArray(R);
if isempty(E), E = PhasorArray(eye(size(A, 1))); end
if ~isa(E, 'PhasorArray'), E = PhasorArray(E); end

outIsReal = all([isreal(A), isreal(B), isreal(Q), isreal(R), isreal(E)]);

%% --- K0 initialisation and stability check ---
if ~options.skipValidate
    [K0, fallbackMsg] = validateOrFallbackK0(K0, A, B, E, outIsReal, T, options.direction);
else
    fallbackMsg = [];
    if isempty(K0)
        % Fallback still needed even when validation is skipped.
        [K0, fallbackMsg] = validateOrFallbackK0(K0, A, B, E, outIsReal, T, options.direction);
    else
        warning('RicHarmonicKleinGen:skipValidate', ...
            'Skipping K0 validation: ensure K0 stabilises E*x''=(A-B*K0)*x or the iteration may diverge.')
    end
end
if ~isa(K0, 'PhasorArray'), K0 = PhasorArray(K0); end

%% --- h initialisation ---

h0 = options.h;
if isempty(h0)
    h0 = max([A.h, B.h, Q.h, R.h, E.h]);
    if ~isempty(K0), h0 = max(h0, K0.h); end
    options.autoUpdateh = true;
end
htrunc = h0;

n_sys = size(A, 1);
h_op  = max([A.h, B.h, Q.h, R.h, E.h]);

%% --- Verbose header ---

if options.verbose >= 1
    rule = repmat('=', 1, 62);
    fprintf('\n%s\n', rule)
    fprintf('  RicHarmonicKleinGen  --  descriptor GPDRE via Kleinman (GHARE in HSS)\n')
    fprintf('%s\n', rule)
    fprintf('  Equation   :  d/dt(E''*X*E) + A''*X*E + E''*X*A - E''*X*B*R^-1*B''*X*E + Q = 0\n')
    fprintf('  Formulation:  d/dt(E''*S*E) + Ak''*S*E + E''*S*Ak + K''*R*K + Q = 0\n')
    fprintf('                Ak = A - B*K,   K = R^-1*B''*S*E\n')
    fprintf('%s\n', repmat('-', 1, 62))
    fprintf('  Problem\n')
    fprintf('    State dimension   n        = %d\n',   n_sys)
    fprintf('    Harmonic order    h_op     = %d\n',   h_op)
    fprintf('    Mass matrix       h_E      = %d\n',   E.h)
    fprintf('    Period            T        = %.6g\n', T)
    fprintf('    System type                = %s\n',   options.systemType)
    fprintf('%s\n', repmat('-', 1, 62))
    fprintf('  Convergence\n')
    fprintf('    thresholdResidual          = %.2e   (relative GHARE residual)\n', thresholdResidual)
    fprintf('    relChangeThreshold         = %.2e   (||Kk-Kk-1||/||Kk||, freeze guard)\n', options.relChangeThreshold)
    fprintf('    maxIter                    = %d\n',   options.maxIter)
    fprintf('    stagnationWindow           = %d     steps\n', options.stagnationWindow)
    fprintf('    stagnationRatio            = %.2f   (min relative residual improvement)\n', options.stagnationRatio)
    fprintf('%s\n', repmat('-', 1, 62))
    fprintf('  Lyapunov solver (lyapG)\n')
    fprintf('    h0 (harmonic truncation)   = %d\n',   h0)
    fprintf('    autoUpdateh                = %d\n',   options.autoUpdateh)
    fprintf('    warmStartFraction          = %.2f\n', options.warmStartFraction)
    if ~isempty(options.maxh)
        fprintf('    maxh                       = %d\n', options.maxh)
    end
    fprintf('    updateMethod               = %s\n',   options.updateMethod)
    fprintf('    reduceThreshold            = %.2e\n', reduceThresh)
    fprintf('%s\n', repmat('-', 1, 62))
    fprintf('  Initial gain K0\n')
    if isempty(fallbackMsg)
        fprintf('    Source                     = user-provided\n')
    else
        fprintf('    Source                     = %s\n', fallbackMsg)
    end
    fprintf('%s\n', rule)
    fprintf('  <a href="matlab:doc RicHarmonicKleinGen">doc RicHarmonicKleinGen</a>  —  full notation, options and output reference\n')
    fprintf('\n')
end
if options.verbose == 1 || options.verbose == 2
    lyapColFmt = '%-19s';  lyapColHdr = 'LyapG';
    fprintf('  Column legend:\n')
    fprintf('    Ric Res    = ||d/dt(E''*S*E) + A''*S*E + E''*S*A - E''*S*B*R^-1*B''*S*E + Q||_F\n')
    fprintf('    Ric relRes = Ric Res / ||Q||_F\n')
    fprintf('    sol relChg = ||Kk - Kk-1||_F / ||Kk||_F\n')
    fprintf('    LyapG      = h=N->M [K] status   e.g. h=6->61 [35] conv\n')
    fprintf([' %4s | %3s | %11s | %11s | %11s | ' lyapColFmt ' | %8s | %9s | %s\n'], ...
        'iter', 'h', 'Ric Res', 'Ric relRes', 'sol relChg', lyapColHdr, 'step (s)', 'total (s)', 'Note')
    fprintf('------|-----|-------------|-------------|-------------|---------------------|----------|-----------|------\n')
end

t_total = tic;

%% --- Pre-allocate history ---

resRicAbs_history = zeros(1, maxIter);
resRic_history    = zeros(1, maxIter);
relChange_history = zeros(1, maxIter);
h_history         = zeros(1, maxIter);
if options.storeIterates
    S_history    = cell(1, maxIter);
    Kk_history   = cell(1, maxIter);
    lyap_info_h  = cell(1, maxIter);
end

%% --- Kleinman loop ---

Kk = K0;

resRicnorm = NaN;
relChange  = NaN;
S_best     = [];
K_best     = K0;
resRic_best = Inf;
kk_best     = 1;
status      = -1;
statusMsg   = '';

Rinv = inv(R);

for kk = 1:maxIter
    t_step = tic;

    %% --- Compute Ak = A - B*K,  Yk = K'*R*K ---
    Ak = A - B * Kk;
    Yk = Kk.' * R * Kk;
    if outIsReal
        Ak = mreal(Ak);
        Yk = mreal(Yk);
    end

    %% --- reduce before lyapG to gate harmonic explosion ---
    Ak_r = reduce(Ak, 'reduceMethod', 'relative', ...
        'reduceThreshold', reduceThresh, 'exclude0Phasor', false);
    QY_r = reduce(Yk + Q, 'reduceMethod', 'relative', ...
        'reduceThreshold', reduceThresh, 'exclude0Phasor', false);

    %% --- Generalized Lyapunov solve ---
    if options.verbose >= 3
        fprintf('\n%s\n  iter %d\n%s\n', repmat('─',1,52), kk, repmat('─',1,52))
    end
    verboseLyap = (options.verbose >= 3);
    [Sk, lyap_info] = lyapG(Ak_r, QY_r, E, ...
        'T',                  T, ...
        'h',                  htrunc, ...
        'autoUpdateh',        options.autoUpdateh, ...
        'maxh',               options.maxh, ...
        'systemType',         options.systemType, ...
        'updateMethod',       options.updateMethod, ...
        'thresholdResidual',  thresholdResidual/100, ...
        'stagnationWindow',   options.stagnationWindow, ...
        'stagnationRatio',    options.stagnationRatio, ...
        'derivativeForm',     options.derivativeForm, ...
        'direction',          options.direction, ...
        'verbose',            verboseLyap);

    if outIsReal, Sk = mreal(Sk); end

    %% --- Derive new gain:  K = R^-1 * B' * S * E ---
    Kk_next = Rinv * B.' * Sk * E;
    if outIsReal, Kk_next = mreal(Kk_next); end

    %% --- Warm h-start for next iteration ---
    if options.autoUpdateh && ~isempty(lyap_info.h_history)
        h_used = lyap_info.h;
        htrunc = max(h0, floor(h_used * options.warmStartFraction));
    else
        h_used = htrunc;
    end

    dt_step = toc(t_step);
    h_history(kk) = h_used;

    if options.storeIterates
        S_history{kk}   = Sk;
        Kk_history{kk}  = Kk_next;
        lyap_info_h{kk} = lyap_info;
    end

    note = '';
    if options.verbose == 1
        lyap_st_words = {'conv','stgn','maxh','fail','auto'};
        lyap_st = lyap_st_words{min(lyap_info.status + 1, 5)};
        lyap_str = sprintf('h=%d %s', h_used, lyap_st);
    elseif options.verbose >= 2
        lyap_str = lyapCompact(lyap_info, h_used);
    else
        lyap_str = '';
    end

    %% --- Full GHARE residual  (always, from iter 1) ---
    XE = Sk * E;   EX = E.' * Sk;                 % shared subexpressions
    if strcmp(options.derivativeForm, 'sandwich')
        derivTerm = E.' * d(Sk, T) * E;
    else
        derivTerm = d(E.' * XE, T);
    end
    % 'forward' flips the sign of the derivative term (covariance type)
    if strcmp(options.direction, 'forward'), derivTerm = -derivTerm; end
    riccati_res = derivTerm + A.' * XE + EX * A ...
                  - EX * B * Rinv * B.' * XE + Q;
    Qnorm      = norm(Q.value, 'fro');
    resRicAbs  = norm(riccati_res.value, 'fro');
    resRicnorm = resRicAbs / (Qnorm + eps);

    resRicAbs_history(kk) = resRicAbs;
    resRic_history(kk)    = resRicnorm;

    %% Best-solution tracking  (from iter 1)
    if resRicnorm < resRic_best || kk == 1
        resRic_best = resRicnorm;
        S_best      = Sk;
        K_best      = Kk_next;
        kk_best     = kk;
    end

    %% --- Relative change + convergence checks ---
    relChange = norm(value(Kk_next - Kk), 'fro') / (norm(value(Kk_next), 'fro') + eps);
    relChange_history(kk) = relChange;

    %% Convergence
    if resRicnorm < thresholdResidual
        status    = 0;
        statusMsg = sprintf('Converged (Ric relRes) at iter %d: Ric relRes = %.2e.', kk, resRicnorm);
        note      = 'converged (Ric Res)';
    elseif relChange < options.relChangeThreshold
        status    = 1;
        statusMsg = sprintf('Converged (sol relChg) at iter %d: ||Kk-Kk-1||/||Kk|| = %.2e.', kk, relChange);
        note      = 'converged (sol relChg)';
    end

    %% Stagnation
    if status == -1 && kk >= options.stagnationWindow
        win        = resRic_history(kk - options.stagnationWindow + 1 : kk);
        rel_improv = (win(1) - min(win)) / (win(1) + eps);
        if rel_improv < options.stagnationRatio
            status    = 2;
            statusMsg = sprintf( ...
                'Stagnated at iter %d (%.1f%% improvement over %d steps). Best: iter %d, resRic=%.2e.', ...
                kk, rel_improv*100, options.stagnationWindow, kk_best, resRic_best);
            note = 'stagnated';
            Sk = S_best;  Kk = K_best;
        end
    end

    %% --- Verbose row print ---
    if options.verbose == 1 || options.verbose == 2
        lFmt = lyapColFmt; %#ok<NODEF>
        fprintf([' %4d | %3d | %11.4e | %11.4e | %11.4e | ' lFmt ' | %8.3f | %9.3f | %s\n'], ...
            kk, h_used, resRicAbs, resRicnorm, relChange, lyap_str, dt_step, toc(t_total), note)
    elseif options.verbose >= 3
        fprintf('\n  -> iter %2d | h=%d | Ric Res=%.2e | Ric relRes=%.2e | sol relChg=%.2e | %.3fs | %s\n', ...
            kk, h_used, resRicAbs, resRicnorm, relChange, dt_step, note)
    end

    if status ~= -1, break, end
    Kk = Kk_next;
end

%% --- Finalise ---

if status == -1
    status    = 3;
    statusMsg = sprintf('Reached maxIter=%d. Best: iter %d, resRic=%.2e.', ...
        maxIter, kk_best, resRic_best);
    Sk = S_best;  Kk = K_best;
    if options.verbose >= 1
        fprintf('  → maxIter reached. Returning best solution (iter %d).\n', kk_best)
    end
end

K = Kk;
S = Sk;

%% --- Pack info ---

nIter = kk;
info.status            = status;
info.statusMsg         = statusMsg;
info.resRicnorm        = resRicnorm;
info.relChange         = relChange;
info.iter              = nIter;
info.h                 = h_history(nIter);
info.resRicAbs_history = resRicAbs_history(1:nIter);
info.resRic_history    = resRic_history(1:nIter);
info.relChange_history = relChange_history(1:nIter);
info.h_history         = h_history(1:nIter);
if options.storeIterates
    info.S_history  = S_history(1:nIter);
    info.K_history  = Kk_history(1:nIter);
    info.lyap_info  = lyap_info_h(1:nIter);
else
    info.S_history  = {};
    info.K_history  = {};
    info.lyap_info  = {};
end

end % RicHarmonicKleinGen

%% =========================================================================
function s = lyapCompact(li, h_used)
%LYAPCOMPACT  One-line summary of a lyapG info struct for verbose=2 column.
st_words = {'conv','stgn','maxh','fail','auto'};
st = st_words{min(li.status + 1, 5)};
if isempty(li.h_history)
    s = sprintf('h=%d %s', h_used, st);
    return
end
h0  = li.h_history(1);
nstep = numel(li.h_history);
reg = li.regime_history(~strcmp(li.regime_history, 'initial'));
if isempty(reg)
    dom = 'init'; %#ok<NASGU>
else
    u = unique(reg);
    [~, idx] = max(cellfun(@(r) sum(strcmp(reg, r)), u));
    dom = u{idx}(1:min(4, numel(u{idx}))); %#ok<NASGU>
end
if h0 == h_used
    s = sprintf('h=%d [%d] %s', h_used, nstep, st);
else
    s = sprintf('h=%d->%d [%d] %s', h0, h_used, nstep, st);
end
end

%% =========================================================================
function [K0, msg] = validateOrFallbackK0(K0, A, B, E, outIsReal, T, direction)
%VALIDATEORFALLBACKK0  Ensure K0 stabilises E*x'=(A-B*K0)*x; DC fallback if not.
%   Floquet check requires the explicit closed loop inv(E)*(A-B*K0): the
%   pointwise inversion is acceptable HERE (diagnostic only, never used in
%   the solve itself).
%   'forward' (dual / Kalman use): the physical loop is the TRANSPOSED
%   descriptor (Acl.', E.') — e.g. the estimator E*e' = (A-L*C)*e — and in
%   LTP the Floquet exponents of (M,E) and (M.',E.') differ, so transpose
%   before checking.
msg = '';

if ~isempty(K0)
    if ~isa(K0, 'PhasorArray'), K0 = PhasorArray(K0); end
    Ak0 = A - B * K0;
    if strcmp(direction, 'forward'), LL = floquetDescriptor(Ak0.', E.', T);
    else,                            LL = floquetDescriptor(Ak0,   E,   T);
    end
    if all(real(LL) <= 0), return, end
    warning('RicHarmonicKleinGen:unstableK0', ...
        'Provided K0 does not stabilise E*x''=(A-B*K0)*x. Attempting DC generalized LQR fallback.');
end

% DC generalized LQR fallback: icare supports a mass matrix E directly and
% returns the descriptor gain K = R^-1*B'*X*E.
A_dc = A.phas(0);
B_dc = B.phas(0);
E_dc = E.phas(0);
Q_dc = eye(size(A, 1));
R_dc = eye(size(B, 2));
try
    [~, K0_lqr, ~] = icare(A_dc, B_dc, Q_dc, R_dc, [], E_dc);
catch
    K0_lqr = [];
end
if isempty(K0_lqr) || any(~isfinite(K0_lqr(:)))
    error('RicHarmonicKleinGen:noStabilizingK0', ...
        'DC generalized LQR fallback failed. Provide a stabilising K0 manually.')
end
K0 = PhasorArray(K0_lqr);
if outIsReal, K0 = mreal(K0); end

Ak0 = A - B * K0;
if strcmp(direction, 'forward'), LL = floquetDescriptor(Ak0.', E.', T);
else,                            LL = floquetDescriptor(Ak0,   E,   T);
end
if any(real(LL) > 0)
    disp(LL)
    error('RicHarmonicKleinGen:noStabilizingK0', ...
        'DC generalized LQR fallback does not stabilise the periodic descriptor system. Provide a stabilising K0 manually.')
end
msg = 'DC generalized LQR fallback used for K0.';
end

%% =========================================================================
function LL = floquetDescriptor(Acl, E, T)
%FLOQUETDESCRIPTOR  Floquet exponents of E*x' = Acl*x via explicit conversion.
%   Pointwise inversion of E is fine for this diagnostic (never enters the
%   harmonic solve); the PhasorInv advisory warning is therefore silenced.
warnState = warning('off', 'phasorArray:PhasorInv:useHmcDivide');
cleaner   = onCleanup(@() warning(warnState));
Aexpl = inv(E) * Acl;
LL = findTrueFloquet(Aexpl, T);
end
