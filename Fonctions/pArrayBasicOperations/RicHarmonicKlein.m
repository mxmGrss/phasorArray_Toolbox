function [K, S, info] = RicHarmonicKlein(A, B, Q, R, K0, T, nvp)
%RICHARMONICKLEIN  Periodic Riccati solver via Kleinman iteration.
%
%   Solves the periodic (LTP) algebraic Riccati equation in the harmonic
%   domain, for both optimal CONTROL (LQR, backward) and optimal
%   OBSERVATION (Kalman filter, forward via duality).
%
%   ======================================================================
%   RECIPE 1 — Optimal control (periodic LQR)
%   ======================================================================
%   System:    dx/dt = A(t)*x + B(t)*u,   cost  J = int( x'Qx + u'Ru )dt
%   Call:
%       [K, S, info] = RicHarmonicKlein(A, B, Q, R)
%   Solves the CONTROL PDRE (backward / cost-to-go, the default):
%       dS/dt + A'*S + S*A - S*B*R^-1*B'*S + Q = 0
%   Returns:
%       K : optimal state feedback,  u = -K(t)*x,  K = R^-1*B'*S
%       S : periodic cost-to-go matrix,  J*(x0,t0) = x0'*S(t0)*x0
%   Closed loop A - B*K is Floquet-stable on exit (status 0/1).
%
%   ======================================================================
%   RECIPE 2 — Optimal observation (periodic Kalman filter)
%   ======================================================================
%   System:    dx/dt = A(t)*x + w,   y = C(t)*x + v
%              with process noise cov W = E[ww'] and measurement noise
%              cov V = E[vv'].
%   Call (note the TRANSPOSED inputs and direction='forward'):
%       [Kd, Y, info] = RicHarmonicKlein(A', C', W, V, direction='forward')
%       L = Kd';                       % observer (Kalman) gain
%   Solves the FILTER PDRE (forward / covariance type):
%       dY/dt = A*Y + Y*A' - Y*C'*V^-1*C*Y + W
%   Returns (after transposition):
%       L : Kalman gain for the observer
%           dxhat/dt = A*xhat + L*( y - C*xhat )
%       Y : periodic state-error covariance,  Y = E[(x-xhat)(x-xhat)']
%   The estimator loop A - L*C is Floquet-stable on exit.
%
%   ======================================================================
%
%   Algorithm (pseudocode):
%
%     given K0 stabilising A-B*K0
%     for k = 1, 2, ...
%
%       % --- Kleinman linearisation ---
%       Ak = A - B*Kk                      % closed-loop system
%       Yk = Kk''*R*Kk                      % quadratic cost term
%
%       % --- Inner Lyapunov solve (harmonic, with h-adaptation) ---
%       h  = h0                            % start from warm-start or h0
%       repeat
%           solve:  dS/dt + Ak''*S + S*Ak + Yk + Q = 0  at truncation h
%           if lyap residual < thresholdResidual  ->  accept Sk, exit
%           if stagnated or h >= maxh            ->  accept best Sk, exit
%           else                                 ->  increase h, retry
%       end
%
%       % --- Update gain ---
%       Kk+1 = R^-1 * B'' * Sk
%
%       % --- Outer convergence checks ---
%       Ric Res    = ||dS/dt + A''*S + S*A - S*B*R^-1*B''*S + Q||_F
%       Ric relRes = Ric Res / ||Q||_F
%       sol relChg = ||Kk - Kk-1||_F / ||Kk||_F
%
%       if Ric relRes < thresholdResidual  ->  converged
%       if sol relChg < relChangeThreshold ->  solution frozen (machine eps)
%       if stagnated over stagnationWindow ->  stagnated, return best Sk
%
%       % --- Warm-start h for next iteration ---
%       h0 = max(h0_init, floor(h_used * warmStartFraction))
%
%     end
%
%   Full Riccati equation solved (in solver variables A, B, Q, R):
%     dS/dt + A''*S + S*A - S*B*R^-1*B''*S + Q = 0     (direction='backward', default)
%     dS/dt = A''*S + S*A - S*B*R^-1*B''*S + Q         (direction='forward')
%
%   Usage:
%     [K, S, info] = RicHarmonicKlein(A, B, Q, R)
%     [K, S, info] = RicHarmonicKlein(A, B, Q, R, K0)
%     [K, S, info] = RicHarmonicKlein(A, B, Q, R, K0, T)
%     [K, S, info] = RicHarmonicKlein(A, B, Q, R, K0, T, 'maxIter', 50, ...)
%     [Kd, Y]      = RicHarmonicKlein(A', C', W, V, [], T, ...
%                        direction='forward', skipValidate=false)   % Kalman
%
%   Inputs:
%     A, B          PhasorArray system and input matrices
%                   (observation: pass A' and C' instead — see RECIPE 2)
%     Q, R          State and control weighting (double or PhasorArray)
%                   (observation: noise covariances W and V)
%     K0            Initial feedback gain (default: []). [] → DC LQR fallback,
%                   which always runs when K0 is empty (independent of
%                   skipValidate), so RicHarmonicKlein(A,B,Q,R,[],T) just works.
%     T             Period (default: 2*pi)
%
%   Options (name-value):
%     maxIter               Max Kleinman iterations (default: 100)
%     h                     Fixed h for lyap; [] → auto (default: [])
%     autoUpdateh           Enable adaptive h inside lyap (default: false)
%     maxh                  Hard upper bound on h (default: [])
%     systemType            'rectangle' (default) or 'square'
%     updateMethod          'adaptive' (default) or 'incremental'
%     thresholdResidual     Convergence on relative Riccati residual (default: 1e-6)
%     relChangeThreshold    Exit when ‖Kk-Kk-1‖/‖Kk‖ < thr (solution frozen); keep
%                           >> thresholdResidual to avoid premature exit (default: 1e-7)
%                           Prefer stagnation detection for early stopping in most cases.
%     reduceThreshold       Relative threshold for reduce() before lyap; 0 → τ_ric/100
%     warmStartFraction     h_next = max(h0, floor(h_prev*f)); f∈(0,1] (default: 1.0).
%                           f<1 steps back to look for a smaller h. Measured 1.5x to
%                           4x slower on five problems for a gain that trunc gives
%                           free afterwards: truncating K from h=94 to h=40 costs
%                           0.4 ms, leaves the closed loop at -0.72071 unchanged and
%                           discards 1.1e-06 of its energy.
%     stagnationWindow      Sliding window for stagnation (default: 5)
%     stagnationRatio       Min relative improvement to avoid stagnation (default: 0.05)
%     verbose               0=silent, 1/2=iteration table with Lyap summary,
%                           3=block-per-iteration with full Lyap trace (default: 0)
%     storeIterates         Store S_history, K_history, lyap_info (default: false)
%     skipValidate          skip the Floquet check of a *provided* K0 (default: true).
%                           Keep true when the provided K0 is known to be
%                           stabilizing; set false to also re-check it and fall
%                           back if it fails. An empty K0 always runs the fallback.
%     direction             'backward' (default) or 'forward' — sign of dS/dt.
%                           Use 'backward' for control (RECIPE 1),
%                           'forward' for observation/Kalman (RECIPE 2).
%
%   Outputs:
%     K     Final feedback gain (PhasorArray).
%           Control: u = -K*x. Observation: Kalman gain L = K'.
%     S     Final Riccati solution (PhasorArray).
%           Control: cost-to-go. Observation: error covariance Y.
%     info  Diagnostics struct — see spec-algorithm-ric-harmonic-klein-v2.md
%
%   See also: lyap, SylvHarmonic, RicHarmonicKlein (v1)

arguments
    A    PhasorArray
    B    PhasorArray
    Q                                                                       % double or PhasorArray
    R                                                                       % double or PhasorArray
    K0                                                   = []               % initial gain; [] → DC LQR fallback
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
    nvp.warmStartFraction    (1,1) double            = 1.0             % ∈ (0,1]
    nvp.stagnationWindow     (1,1) {mustBeInteger, mustBePositive} = 5
    nvp.stagnationRatio      (1,1) double            = 0.05
    nvp.verbose              (1,1) {mustBeInteger, mustBeNonnegative} = 0
    nvp.storeIterates        (1,1) logical           = false
    nvp.skipValidate         (1,1) logical           = true
    nvp.direction            {mustBeMember(nvp.direction, {'backward','forward'})} = 'backward'
    nvp.E                                            = []               % descriptor mass; [] → base HARE
    nvp.derivativeForm       {mustBeMember(nvp.derivativeForm, {'product','sandwich'})} = 'product'
end

omega             = 2*pi / T;   %#ok<NASGU> (used implicitly via lyap)
maxIter           = nvp.maxIter;
thresholdResidual = nvp.thresholdResidual;
reduceThresh      = nvp.reduceThreshold;
if reduceThresh == 0, reduceThresh = thresholdResidual / 100; end

Q = PhasorArray(Q);
R = PhasorArray(R);

% Fast path: use base HARE without generalized E.
useE = ~isempty(nvp.E);
if useE
    E = nvp.E;
    if ~isa(E, 'PhasorArray'), E = PhasorArray(E); end
else
    E = [];
end

outIsReal = all([isreal(A), isreal(B), isreal(Q), isreal(R)]);
if useE, outIsReal = outIsReal && isreal(E); end

%% --- K0 initialisation and stability check ---
% An empty K0 always triggers the DC LQR fallback (it means "find me a
% stabilising start"), regardless of skipValidate. skipValidate only skips the
% Floquet check of a *provided* K0.
if ~nvp.skipValidate || isempty(K0)
    [K0, fallbackMsg] = validateOrFallbackK0(K0, A, B, E, outIsReal, T, nvp.direction);
else
    fallbackMsg = [];
    if useE
        warning('RicHarmonicKlein:skipValidate', ...
            'Skipping K0 validation: ensure K0 stabilises E*x''=(A-B*K0)*x or the iteration may diverge.')
    else
        warning('RicHarmonicKlein:skipValidate', ...
            'Skipping K0 validation: ensure K0 stabilises A-B*K0 or the iteration may diverge.')
    end
end
% With skipValidate=true the validation branch is skipped, so a numeric K0
% would never be wrapped — K0.h would crash later. Wrap here unconditionally.
if ~isempty(K0) && ~isa(K0, 'PhasorArray'), K0 = PhasorArray(K0); end
%% --- h initialisation ---

hOperands = [A.h, B.h, Q.h, R.h];
if useE, hOperands(end+1) = E.h; end

h0 = nvp.h;
if isempty(h0)
    h0 = max(hOperands);
    if ~isempty(K0), h0 = max(h0, K0.h); end
    nvp.autoUpdateh = true;
end
htrunc = h0;

n_sys = size(A, 1);
h_op  = max(hOperands);

%% --- Verbose header ---

if nvp.verbose >= 1
    rule = repmat('=', 1, 62);
    fprintf('\n%s\n', rule)
    if useE
        fprintf('  RicHarmonicKlein  --  descriptor GPDRE via Kleinman (GHARE in HSS)\n')
        fprintf('%s\n', rule)
        fprintf('  Equation   :  d/dt(E''*X*E) + A''*X*E + E''*X*A - E''*X*B*R^-1*B''*X*E + Q = 0\n')
        fprintf('  Formulation:  d/dt(E''*S*E) + Ak''*S*E + E''*S*Ak + K''*R*K + Q = 0\n')
        fprintf('                Ak = A - B*K,   K = R^-1*B''*S*E\n')
    else
        fprintf('  RicHarmonicKlein  --  Periodic Riccati via Kleinman iteration\n')
        fprintf('%s\n', rule)
        fprintf('  Equation   :  dS/dt + A''*S + S*A - S*B*R^-1*B''*S + Q = 0\n')
        fprintf('  Formulation:  dS/dt + Ak''*S + S*Ak + K''*R*K + Q = 0\n')
        fprintf('                Ak = A - B*K,   K = R^-1*B''*S\n')
    end
    fprintf('%s\n', repmat('-', 1, 62))
    fprintf('  Problem\n')
    fprintf('    State dimension   n        = %d\n',   n_sys)
    fprintf('    Harmonic order    h_op     = %d\n',   h_op)
    if useE, fprintf('    Mass matrix       h_E      = %d\n', E.h); end
    fprintf('    Period            T        = %.6g\n', T)
    fprintf('    System type                = %s\n',   nvp.systemType)
    fprintf('%s\n', repmat('-', 1, 62))
    fprintf('  Convergence\n')
    fprintf('    thresholdResidual          = %.2e   (relative Riccati residual)\n', thresholdResidual)
    fprintf('    relChangeThreshold         = %.2e   (||Kk-Kk-1||/||Kk||, freeze guard)\n', nvp.relChangeThreshold)
    fprintf('    maxIter                    = %d\n',   nvp.maxIter)
    fprintf('    stagnationWindow           = %d     steps\n', nvp.stagnationWindow)
    fprintf('    stagnationRatio            = %.2f   (min relative residual improvement)\n', nvp.stagnationRatio)
    fprintf('%s\n', repmat('-', 1, 62))
    fprintf('  Lyapunov solver\n')
    fprintf('    h0 (harmonic truncation)   = %d\n',   h0)
    fprintf('    autoUpdateh                = %d\n',   nvp.autoUpdateh)
    fprintf('    warmStartFraction          = %.2f\n', nvp.warmStartFraction)
    if ~isempty(nvp.maxh)
        fprintf('    maxh                       = %d\n', nvp.maxh)
    end
    fprintf('    updateMethod               = %s\n',   nvp.updateMethod)
    fprintf('    reduceThreshold            = %.2e\n', reduceThresh)
    fprintf('%s\n', repmat('-', 1, 62))
    fprintf('  Initial gain K0\n')
    if isempty(fallbackMsg)
        fprintf('    Source                     = user-provided\n')
    else
        fprintf('    Source                     = %s\n', fallbackMsg)
    end
    fprintf('%s\n', rule)
    fprintf('  <a href="matlab:doc RicHarmonicKlein">doc RicHarmonicKlein</a>  —  full notation, nvp and output reference\n')
    fprintf('\n')
end
if nvp.verbose == 1 || nvp.verbose == 2
    % verbose=1 and verbose=2 produce the same table (detail level unified)
    lyapColFmt = '%-19s';  lyapColHdr = 'Lyap';
    fprintf('  Column legend:\n')
    fprintf('    Ric Res    = ||dS/dt + A''*S + S*A - S*B*R^-1*B''*S + Q||_F\n')
    fprintf('    Ric relRes = Ric Res / ||Q||_F\n')
    fprintf('    sol relChg = ||Kk - Kk-1||_F / ||Kk||_F\n')
    fprintf('    Lyap       = h=N->M [K] status   e.g. h=6->61 [35] conv\n')
    fprintf('                 N=h start, M=h final, K=Lyap h-adaptation steps\n')
    fprintf('                 status: conv stgn maxh fail auto\n')
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
if nvp.storeIterates
    S_history    = cell(1, maxIter);
    Kk_history    = cell(1, maxIter);
    lyap_info_h  = cell(1, maxIter);
end

%% --- Kleinman loop ---

Kk   = K0;
Sk   = [];
Sk_prev = [];

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

    %% --- reduce before lyap to gate harmonic explosion ---
    Ak_r = reduce(Ak, 'reduceMethod', 'relative', ...
        'reduceThreshold', reduceThresh, 'exclude0Phasor', false);
    QY_r = reduce(Yk + Q, 'reduceMethod', 'relative', ...
        'reduceThreshold', reduceThresh, 'exclude0Phasor', false);

    %% --- Lyapunov solve ---
    if nvp.verbose >= 3
        fprintf('\n%s\n  iter %d\n%s\n', repmat('─',1,52), kk, repmat('─',1,52))
    end
    verboseLyap = (nvp.verbose >= 3);
    lyapArgs = { ...
        'T',                  T, ...
        'h',                  htrunc, ...
        'autoUpdateh',        nvp.autoUpdateh, ...
        'maxh',               nvp.maxh, ...
        'systemType',         nvp.systemType, ...
        'updateMethod',       nvp.updateMethod, ...
        'thresholdResidual',  thresholdResidual/100, ...
        'stagnationWindow',   nvp.stagnationWindow, ...
        'stagnationRatio',    nvp.stagnationRatio, ...
        'direction',          nvp.direction, ...
        'verbose',            verboseLyap};
    if useE
        % Descriptor inner solve: d/dt(E'SE) + Ak'SE + E'SAk + (Yk+Q) = 0
        [Sk, lyap_info] = lyapG(Ak_r, QY_r, E, lyapArgs{:}, ...
            'derivativeForm', nvp.derivativeForm);
    else
        % Base inner solve keeps the specialised (blkdiag) Toeplitz path.
        [Sk, lyap_info] = lyap(Ak_r, QY_r, lyapArgs{:});
    end

    if outIsReal, Sk = mreal(Sk); end

    %% --- Derive new gain ---
    % Base:       K = R^-1*B'*S
    % Descriptor: K = R^-1*B'*S*E
    if useE
        Kk_next = Rinv * B.' * Sk * E;
    else
        Kk_next = Rinv * B.' * Sk;
    end
    if outIsReal, Kk_next = mreal(Kk_next); end

    %% --- Warm h-start for next iteration ---
    if nvp.autoUpdateh && ~isempty(lyap_info.h_history)
        h_used = lyap_info.h;
        htrunc = max(h0, floor(h_used * nvp.warmStartFraction));
    else
        h_used = htrunc;
    end

    dt_step = toc(t_step);
    h_history(kk) = h_used;

    if nvp.storeIterates
        S_history{kk}   = Sk;
        Kk_history{kk}  = Kk_next;
        lyap_info_h{kk} = lyap_info;
    end

    %% --- Convergence checks (skip on first iteration — no Sk_prev yet) ---
    note = '';
    if nvp.verbose == 1
        lyap_st_words = {'conv','stgn','maxh','fail','auto'};
        lyap_st = lyap_st_words{min(lyap_info.status + 1, 5)};
        lyap_str = sprintf('h=%d %s', h_used, lyap_st);
    elseif nvp.verbose >= 2
        lyap_str = lyapCompact(lyap_info, h_used);
    else
        lyap_str = '';
    end

    %% --- Riccati residual  (always, from iter 1) ---
    % 'backward': dS/dt + A'S + SA - SBR^-1B'S + Q = 0
    % 'forward' : -dS/dt + A'S + SA - SBR^-1B'S + Q = 0
    if useE
        % GHARE residual. 'product' differentiates the whole congruence
        % E'*S*E; 'sandwich' differentiates S alone with the masses outside.
        XE = Sk * E;   EX = E.' * Sk;                 % shared subexpressions
        if strcmp(nvp.derivativeForm, 'sandwich')
            derivTerm = E.' * d(Sk, T) * E;
        else
            derivTerm = d(E.' * XE, T);
        end
        if strcmp(nvp.direction, 'forward'), derivTerm = -derivTerm; end
        riccati_res = derivTerm + A.' * XE + EX * A ...
                      - EX * B * Rinv * B.' * XE + Q;
    else
        if strcmp(nvp.direction, 'forward'), signD = -1; else, signD = +1; end
        riccati_res = signD * d(Sk, T) + A.' * Sk + Sk * A ...
                      - Sk * B * Rinv * B.' * Sk + Q;
    end
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
    elseif relChange < nvp.relChangeThreshold
        status    = 1;
        statusMsg = sprintf('Converged (sol relChg) at iter %d: ||Kk-Kk-1||/||Kk|| = %.2e.', kk, relChange);
        note      = 'converged (sol relChg)';
    end

    %% Stagnation
    if status == -1 && kk >= nvp.stagnationWindow
        win        = resRic_history(kk - nvp.stagnationWindow + 1 : kk);
        rel_improv = (win(1) - min(win)) / (win(1) + eps);
        if rel_improv < nvp.stagnationRatio
            status    = 2;
            statusMsg = sprintf( ...
                'Stagnated at iter %d (%.1f%% improvement over %d steps). Best: iter %d, resRic=%.2e.', ...
                kk, rel_improv*100, nvp.stagnationWindow, kk_best, resRic_best);
            note = 'stagnated';
            Sk = S_best;  Kk = K_best;
        end
    end

    %% --- Verbose row print ---
    if nvp.verbose == 1 || nvp.verbose == 2
        lFmt = lyapColFmt; %#ok<NODEF>
        fprintf([' %4d | %3d | %11.4e | %11.4e | %11.4e | ' lFmt ' | %8.3f | %9.3f | %s\n'], ...
            kk, h_used, resRicAbs, resRicnorm, relChange, lyap_str, dt_step, toc(t_total), note)
    elseif nvp.verbose >= 3
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
    if nvp.verbose >= 1
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
if nvp.storeIterates
    info.S_history  = S_history(1:nIter);
    info.K_history  = Kk_history(1:nIter);
    info.lyap_info  = lyap_info_h(1:nIter);
else
    info.S_history  = {};
    info.K_history  = {};
    info.lyap_info  = {};
end

end % RicHarmonicKlein

%% =========================================================================
function s = lyapCompact(li, h_used)
%LYAPCOMPACT  One-line summary of a lyap info struct for verbose=2 column.
%   Status words: conv=converged, stgn=stagnated, maxh=hit maxh, fail=failed, auto=auto-exit
%   Regime words: first 4 chars of the dominant regime name (e.g. adap, incr, init)
st_words = {'conv','stgn','maxh','fail','auto'};
st = st_words{min(li.status + 1, 5)};
if isempty(li.h_history)
    s = sprintf('h=%d %s', h_used, st);
    return
end
h0  = li.h_history(1);
nstep = numel(li.h_history);
% dominant non-initial regime
reg = li.regime_history(~strcmp(li.regime_history, 'initial'));
if isempty(reg)
    dom = 'init';
else
    u = unique(reg);
    [~, idx] = max(cellfun(@(r) sum(strcmp(reg, r)), u));
    dom = u{idx}(1:min(4, numel(u{idx})));
end
if h0 == h_used
    s = sprintf('h=%d [%d] %s', h_used, nstep, st);
else
    s = sprintf('h=%d->%d [%d] %s', h0, h_used, nstep, st);
end
end

%% =========================================================================
function [K0, msg] = validateOrFallbackK0(K0, A, B, E, outIsReal, T, direction)
%VALIDATEORFALLBACKK0  Ensure K0 stabilises the closed loop; LQR fallback if not.
%   Covers both formulations: x' = (A-B*K0)*x for E = [], E*x' = (A-B*K0)*x
%   otherwise. Two fallback levels, each checked against the TRUE periodic
%   system: DC LQR (ignores the periodic coupling, so it often fails), then
%   truncated HSS LQR on the lifted pair (captures it, converges in h).
msg = '';

if ~isempty(K0)
    if ~isa(K0, 'PhasorArray'), K0 = PhasorArray(K0); end
    Ak0 = A - B * K0;
    LL  = closedLoopFloquet(Ak0, E, T, direction);
    % Lenient with a user-supplied K0 (marginal stability accepted); the
    % fallbacks below are held to strict stability via stabilises().
    if all(real(LL) <= 0), return, end
    if isempty(E)
        warning('RicHarmonicKlein:unstableK0', ...
            ['Provided K0 does not stabilise A-B*K0 (max real Floquet exponent ', ...
             '%.3e). Attempting DC LQR fallback.'], max(real(LL)));
    else
        warning('RicHarmonicKlein:unstableK0', ...
            ['Provided K0 does not stabilise E*x''=(A-B*K0)*x (max real Floquet ', ...
             'exponent %.3e). Attempting DC generalized LQR fallback.'], max(real(LL)));
    end
end

n = size(A, 1); m = size(B, 2);
K0 = [];

% --- Fallback 1: DC LQR (uses only the DC blocks) ------------------------
try
    if isempty(E)
        [~, K0_dc, ~] = icare(A.phas(0), B.phas(0), eye(n), eye(m));
    else
        [~, K0_dc, ~] = icare(A.phas(0), B.phas(0), eye(n), eye(m), [], E.phas(0));
    end
catch
    K0_dc = [];
end
if ~isempty(K0_dc) && all(isfinite(K0_dc(:)))
    Kc = PhasorArray(K0_dc); if outIsReal, Kc = mreal(Kc); end
    if stabilises(A, B, Kc, E, T, direction)
        K0  = Kc;
        msg = fallbackLabel(E, 'DC LQR fallback used for K0.', ...
                               'DC generalized LQR fallback used for K0.');
    end
end

% --- Fallback 2: truncated HSS LQR (captures the periodicity) ------------
if isempty(K0)
    for hh = 2:2:12
        try
            Alift = T_tb(A, hh) - N_tb(n, hh, T);
            Blift = T_tb(B, hh);
            Qlift = eye(n*(2*hh+1));
            Rlift = eye(m*(2*hh+1));
            if isempty(E)
                Kh = lqr(Alift, Blift, Qlift, Rlift);
            else
                [~, Kh, ~] = icare(Alift, Blift, Qlift, Rlift, [], T_tb(E, hh));
            end
            if isempty(Kh) || any(~isfinite(Kh(:))), continue, end
            Kc = PhasorArray.fromTBMatrix(Kh, m, 'n1');
            if outIsReal, Kc = mreal(Kc); end
        catch
            continue
        end
        if stabilises(A, B, Kc, E, T, direction)
            K0  = Kc;
            msg = sprintf(fallbackLabel(E, ...
                    'Truncated HSS LQR fallback (h=%d) used for K0.', ...
                    'Truncated HSS generalized LQR fallback (h=%d) used for K0.'), hh);
            break
        end
    end
end

if isempty(K0)
    error('RicHarmonicKlein:noStabilizingK0', ...
        ['DC and truncated-HSS LQR fallbacks failed to stabilise the periodic ', ...
         'system. Provide a stabilising K0 manually.'])
end
end

%% =========================================================================
function s = fallbackLabel(E, baseMsg, descriptorMsg)
%FALLBACKLABEL  Pick the message matching the active formulation.
if isempty(E), s = baseMsg; else, s = descriptorMsg; end
end

%% =========================================================================
function tf = stabilises(A, B, K0, E, T, direction)
%STABILISES  True if the closed loop is strictly Floquet-stable.
tf = all(real(closedLoopFloquet(A - B*K0, E, T, direction)) < 0);
end

%% =========================================================================
function LL = closedLoopFloquet(Acl, E, T, direction)
%CLOSEDLOOPFLOQUET  Floquet exponents of the closed loop, both formulations.
%   E = []  : x'   = Acl*x
%   E given : E*x' = Acl*x, converted explicitly. The pointwise inversion of
%             E is acceptable here — this is a diagnostic and never enters the
%             harmonic solve — so the PhasorInv advisory is silenced.
if strcmp(direction, 'forward')
    if isempty(E), Acl = Acl.'; else, E = E.'; Acl = Acl.'; end
end
if ~isempty(E)
    warnState = warning('off', 'phasorArray:PhasorInv:useHmcDivide');
    cleaner   = onCleanup(@() warning(warnState));
    Acl = inv(E) * Acl;
end
LL = findTrueFloquet(Acl, T);
end
