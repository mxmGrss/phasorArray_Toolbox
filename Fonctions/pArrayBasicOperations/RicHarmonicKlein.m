function [K, S, info] = RicHarmonicKlein(A, B, Q, R, K0, T, options)
%RICHARMONICKLEIN  Periodic Riccati solver via Kleinman iteration.
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
%   Full Riccati equation solved:
%     dS/dt + A''*S + S*A - S*B*R^-1*B''*S + Q = 0
%
%   Usage:
%     [K, S, info] = RicHarmonicKlein(A, B, Q, R)
%     [K, S, info] = RicHarmonicKlein(A, B, Q, R, K0)
%     [K, S, info] = RicHarmonicKlein(A, B, Q, R, K0, T)
%     [K, S, info] = RicHarmonicKlein(A, B, Q, R, K0, T, 'maxIter', 50, ...)
%
%   Inputs:
%     A, B          PhasorArray system and input matrices
%     Q, R          State and control weighting (double or PhasorArray)
%     K0            Initial feedback gain; [] → DC LQR fallback (default: [])
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
%     relChangeThreshold    Exit when ‖Sk-Sk-1‖/‖Sk‖ < thr (solution frozen); keep
%                           >> thresholdResidual to avoid premature exit (default: 1e-7)
%                           Prefer stagnation detection for early stopping in most cases.
%     reduceThreshold       Relative threshold for reduce() before lyap; 0 → τ_ric/100
%     warmStartFraction     h_next = max(h0, floor(h_prev*f)); f∈(0,1] (default: 0.7)
%     stagnationWindow      Sliding window for stagnation (default: 5)
%     stagnationRatio       Min relative improvement to avoid stagnation (default: 0.05)
%     verbose               0=silent, 1/2=iteration table with Lyap summary,
%                           3=block-per-iteration with full Lyap trace (default: 0)
%     storeIterates         Store S_history, K_history, lyap_info (default: false)
%
%   Outputs:
%     K     Final feedback gain (PhasorArray)
%     S     Final Riccati solution (PhasorArray)
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
end

omega             = 2*pi / T;   %#ok<NASGU> (used implicitly via lyap)
maxIter           = options.maxIter;
thresholdResidual = options.thresholdResidual;
reduceThresh      = options.reduceThreshold;
if reduceThresh == 0, reduceThresh = thresholdResidual / 100; end

Q = PhasorArray(Q);
R = PhasorArray(R);

outIsReal = all([isreal(A), isreal(B), isreal(Q), isreal(R)]);

%% --- K0 initialisation and stability check ---

[K0, fallbackMsg] = validateOrFallbackK0(K0, A, B, outIsReal,T);

%% --- h initialisation ---

h0 = options.h;
if isempty(h0)
    h0 = max([A.h, B.h, Q.h, R.h]);
    if ~isempty(K0), h0 = max(h0, K0.h); end
    options.autoUpdateh = true;
end
htrunc = h0;

n_sys = size(A, 1);
h_op  = max([A.h, B.h, Q.h, R.h]);

%% --- Verbose header ---

if options.verbose >= 1
    rule = repmat('=', 1, 62);
    fprintf('\n%s\n', rule)
    fprintf('  RicHarmonicKlein  --  Periodic Riccati via Kleinman iteration\n')
    fprintf('%s\n', rule)
    fprintf('  Equation   :  dS/dt + A''*S + S*A - S*B*R^-1*B''*S + Q = 0\n')
    fprintf('  Formulation:  dS/dt + Ak''*S + S*Ak + K''*R*K + Q = 0\n')
    fprintf('                Ak = A - B*K,   K = R^-1*B''*S\n')
    fprintf('%s\n', repmat('-', 1, 62))
    fprintf('  Problem\n')
    fprintf('    State dimension   n        = %d\n',   n_sys)
    fprintf('    Harmonic order    h_op     = %d\n',   h_op)
    fprintf('    Period            T        = %.6g\n', T)
    fprintf('    System type                = %s\n',   options.systemType)
    fprintf('%s\n', repmat('-', 1, 62))
    fprintf('  Convergence\n')
    fprintf('    thresholdResidual          = %.2e   (relative Riccati residual)\n', thresholdResidual)
    fprintf('    relChangeThreshold         = %.2e   (||Sk-Sk-1||/||Sk||, freeze guard)\n', options.relChangeThreshold)
    fprintf('    maxIter                    = %d\n',   options.maxIter)
    fprintf('    stagnationWindow           = %d     steps\n', options.stagnationWindow)
    fprintf('    stagnationRatio            = %.2f   (min relative residual improvement)\n', options.stagnationRatio)
    fprintf('%s\n', repmat('-', 1, 62))
    fprintf('  Lyapunov solver\n')
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
    fprintf('  <a href="matlab:doc RicHarmonicKlein">doc RicHarmonicKlein</a>  —  full notation, options and output reference\n')
    fprintf('\n')
end
if options.verbose == 1 || options.verbose == 2
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
if options.storeIterates
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
    if options.verbose >= 3
        fprintf('\n%s\n  iter %d\n%s\n', repmat('─',1,52), kk, repmat('─',1,52))
    end
    verboseLyap = (options.verbose >= 3);
    [Sk, lyap_info] = lyap(Ak_r, QY_r, ...
        'T',                  T, ...
        'h',                  htrunc, ...
        'autoUpdateh',        options.autoUpdateh, ...
        'maxh',               options.maxh, ...
        'systemType',         options.systemType, ...
        'updateMethod',       options.updateMethod, ...
        'thresholdResidual',  thresholdResidual/100, ...
        'stagnationWindow',   options.stagnationWindow, ...
        'stagnationRatio',    options.stagnationRatio, ...
        'verbose',            verboseLyap);

    if outIsReal, Sk = mreal(Sk); end

    %% --- Derive new gain ---
    Kk_next = Rinv * B.' * Sk;
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

    %% --- Convergence checks (skip on first iteration — no Sk_prev yet) ---
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

    %% --- Riccati residual  (always, from iter 1) ---
    riccati_res = d(Sk, T) + A.' * Sk + Sk * A ...
                  - Sk * B * Rinv * B.' * Sk + Q;
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
function [K0, msg] = validateOrFallbackK0(K0, A, B, outIsReal, T)
%VALIDATEORFALLBACKK0  Ensure K0 stabilises A-B*K0; fall back to DC LQR if needed.
msg = '';

    T_tmp = T;
if ~isempty(K0)
    if ~isa(K0, 'PhasorArray'), K0 = PhasorArray(K0); end
    Ak0 = A - B * K0;
    LL  = findTrueFloquet(Ak0,T_tmp);
    if all(real(LL) <= 0), return, end
    warning('RicHarmonicKlein:unstableK0', ...
        'Provided K0 does not stabilise A-B*K0. Attempting DC LQR fallback.');
    LL
end

% DC LQR fallback
A_dc = A.phas(0);
B_dc = B.phas(0);
Q_dc = PhasorArray(eye(size(A,1))).phas(0);     % placeholder; user Q used in loop
R_dc = eye(size(B,2));
try
    [~, K0_lqr, ~] = icare(A_dc, B_dc, Q_dc, R_dc);
catch
    K0_lqr = [];
end
if isempty(K0_lqr) || any(~isfinite(K0_lqr(:)))
    error('RicHarmonicKlein:noStabilizingK0', ...
        'DC LQR fallback failed. Provide a stabilising K0 manually.')
end
K0 = PhasorArray(K0_lqr);
if outIsReal, K0 = mreal(K0); end

T_tmp = 2*pi;
Ak0 = A - B * K0;
LL  = findTrueFloquet(Ak0,T_tmp);
if any(real(LL) > 0)
    LL
    error('RicHarmonicKlein:noStabilizingK0', ...
        'DC LQR fallback does not stabilise the periodic system. Provide a stabilising K0 manually.')
end
msg = 'DC LQR fallback used for K0.';
end
