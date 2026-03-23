function [K, S, info] = RicHarmonicKlein(A, B, Q, R, K0, T, options)
%RICHARMONICKLEIN  Periodic Riccati solver via Kleinman iteration.
%
%   Lyapunov formulation (Kleinman):
%     dS/dt + Ak'·S + S·Ak + Yk + Q = 0,   Ak = A - B·K,  Yk = K'·R·K
%   until  K = R⁻¹·B'·S  satisfies the full Riccati equation:
%     dS/dt + A'·S + S·A - S·B·R⁻¹·B'·S + Q = 0
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
%     relChangeThreshold    Convergence on relative change ‖Sk-Sk-1‖/‖Sk‖ (default: 1e-3)
%     reduceThreshold       Relative threshold for reduce() before lyap; 0 → τ_ric/100
%     warmStartFraction     h_next = max(h0, floor(h_prev*f)); f∈(0,1] (default: 0.7)
%     stagnationWindow      Sliding window for stagnation (default: 5)
%     stagnationRatio       Min relative improvement to avoid stagnation (default: 0.05)
%     verbose               0=silent, 1=outer table only, 2=outer table + compact lyap
%                           summary column, 3=block-per-iteration with full lyap table
%                           (default: 0)
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
    options.relChangeThreshold   (1,1) double            = 1e-3
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

[K0, fallbackMsg] = validateOrFallbackK0(K0, A, B, outIsReal);

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
    fprintf('\nRicHarmonicKlein — Kleinman iteration\n')
    fprintf('  Equation:  dS/dt + A''·S + S·A - S·B·R⁻¹·B''·S + Q = 0\n')
    fprintf('  n = %d,  h_op = %d,  tau_ric = %.2e,  systemType = %s\n', ...
        n_sys, h_op, thresholdResidual, options.systemType)
    if ~isempty(fallbackMsg)
        fprintf('  K0: %s\n', fallbackMsg)
    end
    fprintf('\n')
end
if options.verbose == 1 || options.verbose == 2
    % Continuous table — lyap column width differs per level
    if options.verbose == 1
        lyapColFmt = '%-14s';  lyapColHdr = 'Lyap';
    else
        lyapColFmt = '%-22s';  lyapColHdr = 'Lyap (compact)';
    end
    fprintf([' %4s | %4s | %11s | %11s | ' lyapColFmt ' | %8s | %9s | %s\n'], ...
        'iter', 'h', 'Res Riccati', 'Rel dS', lyapColHdr, 'Step (s)', 'Total (s)', 'Note')
    sep = repmat('-', 1, 14 + 2 + (options.verbose-1)*8);  % +2 for spaces flanking column
    fprintf('------|------|-------------|-------------|%s|----------|-----------|------\n', sep)
end

t_total = tic;

%% --- Pre-allocate history ---

resRic_history    = zeros(1, maxIter);
relChange_history = zeros(1, maxIter);
h_history         = zeros(1, maxIter);
if options.storeIterates
    S_history    = cell(1, maxIter);
    K_history    = cell(1, maxIter);
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
    if kk == 1
        Ak = A - B * Kk;
        Yk = Kk.' * R * Kk;
    else
        Kk = Rinv * B.' * Sk;
        if outIsReal, Kk = mreal(Kk); end
        Ak = A - B * Kk;
        Yk = Kk.' * R * Kk;
    end
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
        'thresholdResidual',  thresholdResidual, ...
        'stagnationWindow',   options.stagnationWindow, ...
        'stagnationRatio',    options.stagnationRatio, ...
        'verbose',            verboseLyap);

    if outIsReal, Sk = mreal(Sk); end

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
        K_history{kk}   = Kk;
        lyap_info_h{kk} = lyap_info;
    end

    %% --- Convergence checks (skip on first iteration — no Sk_prev yet) ---
    note = '';
    if options.verbose == 1
        lyap_str = sprintf('h=%d s=%d', h_used, lyap_info.status);
    elseif options.verbose >= 2
        lyap_str = lyapCompact(lyap_info, h_used);
    else
        lyap_str = '';
    end

    if kk > 1
        %% Riccati residual
        riccati_res = d(Sk, T) + A.' * Sk + Sk * A ...
                      - Sk * B * Rinv * B.' * Sk + Q;
        Qnorm      = norm(Q.value, 'fro');
        resRicnorm = norm(riccati_res.value, 'fro') / (Qnorm + eps);

        %% Relative change in S
        relChange = norm(value(Sk - Sk_prev), 'fro') / (norm(value(Sk), 'fro') + eps);

        resRic_history(kk)    = resRicnorm;
        relChange_history(kk) = relChange;

        %% Best-solution tracking
        if resRicnorm < resRic_best
            resRic_best = resRicnorm;
            S_best      = Sk;
            K_best      = Kk;
            kk_best     = kk;
        end

        %% Convergence
        if resRicnorm < thresholdResidual
            status    = 0;
            statusMsg = sprintf('Converged (resRic) at iter %d: resRic=%.2e.', kk, resRicnorm);
            note      = 'converged (resRic)';
        elseif relChange < options.relChangeThreshold
            status    = 1;
            statusMsg = sprintf('Converged (relChange) at iter %d: relChange=%.2e.', kk, relChange);
            note      = 'converged (relChange)';
        end

        %% Stagnation
        if status == -1 && kk >= options.stagnationWindow
            win       = resRic_history(kk - options.stagnationWindow + 1 : kk);
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
    else
        resRic_history(kk)    = NaN;
        relChange_history(kk) = NaN;
    end

    if options.verbose == 1 || options.verbose == 2
        lFmt = lyapColFmt; %#ok<NODEF>
        if kk == 1
            fprintf([' %4d | %4d | %11s | %11s | ' lFmt ' | %8.3f | %9.3f | %s\n'], ...
                kk, h_used, '—', '—', lyap_str, dt_step, toc(t_total), note)
        else
            fprintf([' %4d | %4d | %11.4e | %11.4e | ' lFmt ' | %8.3f | %9.3f | %s\n'], ...
                kk, h_used, resRicnorm, relChange, lyap_str, dt_step, toc(t_total), note)
        end
    elseif options.verbose >= 3
        if kk == 1
            fprintf('\n  → iter %2d | h=%d | resRic=— | relΔS=— | %.3fs | %s\n', ...
                kk, h_used, dt_step, note)
        else
            fprintf('\n  → iter %2d | h=%d | resRic=%.2e | relΔS=%.2e | %.3fs | %s\n', ...
                kk, h_used, resRicnorm, relChange, dt_step, note)
        end
    end

    if status ~= -1, break, end
    Sk_prev = Sk;
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
info.resRic_history    = resRic_history(1:nIter);
info.relChange_history = relChange_history(1:nIter);
info.h_history         = h_history(1:nIter);
if options.storeIterates
    info.S_history  = S_history(1:nIter);
    info.K_history  = K_history(1:nIter);
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
st_codes = {'cv','st','mh','fh','au'};
st = st_codes{min(li.status + 1, 5)};
if isempty(li.h_history)
    s = sprintf('h=%d %s', h_used, st);
    return
end
h0  = li.h_history(1);
nst = numel(li.h_history);
% dominant non-initial regime
reg = li.regime_history(~strcmp(li.regime_history, 'initial'));
if isempty(reg)
    dom = 'ini';
else
    u = unique(reg);
    [~, idx] = max(cellfun(@(r) sum(strcmp(reg, r)), u));
    dom = u{idx}(1:min(3, numel(u{idx})));
end
if h0 == h_used
    s = sprintf('h=%d %dst %s %s', h_used, nst, dom, st);
else
    s = sprintf('h:%d>%d %dst %s %s', h0, h_used, nst, dom, st);
end
end

%% =========================================================================
function [K0, msg] = validateOrFallbackK0(K0, A, B, outIsReal)
%VALIDATEORFALLBACKK0  Ensure K0 stabilises A-B*K0; fall back to DC LQR if needed.
msg = '';

if ~isempty(K0)
    if ~isa(K0, 'PhasorArray'), K0 = PhasorArray(K0); end
    T_tmp = 2*pi;
    Ak0 = A - B * K0;
    LL  = HmqNEig(Ak0, max(Ak0.h, 1), T_tmp, 'fundamental');
    if all(real(LL) <= 0), return, end
    warning('RicHarmonicKlein:unstableK0', ...
        'Provided K0 does not stabilise A-B*K0. Attempting DC LQR fallback.');
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
LL  = HmqNEig(Ak0, max(Ak0.h, 1), T_tmp, 'fundamental');
if any(real(LL) > 0)
    error('RicHarmonicKlein:noStabilizingK0', ...
        'DC LQR fallback does not stabilise the periodic system. Provide a stabilising K0 manually.')
end
msg = 'DC LQR fallback used for K0.';
end
