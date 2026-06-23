function [res, info] = lyap(o1, o2, o3, options)
%LYAP  Periodic Lyapunov / Sylvester solver for PhasorArray.
%
%   Lyapunov mode:  [P, info] = lyap(A, Q, ...)
%     Solves:  dP/dt + A'(t)P + PA(t) + Q = 0
%     A and Q must be same-size square PhasorArray (or convertible).
%
%   Sylvester mode: [M, info] = lyap(A, B, C, ...)
%     Solves:  dM/dt + A(t)M + MB(t) + C = 0
%     A and B must be square; rows(A)==rows(C), cols(B)==cols(C).
%
%   Descriptor (generalized) mode: pass a mass matrix via name-value 'E'
%   (Lyapunov) or 'Ea'/'Eb' (Sylvester), e.g.  lyap(A, Q, 'E', E)  solves
%   d/dt(E'PE) + A'PE + E'PA + Q = 0, with 'derivativeForm'
%   ('product'|'sandwich'). Delegates to the lyapG engine. E = [] (default)
%   is the standard equation above.
%
%   Options (name-value):
%     T                   Period of the periodic system (default: 2*pi)
%     h                   Fixed harmonic truncation order; [] = auto (default: [])
%     thresholdResidual   Convergence threshold on relative residual (default: 1e-6)
%     autoUpdateh         Enable adaptive h-refinement loop (default: false)
%     maxh                Hard upper bound on h during adaptive loop; [] = h0*20 (default: [])
%     stagnationWindow    Look-back window for stagnation detection (default: 5)
%     stagnationRatio     Min relative improvement to avoid stagnation (default: 0.05)
%     verbose             Print iteration table to console (default: false)
%     storeResidualPhasor Store full residual PhasorArray in info (default: false)
%     systemType          'rectangle' (default: keeps every Cauchy-product output
%                         harmonic, hOut=h+h_op) or 'square' (Galerkin, hOut=h —
%                         border convolution terms neglected, controls cost explosion)
%     updateMethod        'adaptive' (regime-based jump, default) or 'incremental' (h+1)
%     direction           'backward' (default) solves  dP/dt + A'P + PA + Q = 0
%                         (cost-to-go / observability type); 'forward' solves
%                         dP/dt = A'P + PA + Q (covariance type — for the usual
%                         covariance equation dP/dt = AP + PA' + Q, pass A').
%                         Same convention in Sylvester mode (sign of dM/dt).
%
%   Returns:
%     res   PhasorArray solution
%     info  Solver diagnostics struct — all fields always present:
%       .status          0=converged, 1=stagnated, 2=maxh_reached, 3=fixed_h, 4=algebraic_unreachable
%       .statusMsg       Human-readable exit message
%       .resrelnorm      Final relative residual norm  ||res||_F / (||C||_F + eps)
%       .resnorm         Final absolute residual norm
%       .h               Final harmonic order used
%       .h_history       Vector of h values tried  ([] if autoUpdateh=false)
%       .resrel_history  Corresponding relative residuals ([] if autoUpdateh=false)
%       .res_history     Absolute residual history ([] if autoUpdateh=false)
%       .time_history    Step solve time history in seconds ([] if autoUpdateh=false)
%       .regime_history  Regime per step: 'initial'/'exponential'/'algebraic'/'stagnated'
%       .s_alg_history   Algebraic slope estimates ([] if initial regime only)
%       .s_exp_history   Exponential slope estimates ([] if initial regime only)
%       .resPsym         ||P - P'||_F   (NaN if solution is not square)
%       .residualPhasor  Residual as PhasorArray ([] unless storeResidualPhasor=true)
%
%   See also: SylvHarmonic, RicHarmonicKlein

arguments
    o1                                                        % A (both modes)
    o2                                                        % Q (Lyapunov) or B (Sylvester)
    o3                               = []                     % [] → Lyapunov; provided → Sylvester
    options.T                   (1,1) double  = 2*pi
    options.h                                 = []            % [] triggers autoUpdateh
    options.thresholdResidual   (1,1) double  = 1e-6
    options.autoUpdateh         (1,1) logical = false
    options.maxh                              = []            % default: h0 * 20
    options.stagnationWindow    (1,1) {mustBeInteger, mustBePositive} = 5
    options.stagnationRatio     (1,1) double  = 0.05
    options.verbose             (1,1) logical = false
    options.storeResidualPhasor (1,1) logical = false
    options.systemType          {mustBeMember(options.systemType,  {'rectangle','square'})}   = 'rectangle'
    options.updateMethod        {mustBeMember(options.updateMethod,{'adaptive','incremental'})} = 'adaptive'
    options.direction           {mustBeMember(options.direction,   {'backward','forward'})}      = 'backward'
    options.E                                 = []            % descriptor mass (Lyapunov); [] → standard
    options.Ea                                = []            % descriptor left mass (Sylvester)
    options.Eb                                = []            % descriptor right mass (Sylvester)
    options.derivativeForm      {mustBeMember(options.derivativeForm,{'product','sandwich'})}  = 'product'
end

% --- Descriptor (generalized) mode: delegate to the lyapG engine -----------
% E / Ea / Eb supplied via name-value → solve the descriptor equation
%   Lyapunov : d/dt(E'PE) + A'PE + E'PA + Q = 0
%   Sylvester: d/dt(Ea X Eb) + A X Eb + Ea X B + C = 0
% The standard (no E) path falls through to the base solver below.
if ~isempty(options.E) || ~isempty(options.Ea) || ~isempty(options.Eb)
    fwd = {'T',options.T, 'h',options.h, 'thresholdResidual',options.thresholdResidual, ...
           'autoUpdateh',options.autoUpdateh, 'maxh',options.maxh, ...
           'stagnationWindow',options.stagnationWindow, 'stagnationRatio',options.stagnationRatio, ...
           'verbose',options.verbose, 'storeResidualPhasor',options.storeResidualPhasor, ...
           'systemType',options.systemType, 'updateMethod',options.updateMethod, ...
           'derivativeForm',options.derivativeForm, 'direction',options.direction};
    if isempty(o3)
        [res, info] = lyapG(o1, o2, options.E, [], [], fwd{:});   % Lyapunov descriptor
    else
        [res, info] = lyapG(o1, o2, o3, options.Ea, options.Eb, fwd{:}); % Sylvester descriptor
    end
    return
end

T                 = options.T;
omega             = 2*pi / T;
thresholdResidual = options.thresholdResidual;
autoUpdateh       = options.autoUpdateh;

%% --- Mode dispatch, validation, and h initialisation ---

isLyapunov = isempty(o3);

if isLyapunov
    % Lyapunov mode: dP/dt + A'P + PA + Q = 0
    if ~isa(o1,'PhasorArray'), o1 = PhasorArray(o1); end
    if ~isa(o2,'PhasorArray'), o2 = PhasorArray(o2); end
    if size(o1,1) ~= size(o1,2)
        error('PhasorArray:lyap:nonSquareA', 'A must be square (got %dx%d).', size(o1,1), size(o1,2))
    end
    if size(o1,1)~=size(o2,1) || size(o1,2)~=size(o2,2)
        error('PhasorArray:lyap:dimensionMismatch', 'A (%dx%d) and Q (%dx%d) must be the same size.', size(o1,1), size(o1,2), size(o2,1), size(o2,2))
    end
    h = options.h;
    if isempty(h)
        h = max(o1.h, o2.h);
        autoUpdateh = true;
    end
    % Remap to Sylvester: dP/dt + (A').P + P.(A) + Q = 0  →  A'→o1, A→o2, Q→o3
    o3 = o2; o2 = o1; o1 = o1.';

else
    % Sylvester mode: dM/dt + AM + MB + C = 0
    if ~isa(o1,'PhasorArray'), o1 = PhasorArray(o1); end
    if ~isa(o2,'PhasorArray'), o2 = PhasorArray(o2); end
    if ~isa(o3,'PhasorArray'), o3 = PhasorArray(o3); end
    if size(o1,1) ~= size(o1,2)
        error('PhasorArray:lyap:nonSquareA', 'Sylvester: A must be square (got %dx%d).', size(o1,1), size(o1,2))
    end
    if size(o2,1) ~= size(o2,2)
        error('PhasorArray:lyap:nonSquareB', 'Sylvester: B must be square (got %dx%d).', size(o2,1), size(o2,2))
    end
    if size(o1,1) ~= size(o3,1)
        error('PhasorArray:lyap:dimensionMismatch', 'Sylvester: rows(A)=%d must equal rows(C)=%d.', size(o1,1), size(o3,1))
    end
    if size(o2,2) ~= size(o3,2)
        error('PhasorArray:lyap:dimensionMismatch', 'Sylvester: cols(B)=%d must equal cols(C)=%d.', size(o2,2), size(o3,2))
    end
    h = options.h;
    if isempty(h)
        h = max([o1.h, o2.h, o3.h]);
        autoUpdateh = true;
    end
end

%% --- Initial solve ---

t_start  = tic;
t_step   = tic;
res      = PhasorArray(SylvHarmonic(o1, o2, o3, h, omega, options.systemType, options.direction));
dt_step  = toc(t_step);
[resnorm, resrelnorm, resPhasor] = computeResidual(res, o1, o2, o3, T, options.direction);

%% --- Fixed-h: early return ---

if ~autoUpdateh
    info = packInfo(3, sprintf('Fixed h=%d.', h), ...
        resrelnorm, resnorm, h, [], [], [], [], {}, [], [], res, resPhasor, options);
    return
end

%% --- Adaptive-h loop ---

maxh = options.maxh;
if isempty(maxh), maxh = max(h * 20, h + 20); end  % h+20 guards against h=0
capacity = maxh - h + 1;

h_history      = zeros(1, capacity);
resrel_history = zeros(1, capacity);
res_history    = zeros(1, capacity);
time_history   = zeros(1, capacity);
h_history(1)      = h;
resrel_history(1) = resrelnorm;
res_history(1)    = resnorm;
time_history(1)   = dt_step;
regime_history = {'initial'};   % cell — grows with end+1
s_alg_history  = [];            % grows with end+1
s_exp_history  = [];
nIter = 1;

res_best        = res;
resnorm_best    = resnorm;
resrelnorm_best = resrelnorm;
resPhasor_best  = resPhasor;
h_best          = h;

stagnationWindow    = options.stagnationWindow;
stagnationRatio     = options.stagnationRatio;
algebraic_hit_count = 0;

n_sys      = size(o1, 1) * size(o2, 1);  % n1*n2: vectorised solution size
h_op       = max([o1.h, o2.h, o3.h]);   % spectral width of operators (fixed)
isSquare   = strcmp(options.systemType, 'square');
hOut_of    = @(hh) hh * isSquare + (hh + h_op) * ~isSquare;  % hOut as function of h
if options.verbose
    n1 = size(o3, 1);  n2 = size(o3, 2);
    fprintf('\nPhasorArray.lyap — adaptive h refinement\n')
    if isLyapunov
        fprintf('  Equation:  dP/dt + A''·P + P·A + Q = 0   [Lyapunov, %dx%d]\n', n1, n2)
        fprintf('  Each step solves:  E · vec(X) = vec(Q)\n')
    else
        fprintf('  Equation:  dM/dt + A·M + M·B + C = 0   [Sylvester, %dx%d]\n', n1, n2)
        fprintf('  Each step solves:  E · vec(X) = vec(C)\n')
    end
    if isLyapunov, hop_label = 'max(hA,hQ)'; else, hop_label = 'max(hA,hB,hC)'; end
    if isSquare
        fprintf('  where  n = n1·n2 = %d·%d = %d,  h_op = %s = %d,  hOut = h  [square mode]\n', n1, n2, n_sys, hop_label, h_op)
        fprintf('    E  :  [n·(2·h+1)] × [n·(2·h+1)]   (square — border harmonics discarded)\n')
    else
        fprintf('  where  n = n1·n2 = %d·%d = %d,  h_op = %s = %d,  hOut = h + h_op\n', n1, n2, n_sys, hop_label, h_op)
        fprintf('    E  :  [n·(2·hOut+1)] × [n·(2·h+1)]   (rectangle — full residual L2 coverage)\n')
    end
    fprintf('    X  :  solution truncated to h harmonics\n')
    fprintf('\n')
    fprintf('%4s | %4s | %11s | %12s | %-12s | %8s | %9s | %s\n', ...
        'h', 'hOut', 'Res norm', 'Rel res norm', 'Regime', 'Step (s)', 'Total (s)', 'Note')
    fprintf('-----|------|-------------|--------------|--------------|----------|-----------|------\n')
    note0 = '';
    if resrelnorm <= thresholdResidual, note0 = 'converged'; end
    fprintf('%4d | %4d | %11.4e | %12.4e | %-12s | %8.3f | %9.3f | %s\n', ...
        h, hOut_of(h), resnorm, resrelnorm, 'initial', dt_step, toc(t_start), note0)
end

%% --- Check initial-solve convergence before entering the loop ---

if resrelnorm <= thresholdResidual
    status    = 0;
    statusMsg = sprintf('Converged at h=%d (initial solve, resrel=%.2e).', h, resrelnorm);
else
    status    = -1;
    statusMsg = '';
end

while status == -1 && h < maxh
    %% --- Adaptive step selection ---
    regime = 'initial';

    if strcmp(options.updateMethod, 'incremental') || h < h_op * 1.1 || nIter <= 1
        h = h + 1;
    else  % adaptive
        idx_start = find(h_history(1:nIter) >= h_op * 1.1, 1);
        if isempty(idx_start) || idx_start >= nIter
            h = h + 1;
        else
            h1 = h_history(idx_start);  e1 = resrel_history(idx_start);
            h2 = h_history(nIter);      e2 = resrel_history(nIter);

            s_exp = (log(e2+eps) - log(e1+eps)) / (h2 - h1 + eps);
            s_alg = (log(e2+eps) - log(e1+eps)) / (log(h2+eps) - log(h1+eps));

            h_exp = h2 + ceil((log(thresholdResidual+eps) - log(e2+eps)) / (s_exp - eps));
            h_alg = ceil(h2 * (thresholdResidual / (e2+eps))^(1 / (s_alg - eps)));

            s_alg_history(end+1) = s_alg; %#ok<AGROW>
            s_exp_history(end+1) = s_exp; %#ok<AGROW>

            if s_alg < -0.1 && s_alg > -1.5
                deltah = h_alg - h2;
                regime = 'algebraic';
            elseif s_exp < -1e-4
                deltah = h_exp - h2;
                regime = 'exponential';
            else
                deltah = 1;
                regime = 'stagnated';
            end

            deltah = ceil(deltah * 0.8);
            deltah = max(1, deltah);
            deltah = min(deltah, 50);
            deltah = min(deltah, ceil(h * 0.5));
            h      = min(h2 + deltah, maxh);

            % Algebraic early exit: target h is beyond maxh
            if strcmp(regime, 'algebraic') && h_alg > maxh
                algebraic_hit_count = algebraic_hit_count + 1;
                if algebraic_hit_count >= 2
                    status    = 4;
                    statusMsg = sprintf( ...
                        'Algebraic convergence too slow (slope=%.2f). Target h=%d unreachable (maxh=%d). Best: h=%d, resrel=%.2e.', ...
                        s_alg, h_alg, maxh, h_best, resrelnorm_best);
                    res        = res_best;  resnorm    = resnorm_best;
                    resrelnorm = resrelnorm_best;
                    resPhasor  = resPhasor_best;  h = h_best;
                    if options.verbose
                        fprintf('%4d | %4d | %11.4e | %12.4e | %-12s | %8s | %9s | unreachable (slope=%.2f, target h=%d)\n', ...
                            h, hOut_of(h), resnorm, resrelnorm, regime, '-', '-', s_alg, h_alg)
                    end
                    break
                end
            else
                algebraic_hit_count = 0;
            end
        end
    end

    if status == 4, break, end

    %% --- Solve at new h ---
    nIter   = nIter + 1;
    t_step  = tic;
    res     = PhasorArray(SylvHarmonic(o1, o2, o3, h, omega, options.systemType, options.direction));
    dt_step = toc(t_step);
    [resnorm, resrelnorm, resPhasor] = computeResidual(res, o1, o2, o3, T, options.direction);

    h_history(nIter)      = h;
    resrel_history(nIter) = resrelnorm;
    res_history(nIter)    = resnorm;
    time_history(nIter)   = dt_step;
    regime_history{nIter} = regime; %#ok<AGROW>

    if resrelnorm < resrelnorm_best
        res_best        = res;
        resnorm_best    = resnorm;
        resrelnorm_best = resrelnorm;
        resPhasor_best  = resPhasor;
        h_best          = h;
    end

    note = '';

    % Convergence check (inside loop so note appears in the same table row)
    if resrelnorm <= thresholdResidual
        status    = 0;
        statusMsg = sprintf('Converged at h=%d (resrel=%.2e).', h, resrelnorm);
        note = 'converged';
        if options.verbose
            fprintf('%4d | %4d | %11.4e | %12.4e | %-12s | %8.3f | %9.3f | %s\n', ...
                h, hOut_of(h), resnorm, resrelnorm, regime, dt_step, toc(t_start), note)
        end
        break
    end

    % Stagnation check
    if nIter >= stagnationWindow
        window     = resrel_history(nIter - stagnationWindow + 1 : nIter);
        rel_improv = (window(1) - min(window)) / (window(1) + eps);
        if rel_improv < stagnationRatio
            status    = 1;
            statusMsg = sprintf('Stagnated at h=%d (%.1f%% improvement over %d steps). Best: h=%d, resrel=%.2e.', ...
                h, rel_improv*100, stagnationWindow, h_best, resrelnorm_best);
            note = 'stagnated';
        end
    end

    if options.verbose
        fprintf('%4d | %4d | %11.4e | %12.4e | %-12s | %8.3f | %9.3f | %s\n', ...
            h, hOut_of(h), resnorm, resrelnorm, regime, dt_step, toc(t_start), note)
    end

    if status == 1
        res        = res_best;  resnorm    = resnorm_best;
        resrelnorm = resrelnorm_best;
        resPhasor  = resPhasor_best;  h = h_best;
        break
    end
end

%% --- Finalise status (only maxh case remains) ---

if status == -1
    status    = 2;
    statusMsg = sprintf('Reached maxh=%d without convergence. Best: h=%d, resrel=%.2e.', ...
        maxh, h_best, resrelnorm_best);
    res        = res_best;
    resnorm    = resnorm_best;
    resrelnorm = resrelnorm_best;
    resPhasor  = resPhasor_best;
    h          = h_best;
    if options.verbose, fprintf('  → maxh reached. Returning best solution (h=%d).\n', h_best), end
end

h_history      = h_history(1:nIter);
resrel_history = resrel_history(1:nIter);
res_history    = res_history(1:nIter);
time_history   = time_history(1:nIter);
regime_history = regime_history(1:nIter);
info = packInfo(status, statusMsg, resrelnorm, resnorm, h, ...
    h_history, resrel_history, res_history, time_history, ...
    regime_history, s_alg_history, s_exp_history, res, resPhasor, options);

end % lyap

%% =========================================================================
function [resnorm, resrelnorm, resPhasor] = computeResidual(res, o1, o2, o3, T, direction)
%COMPUTERESIDUAL  Evaluate Sylvester residual and its norms.
% 'backward': dX/dt + AX + XB + C = 0 ; 'forward': -dX/dt + AX + XB + C = 0.
if strcmp(direction, 'forward')
    signD = -1;
else
    signD = +1;
end
resPhasor  = signD*res.d(T) + o1*res + res*o2 + o3;
resnorm    = norm(resPhasor.value, 'fro');
Cnorm      = norm(o3.value, 'fro');
resrelnorm = resnorm / (Cnorm + eps);
end

%% =========================================================================
function info = packInfo(status, statusMsg, resrelnorm, resnorm, h, ...
        h_history, resrel_history, res_history, time_history, ...
        regime_history, s_alg_history, s_exp_history, res, resPhasor, options)
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

if options.storeResidualPhasor
    info.residualPhasor = resPhasor;
else
    info.residualPhasor = [];
end
end
