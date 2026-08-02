function [best, trace] = adaptiveHSolve(solveAtH, h0, cfg)
%ADAPTIVEHSOLVE  Adaptive harmonic-order refinement driver for harmonic solvers.
%
%   SYNTAX
%     [best, trace] = adaptiveHSolve(solveAtH, h0, cfg)
%
%   DESCRIPTION
%   Repeatedly solves a harmonic equation at increasing truncation order h until
%   the relative residual falls below cfg.thresholdResidual, the residual
%   stagnates, or h reaches cfg.maxh. Shared by PhasorArray/lyap,
%   PhasorArray/lyapG, PhasorArray/mlHmcDivide, PhasorArray/mrHmcDivide and
%   PhasorArray/place, which differ only in the equation being solved (supplied
%   as the solveAtH callback) and in the info struct they finally publish.
%
%   STEP-SIZE STRATEGY
%   The relative residual e(h) of a harmonic truncation typically decays either
%   exponentially, e ~ exp(s_exp*h), or algebraically, e ~ h^s_alg, depending on
%   the smoothness of the underlying periodic solution. From the two samples
%   (h1,e1) and (h2,e2), with h1 the first sample past the operator's spectral
%   width, the two slopes are estimated as
%
%       s_exp = [log(e2) - log(e1)] / (h2 - h1)
%       s_alg = [log(e2) - log(e1)] / [log(h2) - log(h1)]
%
%   and extrapolated to the target residual to give a candidate next order:
%
%       h_exp = h2 + [log(thr) - log(e2)] / s_exp
%       h_alg = h2 * (thr / e2)^(1/s_alg)
%
%   An algebraic regime is declared for -1.5 < s_alg < -0.1 (slow decay,
%   trust the power law); otherwise an exponential regime for s_exp < -1e-4.
%   Failing both, the step falls back to +1. The jump is damped by 0.8 and
%   clamped to [1, min(50, ceil(h/2))] so a single extrapolation cannot
%   overshoot into an intractable problem size. Below h = 1.1*hOp the asymptotic
%   regime has not started and stepping is +1 — unless cfg.maxUnitSteps of them
%   have gone by, which a wide operator with a narrow solution would otherwise
%   never escape.
%
%   INPUTS
%     solveAtH - Function handle, @(h) -> [sol, resnorm, resrelnorm, resPhasor].
%                Closes over the caller's operands. Called once per order.
%     h0       - Initial harmonic order (non-negative integer).
%     cfg      - Configuration struct with fields:
%                  .thresholdResidual  Target relative residual. The callback
%                                      owns the normalisation; every solver
%                                      here divides by the right-hand side,
%                                      never by the solution.
%                  .maxh               Upper bound on h; [] = max(h0*20, h0+20).
%                  .stagnationWindow   Look-back length for stagnation detection,
%                                      which watches the residual and aborts.
%                                      Unrelated to maxUnitSteps below, which
%                                      watches the stepper and pushes on.
%                  .stagnationRatio    Min relative improvement over that window.
%                  .updateMethod       'adaptive' or 'incremental'.
%                  .verbose            Print the iteration table.
%                  .hOp                Spectral width of the harmonic operator.
%                  .maxUnitSteps       Unit steps tolerated before extrapolation
%                                      is forced despite the hOp gate. Default 5.
%                  .hOutFcn            @(h) -> equation order, for the table.
%                  .preamble           Text printed above the table (verbose).
%                  .label              Order name used in messages, e.g. "h".
%
%   OUTPUTS
%     best  - Struct with the best iterate found:
%               .sol .resnorm .resrelnorm .resPhasor .h
%     trace - Struct with the diagnostics the callers publish:
%               .status     0=converged, 1=stagnated, 2=maxh reached,
%                           4=algebraic convergence unreachable
%               .statusMsg  Human-readable exit message
%               .h_history .resrel_history .res_history .time_history
%               .regime_history .s_alg_history .s_exp_history
%
%   Status 3 (fixed h) is never produced here — callers handle the fixed-h path
%   themselves and do not enter this driver.
%
%   EXAMPLES
%     solve = @(hh) solveMyEquation(A, B, hh);
%     cfg   = struct('thresholdResidual', 1e-8, 'maxh', 40, ...
%                    'stagnationWindow', 5, 'stagnationRatio', 0.05, ...
%                    'updateMethod', 'adaptive', 'verbose', true, ...
%                    'hOp', A.h + B.h, 'hOutFcn', @(hh) hh, ...
%                    'preamble', "", 'label', "h");
%     [best, trace] = adaptiveHSolve(solve, 4, cfg);
%
%   See also: PhasorArray/lyap, PhasorArray/lyapG, PhasorArray/mlHmcDivide,
%             PhasorArray/mrHmcDivide, PhasorArray/place, warnIfNotConverged

arguments
    solveAtH (1,1) function_handle
    h0       (1,1) double {mustBeInteger, mustBeNonnegative}
    cfg      (1,1) struct
end

thresholdResidual = cfg.thresholdResidual;
stagnationWindow  = cfg.stagnationWindow;
stagnationRatio   = cfg.stagnationRatio;
hOp               = cfg.hOp;
label             = cfg.label;
if ~isfield(cfg, 'maxUnitSteps') || isempty(cfg.maxUnitSteps)
    cfg.maxUnitSteps = 5;    % unit steps tolerated before forcing extrapolation
end

h    = h0;
maxh = cfg.maxh;
if isempty(maxh), maxh = max(h*20, h + 20); end   % h+20 guards against h = 0

% One iteration raises h by at least 1, so the number of steps is bounded by
% the number of admissible orders.
capacity = max(maxh - h + 1, 1);
h_history      = zeros(1, capacity);
resrel_history = zeros(1, capacity);
res_history    = zeros(1, capacity);
time_history   = zeros(1, capacity);
regime_history = {'initial'};   % cell — grows with end+1
s_alg_history  = [];            % grows with end+1
s_exp_history  = [];

%% --- Initial solve at h0 ---

t_start = tic;
t_step  = tic;
[sol, resnorm, resrelnorm, resPhasor] = solveAtH(h);
dt_step = toc(t_step);

nIter             = 1;
h_history(1)      = h;
resrel_history(1) = resrelnorm;
res_history(1)    = resnorm;
time_history(1)   = dt_step;

best = struct('sol', sol, 'resnorm', resnorm, 'resrelnorm', resrelnorm, ...
              'resPhasor', resPhasor, 'h', h);

algebraic_hit_count = 0;
algebraic_streak_h0 = Inf;      % h at which the current algebraic streak began

if cfg.verbose
    fprintf('%s', cfg.preamble);
    fprintf('%4s | %4s | %11s | %12s | %-12s | %8s | %9s | %s\n', ...
        label, 'hOut', 'Res norm', 'Rel res norm', 'Regime', 'Step (s)', 'Total (s)', 'Note')
    fprintf('-----|------|-------------|--------------|--------------|----------|-----------|------\n')
    note0 = '';
    if resrelnorm <= thresholdResidual, note0 = 'converged'; end
    printRow(h, resnorm, resrelnorm, 'initial', sprintf('%8.3f', dt_step), ...
        sprintf('%9.3f', toc(t_start)), note0);
end

%% --- Check initial-solve convergence before entering the loop ---

if resrelnorm <= thresholdResidual
    status    = 0;
    statusMsg = sprintf('Converged at %s=%d (initial solve, resrel=%.2e).', label, h, resrelnorm);
else
    status    = -1;
    statusMsg = '';
end

%% --- Refinement loop ---

while status == -1 && h < maxh
    %% --- Adaptive step selection ---
    regime = 'initial';

    % The 1.1*hOp gate assumes convergence well above hOp. A wide operator with
    % a narrow solution never opens it and the loop crawls at +1, so it also
    % opens after maxUnitSteps unit steps; the jump stays damped and clamped.
    forceExtrapolation = nIter - 1 >= cfg.maxUnitSteps && ...
            all(strcmp(regime_history(max(1,nIter-cfg.maxUnitSteps+1):nIter), 'initial'));

    if strcmp(cfg.updateMethod, 'incremental') || (h < hOp*1.1 && ~forceExtrapolation) || nIter <= 1
        h = h + 1;
    else
        idx_start = find(h_history(1:nIter) >= hOp*1.1, 1);
        if isempty(idx_start) && forceExtrapolation
            % Nothing reached the asymptotic band; fit on the last maxUnitSteps.
            idx_start = nIter - cfg.maxUnitSteps + 1;
        end
        if isempty(idx_start) || idx_start >= nIter
            % Not enough asymptotic history to fit a slope — fall back to +1.
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

            % Damp and clamp the extrapolated jump (h here is still the old h2).
            deltah = ceil(deltah * 0.8);
            deltah = max(1, deltah);
            deltah = min(deltah, 50);
            deltah = min(deltah, ceil(h * 0.5));
            h      = min(h2 + deltah, maxh);

            % Algebraic early exit: the extrapolated target lies beyond maxh.
            % Two consecutive hits are not enough on their own. A geometric decay
            % looks algebraic at first: on 1/(1+0.95cos) the slope between h=3 and
            % h=5 reads -0.91, squarely inside the algebraic band, and predicts a
            % target of 1e15 -- while the residual actually reaches 1e-15 at
            % h=160. The loop has to climb far enough for the true regime to show
            % (there the slope is -18.7, well outside the band), so the exit also
            % waits for h to have quadrupled since the streak began.
            if strcmp(regime, 'algebraic') && h_alg > maxh
                algebraic_hit_count = algebraic_hit_count + 1;
                if algebraic_hit_count == 1
                    algebraic_streak_h0 = h2;
                end
                if algebraic_hit_count >= 2 && h2 >= 4 * algebraic_streak_h0
                    status    = 4;
                    statusMsg = sprintf( ...
                        'Algebraic convergence too slow (slope=%.2f). Target %s=%d unreachable (max%s=%d). Best: %s=%d, resrel=%.2e.', ...
                        s_alg, label, h_alg, label, maxh, label, best.h, best.resrelnorm);
                    if cfg.verbose
                        printRow(best.h, best.resnorm, best.resrelnorm, regime, '       -', ...
                            '        -', sprintf('unreachable (slope=%.2f, target %s=%d)', s_alg, label, h_alg));
                    end
                    break
                end
            else
                algebraic_hit_count = 0;
                algebraic_streak_h0 = Inf;
            end
        end
    end

    % Progress guard: every branch above is meant to raise h by at least one.
    % Enforcing it here makes a non-terminating loop structurally impossible,
    % including for step rules added later.
    if h <= h_history(nIter)
        h = h_history(nIter) + 1;
    end

    %% --- Solve at the new order ---
    nIter   = nIter + 1;
    t_step  = tic;
    [sol, resnorm, resrelnorm, resPhasor] = solveAtH(h);
    dt_step = toc(t_step);

    h_history(nIter)      = h;
    resrel_history(nIter) = resrelnorm;
    res_history(nIter)    = resnorm;
    time_history(nIter)   = dt_step;
    regime_history{nIter} = regime;

    if resrelnorm < best.resrelnorm
        best = struct('sol', sol, 'resnorm', resnorm, 'resrelnorm', resrelnorm, ...
                      'resPhasor', resPhasor, 'h', h);
    end

    note = '';

    % Convergence check sits inside the loop so its note lands on the same row.
    if resrelnorm <= thresholdResidual
        status    = 0;
        statusMsg = sprintf('Converged at %s=%d (resrel=%.2e).', label, h, resrelnorm);
        note      = 'converged';
        % This iterate is the converged one — publish it even if a marginally
        % smaller residual was seen at a lower order.
        best = struct('sol', sol, 'resnorm', resnorm, 'resrelnorm', resrelnorm, ...
                      'resPhasor', resPhasor, 'h', h);
        if cfg.verbose
            printRow(h, resnorm, resrelnorm, regime, sprintf('%8.3f', dt_step), ...
                sprintf('%9.3f', toc(t_start)), note);
        end
        break
    end

    % Stagnation: too little improvement across the look-back window.
    if nIter >= stagnationWindow
        window     = resrel_history(nIter - stagnationWindow + 1 : nIter);
        rel_improv = (window(1) - min(window)) / (window(1) + eps);
        if rel_improv < stagnationRatio
            status    = 1;
            statusMsg = sprintf('Stagnated at %s=%d (%.1f%% improvement over %d steps). Best: %s=%d, resrel=%.2e.', ...
                label, h, rel_improv*100, stagnationWindow, label, best.h, best.resrelnorm);
            note      = 'stagnated';
        end
    end

    if cfg.verbose
        printRow(h, resnorm, resrelnorm, regime, sprintf('%8.3f', dt_step), ...
            sprintf('%9.3f', toc(t_start)), note);
    end

    if status == 1, break, end
end

%% --- Finalise: only the maxh exit remains ---

if status == -1
    status    = 2;
    statusMsg = sprintf('Reached max%s=%d without convergence. Best: %s=%d, resrel=%.2e.', ...
        label, maxh, label, best.h, best.resrelnorm);
    if cfg.verbose
        fprintf('  → max%s reached. Returning best solution (%s=%d).\n', label, label, best.h)
    end
end

trace.status         = status;
trace.statusMsg      = statusMsg;
trace.h_history      = h_history(1:nIter);
trace.resrel_history = resrel_history(1:nIter);
trace.res_history    = res_history(1:nIter);
trace.time_history   = time_history(1:nIter);
trace.regime_history = regime_history(1:nIter);
trace.s_alg_history  = s_alg_history;
trace.s_exp_history  = s_exp_history;

    function printRow(hh, rn, rrn, reg, stepStr, totalStr, note)
        %PRINTROW  Emit one line of the verbose refinement table.
        fprintf('%4d | %4d | %11.4e | %12.4e | %-12s | %s | %s | %s\n', ...
            hh, cfg.hOutFcn(hh), rn, rrn, reg, stepStr, totalStr, note);
    end
end
