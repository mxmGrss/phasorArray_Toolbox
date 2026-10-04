function [best, trace] = adaptiveHSolve(solveAtH, h0, cfg)
%ADAPTIVEHSOLVE Refine the harmonic order until convergence or a stopping guard.
%   [best,trace] = adaptiveHSolve(solveAtH,h0,cfg) calls solveAtH(h), which
%   returns [sol,resnorm,resrelnorm,resPhasor]. The callback defines residual
%   normalization; toolbox solvers normalize by the right-hand side.
%
%   Required cfg fields:
%     thresholdResidual - Relative residual target.
%     maxh              - Order ceiling; [] uses max(20*h0,h0+20).
%     stagnationWindow  - Number of residual samples used to test stagnation.
%     stagnationRatio   - Minimum relative improvement over that window.
%     updateMethod      - 'adaptive' or 'incremental'.
%     hOp               - Spectral width of the harmonic operator.
%     verbose           - Print refinement diagnostics.
%     hOutFcn           - @(h) equation order, for the diagnostic table.
%     preamble, label   - Table heading and order label.
%   Optional cfg fields:
%     maxUnitSteps      - Unit steps before forced extrapolation (default 5).
%     targetResidual    - Target for the final order estimate (default 1e-12).
%
%   Adaptive stepping fits exponential or algebraic residual decay:
%     s_exp = log(e2/e1)/(h2-h1), s_alg = log(e2/e1)/log(h2/h1).
%   Algebraic fits use -1.5 < s_alg < -0.1; exponential fits use s_exp < -1e-4.
%   Extrapolated jumps are damped by 0.8 and bounded by min(50,ceil(h/2)).
%   Otherwise the order increases by one. Fits normally use samples at
%   h >= 1.1*hOp; maxUnitSteps also permits extrapolation below this band.
%   Stagnation stopping requires the entire residual window at h >= 1.1*hOp.
%
%   best contains sol, resnorm, resrelnorm, resPhasor and h. A converged
%   iterate is returned at the target; other exits return the best residual.
%   trace contains status/statusMsg, h_history, res_history, resrel_history,
%   time_history, regime_history, s_alg_history and s_exp_history, together
%   with hForTargetResidual and targetResidual. The order estimate is an
%   extrapolation, not an accuracy guarantee.
%   Status: 0=converged, 1=stagnated, 2=maxh reached, 4=unreachable algebraic
%   target. Fixed-order status 3 is handled by callers outside this driver.
%
%   See also PhasorArray/lyap, PhasorArray/lyapG, PhasorArray/mlHmcDivide,
%            PhasorArray/mrHmcDivide, PhasorArray/place.

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

    % Forced extrapolation controls step size, independently of stopping guards.
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

            % Require two consecutive post-band fits and a fourfold order
            % increase before declaring an algebraic target unreachable.
            if strcmp(regime, 'algebraic') && h_alg > maxh && h2 >= 1.1*hOp
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

    % Test stagnation only after the full window resolves the operator band.
    if nIter >= stagnationWindow && ...
            h_history(nIter - stagnationWindow + 1) >= 1.1*hOp
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

% Order that would reach a near-zero residual, extrapolated from the last two
% samples. Free: the same closed form the stepper already uses, evaluated at a
% different target. Taken at the exit rather than mid-loop, where the fit is
% asymptotic and therefore worth something -- early on it is not, which is why
% the algebraic exit above waits before believing it.
hTarget = NaN;
targetResidual = 1e-12;
if isfield(cfg, 'targetResidual') && ~isempty(cfg.targetResidual)
    targetResidual = cfg.targetResidual;
end
if nIter >= 2
    hA = h_history(nIter-1);      eA = resrel_history(nIter-1);
    hB = h_history(nIter);        eB = resrel_history(nIter);
    if eB <= targetResidual
        hTarget = hB;                                   % already there
    elseif eB > 0 && eB < eA && hB > hA
        sExp = (log(eB) - log(eA)) / (hB - hA);
        sAlg = (log(eB) - log(eA)) / (log(hB + eps) - log(hA + eps));
        if sAlg < -0.1 && sAlg > -1.5
            hTarget = ceil(hB * (targetResidual / eB)^(1 / sAlg));
        elseif sExp < -1e-4
            hTarget = hB + ceil((log(targetResidual) - log(eB)) / sExp);
        end
        if ~isfinite(hTarget) || hTarget > 1e6
            hTarget = Inf;      % the fit blew up: no usable prediction
        end
    end
end
trace.hForTargetResidual = hTarget;
trace.targetResidual     = targetResidual;

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
