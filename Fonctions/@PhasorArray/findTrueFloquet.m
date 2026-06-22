function [true_mu, true_l, stats, ev_m, ev_l] = findTrueFloquet(o1, T, varg)
% FINDTRUEFLOQUET  Recover true Floquet exponents and multipliers with convergence tracking.
%
%   SYNTAX:
%       [true_mu, true_l, stats]             = findTrueFloquet(o1, T)
%       [true_mu, true_l, stats, ev_m, ev_l] = findTrueFloquet(o1, T, Name, Value)
%
%   INPUTS:
%       o1   - (PhasorArray) Square periodic matrix object.
%       T    - (double)      Period of the system. Default: 2*pi.
%
%   NAME-VALUE:
%       hinit              (int)    Initial harmonic truncation order. Default: o1.h.
%       hmax               (int)    Maximum harmonic truncation order. Default: 10*o1.h.
%       thresholdResidual  (double) Convergence threshold on the Hill trace residual.
%                                   Default: 1e-3.
%       stagnationWindow   (int)    Look-back window for stagnation detection. Default: 5.
%       stagnationRatio    (double) Min relative improvement over the window to avoid
%                                   stagnation (0.05 = 5%). Default: 0.05.
%       verbose            (int)    Verbosity level:
%                                     0 = silent (default)
%                                     1 = header + iteration table
%                                     2 = same as 1 (detail level unified)
%                                     3 = detailed per-iteration trace with separators
%
%   OUTPUTS:
%       true_mu  - (vector) True Floquet exponents (nx elements).
%       true_l   - (vector) True Floquet multipliers exp(T*mu).
%       stats    - (struct) Convergence diagnostics — all fields always present:
%                    .status          : 0=converged, 1=stagnated, 2=hmax_reached.
%                    .statusMsg       : Human-readable exit summary.
%                    .hhistory        : Truncation order at each iteration.
%                    .resTrace        : Hill trace residual at each iteration.
%                    .hBest           : h at which the best residual was achieved.
%                    .resBest         : Best (lowest) trace residual found.
%       ev_m     - (vector) Full HSS spectrum — Floquet exponents, (2h+1)*nx elements.
%       ev_l     - (vector) Full HSS spectrum — Floquet multipliers exp(T*ev_m).
%
%   ALGORITHM:
%       Iteratively increases the Hill/HSS harmonic truncation h from hinit to hmax,
%       recomputing the Floquet exponents at each step and checking the
%       Floquet-Liouville trace identity (Abel 1827, Liouville 1838):
%           sum(mu_true) = tr(A_0)   where A_0 is the dc Fourier harmonic of A(t).
%       This identity follows from det(Phi(T)) = exp(integral_0^T tr(A(t)) dt)
%       applied over one period, giving sum(mu_i) = (1/T)*integral_0^T tr(A(t)) dt = tr(A_0).
%       "Hill" in HmqNEig refers to Hill's method (Hill 1877) for forming the
%       block-Toeplitz harmonic matrix; the trace criterion is Liouville's, not Hill's.
%       Convergence is declared when:
%           |sum(mu) - tr(A0)| / |tr(A0)| < thresholdResidual
%       When tr(A0) ≈ 0 (nilpotent dc part), the absolute criterion is used instead:
%           |sum(mu) - tr(A0)| < thresholdResidual
%       The best solution (lowest residual) is returned on stagnation or hmax exit.
%
%   EXIT CODES (stats.status):
%       0  Converged : residual < thresholdResidual.
%       1  Stagnated : improvement < stagnationRatio over stagnationWindow steps.
%       2  hmax      : reached hmax without convergence or stagnation.
%
%   INTERNALS:
%       findTrueFloquet calls HmqNEig to form the Hill matrix spectrum, then
%       calls findTruelm at each h-step to cluster that spectrum into nx modes.
%       findTruelm is an implementation detail — end users should call findTrueFloquet.
%
%   See also: findTruelm (internal clustering step), HmqNEig.

arguments
    o1 PhasorArray
    T    (1,1) double                              = 2*pi
    varg.hinit             (1,1) {mustBeInteger, mustBePositive} = o1.h
    varg.hmax              (1,1) {mustBeInteger, mustBePositive} = o1.h * 10
    varg.thresholdResidual (1,1) double            = 1e-3
    varg.stagnationWindow  (1,1) {mustBeInteger, mustBePositive} = 5
    varg.stagnationRatio   (1,1) double            = 0.05
    varg.verbose           (1,1) {mustBeInteger, mustBeNonnegative} = 0
    varg.method            {mustBeMember(varg.method, {'robust','kde','gap','dbscan','kmeans'})} = 'robust'
end
nx = o1.size(1);
assert(o1.size(1) == o1.size(2), 'PhasorArray must be square for eigenvalue analysis.');
assert(varg.hmax >= varg.hinit, ...
    'findTrueFloquet: hmax (%d) must be >= hinit (%d).', varg.hmax, varg.hinit);

capacity = varg.hmax - varg.hinit + 1;  % max iterations possible

%% --- Trace of A0 (computed once as convergence reference) ---
%   A(t) = sum_k A_k * exp(j*k*omega*t),  A0 = A_{k=0}  (dc harmonic).
%   Floquet-Liouville trace identity:  sum(mu_true) = tr(A0).
trA0     = trace(o1.phas(0));
useAbsCriterion = abs(trA0) < 1e-12 * nx;   % guard: nilpotent dc part
if useAbsCriterion && varg.verbose >= 1
    warning('findTrueFloquet:zeroTrace', ...
        'tr(A0) ≈ 0 (|tr|=%.2e): switching to absolute trace criterion.', abs(trA0))
end

%% --- Verbose header ---

if varg.verbose >= 1
    rule = repmat('=', 1, 60);
    fprintf('\n%s\n', rule)
    fprintf('  findTrueFloquet  --  Floquet analysis via Hill HSS\n')
    fprintf('%s\n', rule)
    fprintf('  Method     : Hill HSS truncation (HmqNEig) + findTruelm (solver: %s)\n', varg.method)
    if useAbsCriterion
        fprintf('  Criterion  : |sum(mu) - tr(A0)| < thr  [absolute — Liouville, tr(A0)~0]\n')
    else
        fprintf('  Criterion  : |sum(mu) - tr(A0)| / |tr(A0)| < thr  [Floquet-Liouville]\n')
    end
    fprintf('%s\n', repmat('-', 1, 60))
    fprintf('  Problem\n')
    fprintf('    State dimension   nx               = %d\n',   nx)
    fprintf('    Period            T                = %.6g\n', T)
    fprintf('    tr(A0)                             = %.6g%+.6gi\n', real(trA0), imag(trA0))
    fprintf('    hinit                              = %d\n',   varg.hinit)
    fprintf('    hmax                               = %d\n',   varg.hmax)
    fprintf('    thresholdResidual                  = %.2e\n', varg.thresholdResidual)
    fprintf('%s\n', repmat('-', 1, 60))
    fprintf('  Stagnation\n')
    fprintf('    stagnationWindow                   = %d\n',   varg.stagnationWindow)
    fprintf('    stagnationRatio                    = %.2f   (min relative improvement)\n', varg.stagnationRatio)
    fprintf('%s\n', rule)
end
if varg.verbose == 1 || varg.verbose == 2
    if useAbsCriterion
        resLabel = 'traceResAbs';
    else
        resLabel = 'traceRes   ';
    end
    fprintf('  Column legend:\n')
    if useAbsCriterion
        fprintf('    traceResAbs = |sum(mu) - tr(A0)|              (absolute, tr(A0)~0)\n')
    else
        fprintf('    traceRes    = |sum(mu) - tr(A0)| / |tr(A0)|  (Floquet-Liouville)\n')
    end
    fprintf('\n')
    fprintf(' %4s | %3s | %11s | %8s | %9s | %s\n', ...
        'iter', 'h', resLabel, 'step(s)', 'total(s)', 'Note')
    fprintf('------|-----|-------------|----------|-----------|------\n')
end

t_total = tic;

%% --- Pre-allocate history ---

hHist  = zeros(1, capacity);
reHist = zeros(1, capacity);

%% --- Initial computation (h = hinit) ---

ev_m = o1.HmqNEig(varg.hinit, T);
ev_l = exp(T * ev_m);
[true_mu, true_l, stats] = findTruelm(ev_m, T, nx, 'method', varg.method, 'trA0', trA0,'verbose',varg.verbose>1);

resCur   = computeResidual(true_mu, trA0, useAbsCriterion);
h        = varg.hinit;
nIter    = 1;
hHist(1) = h;
reHist(1) = resCur;

% Best-solution tracking
true_mu_best = true_mu;
true_l_best  = true_l;
ev_m_best    = ev_m;
ev_l_best    = ev_l;
resBest      = resCur;
hBest        = h;

note = '';
if resCur <= varg.thresholdResidual, note = 'converged'; end

if varg.verbose == 1 || varg.verbose == 2
    fprintf(' %4d | %3d | %11.4e | %8.3f | %9.3f | %s\n', ...
        1, h, resCur, 0, toc(t_total), note)
elseif varg.verbose >= 3
    fprintf('\n%s\n  iter 1  (h=%d)\n%s\n', repmat('-', 1, 52), h, repmat('-', 1, 52))
    fprintf('  traceRes = %.4e\n', resCur)
end

status = -1;
if resCur <= varg.thresholdResidual, status = 0; end

%% --- Iterative h-refinement ---

while status == -1 && h < varg.hmax
    t_step = tic;
    h      = h + 1;
    nIter  = nIter + 1;

    ev_m = o1.HmqNEig(h, T);
    ev_l = exp(T * ev_m);
    [true_mu, true_l, ~] = findTruelm(ev_m, T, nx, 'method', varg.method, 'trA0', trA0);

    resCur           = computeResidual(true_mu, trA0, useAbsCriterion);
    hHist(nIter)     = h;
    reHist(nIter)    = resCur;

    % Best-solution tracking
    if resCur < resBest
        true_mu_best = true_mu;
        true_l_best  = true_l;
        ev_m_best    = ev_m;
        ev_l_best    = ev_l;
        resBest      = resCur;
        hBest        = h;
    end

    note = '';

    % Convergence
    if resCur <= varg.thresholdResidual
        status = 0;
        note   = 'converged';
    end

    % Stagnation (only after filling window)
    if status == -1 && nIter >= varg.stagnationWindow
        win       = reHist(nIter - varg.stagnationWindow + 1 : nIter);
        relImprov = (win(1) - min(win)) / (win(1) + eps);
        if relImprov < varg.stagnationRatio
            status = 1;
            note   = sprintf('stagnated (%.1f%%)', relImprov * 100);
        end
    end

    if varg.verbose == 1 || varg.verbose == 2
        fprintf(' %4d | %3d | %11.4e | %8.3f | %9.3f | %s\n', ...
            nIter, h, resCur, toc(t_step), toc(t_total), note)
    elseif varg.verbose >= 3
        fprintf('\n%s\n  iter %d  (h=%d)\n%s\n', ...
            repmat('-', 1, 52), nIter, h, repmat('-', 1, 52))
        fprintf('  traceRes = %.4e | step = %.3f s | total = %.3f s | %s\n', ...
            resCur, toc(t_step), toc(t_total), note)
    end
end

%% --- hmax exit (no convergence, no stagnation) ---

if status == -1
    status = 2;
end

%% --- Return best solution on non-converged exits ---

if status ~= 0
    true_mu = true_mu_best;
    true_l  = true_l_best;
    ev_m    = ev_m_best;
    ev_l    = ev_l_best;
end

%% --- Pack stats ---

hHist  = hHist(1:nIter);
reHist = reHist(1:nIter);

stats.status    = status;
stats.statusMsg = buildStatusMsg(status, h, varg.hmax, reHist(end), varg.thresholdResidual, ...
    nIter, varg.stagnationWindow, hBest, resBest);
stats.hhistory  = hHist;
stats.resTrace  = reHist;
stats.hBest     = hBest;
stats.resBest   = resBest;

%% --- Verbose footer ---

if varg.verbose >= 1
    rule = repmat('=', 1, 60);
    fprintf('\n  %s\n', stats.statusMsg)
    fprintf('%s\n\n', rule)
end

end

%% =========================================================================

function res = computeResidual(true_mu, trA0, useAbsCriterion)
% computeResidual  Floquet-Liouville trace residual (relative or absolute).
traceErr = abs(sum(true_mu) - trA0);
if useAbsCriterion
    res = traceErr;
else
    res = traceErr / abs(trA0);
end
end

%% =========================================================================

function msg = buildStatusMsg(status, h, hmax, resLast, thr, nIter, stagWin, hBest, resBest)
% buildStatusMsg  One-line exit summary.
switch status
    case 0
        msg = sprintf('Converged  h=%d | traceRes=%.4e (< %.2e)', h, resLast, thr);
    case 1
        msg = sprintf( ...
            'Stagnated at h=%d (%d iters). Best: h=%d, traceRes=%.4e. Try larger hmax.', ...
            h, nIter, hBest, resBest);
        if nIter < stagWin
            msg = [msg, sprintf(' [window=%d not yet full]', stagWin)];
        end
    case 2
        msg = sprintf( ...
            'hmax reached (h=%d). Best: h=%d, traceRes=%.4e >= %.2e — increase hmax.', ...
            hmax, hBest, resBest, thr);
    otherwise
        msg = sprintf('Unknown status %d at h=%d.', status, h);
end
end
