function [r, info] = mlHmcDivide(A, B, options)
%MLHMCDIVIDE  Harmonic left-division: solves A(t)*X(t) = B(t) in the harmonic domain.
%
%   [X, info] = mlHmcDivide(A, B)
%     Solves A(t)*X(t) = B(t) via overdetermined block-Toeplitz least-squares.
%     Equivalent to lyap(A, 0, -B, "T", Inf): as T→∞, ω→0 and the differential
%     term in the Sylvester equation vanishes, leaving A(t)*X(t) = B(t).
%
%   Three guard cases are handled before the main solve:
%     G1 — A scalar constant (h=0, 1×1): X = B / a
%     G2 — A scalar periodic, B matrix: recurse on each B(i,j)
%     G3 — A matrix, B scalar: error (ill-defined)
%
%   Options (name-value):
%     h                   Fixed hIn; [] = max(A.h, B.h) + A.h (default: [])
%     thresholdResidual   Convergence threshold on relative residual (default: 5e-4)
%     autoUpdateh         Enable adaptive hIn-refinement loop (default: false)
%     maxh                Hard upper bound on hIn; [] = max(hIn*20, hIn+20) (default: [])
%     stagnationWindow    Look-back window for stagnation detection (default: 15)
%     stagnationRatio     Min relative improvement to avoid stagnation (default: 0.02)
%     updateMethod        'adaptive' (regime-based jump) | 'incremental' (hIn+1) (default: 'adaptive')
%     verbose             Print iteration table to console (default: false)
%     storeResidualPhasor Store full residual PhasorArray in info (default: false)
%
%   Returns:
%     r     PhasorArray solution X
%     info  Solver diagnostics struct — all fields always present:
%       .status          0=converged, 1=stagnated, 2=maxh_reached,
%                        3=fixed_h, 4=algebraic_unreachable
%       .statusMsg       Human-readable exit message
%       .resrelnorm      Final relative residual norm  ||AX-B||_F / (||B||_F + eps)
%       .resnorm         Final absolute residual norm
%       .h               Final hIn used
%       .h_history       Vector of hIn values tried  ([] if autoUpdateh=false)
%       .res_history     Absolute residual history   ([] if autoUpdateh=false)
%       .resrel_history  Relative residual history   ([] if autoUpdateh=false)
%       .regime_history  'initial'/'exponential'/'algebraic'/'stagnated' per iter
%       .s_alg_history   Algebraic slope estimates   ([] if initial regime only)
%       .s_exp_history   Exponential slope estimates ([] if initial regime only)
%       .ressym          ||(r - r')/2||_F  (NaN if r not square)
%       .resasym         ||(r + r')/2||_F  (NaN if r not square)
%       .residualPhasor  Residual PhasorArray ([] unless storeResidualPhasor=true)
%
%   See also: lyap, mrHmcDivide, plotHmcConvergence

arguments
    A PhasorArray
    B PhasorArray
    options.h                    = []
    options.thresholdResidual    (1,1) double  = 5e-4
    options.autoUpdateh          (1,1) logical = false
    options.maxh                 = []
    options.stagnationWindow     (1,1) {mustBeInteger, mustBePositive} = 15
    options.stagnationRatio      (1,1) double  = 0.02
    options.updateMethod         {mustBeMember(options.updateMethod, ...
                                  {'adaptive','incremental'})} = 'adaptive'
    options.verbose              (1,1) logical = false
    options.storeResidualPhasor  (1,1) logical = false
end

thresholdResidual = options.thresholdResidual;
C                 = namedargs2cell(options);

%% --- Guard G1: A is a constant scalar ---

if isscalar(A) && A.h == 0
    scalarVal = A{1,1,0};
    if scalarVal == 0
        error('PhasorArray:mlHmcDivide:singularScalar', 'Division by zero scalar.')
    end
    r    = PhasorArray(B.value / scalarVal);
    info = packInfo(3, 'Scalar constant — direct division.', ...
        0, 0, 0, [], [], [], [], [], [], r, [], options);
    return
end

%% --- Guard G2: A is scalar periodic, B is matrix ---

if isscalar(A) && ~isscalar(B)
    [nB, mB] = size(B, [1 2]);
    r = PhasorArray(nB, mB);
    worst_resrel = 0;
    worst_resnorm = 0;
    for ii = 1:nB
        for jj = 1:mB
            [xij, ~] = mlHmcDivide(A, B{ii,jj}, C{:});
            r{ii,jj} = xij;
        end
    end
    resPhasor    = A*r - B;
    resnorm      = norm(resPhasor.value, 'fro');
    Bnorm        = norm(B.value, 'fro');
    resrelnorm   = resnorm / (Bnorm + eps);
    hFinal       = (size(r.value,3)-1)/2;
    info = packInfo(3, 'Scalar-periodic × matrix — element-wise recursion.', ...
        resrelnorm, resnorm, hFinal, [], [], [], [], [], [], r, resPhasor, options);
    return
end

%% --- Guard G3: A matrix, B scalar ---

if ~isscalar(A) && isscalar(B)
    error('PhasorArray:mlHmcDivide:underdetermined', ...
        ['A\\B with matrix A (%d×%d) and scalar B is ill-defined. ' ...
         'Use A\\(eye(%d)*B) to broadcast B to an identity-scaled matrix.'], ...
        size(A,1), size(A,2), size(A,1))
end

%% --- Initialise hIn ---

if isempty(options.h)
    hIn = max(A.h, B.h) + A.h;
    options.autoUpdateh = true;
else
    hIn = options.h;
end

Bnorm = norm(B.value, 'fro');

%% --- Initial solve ---

[r, resnorm, resrelnorm, resPhasor] = solveAtH(A, B, hIn, Bnorm);

%% --- Fixed-h: early return ---

if ~options.autoUpdateh
    info = packInfo(3, sprintf('Fixed hIn=%d.', hIn), ...
        resrelnorm, resnorm, hIn, [], [], [], [], [], [], r, resPhasor, options);
    return
end

%% --- Adaptive-h loop ---

maxhIn = options.maxh;
if isempty(maxhIn), maxhIn = max(hIn * 20, hIn + 20); end
capacity = maxhIn - hIn + 1;

h_history      = zeros(1, capacity);
res_history    = zeros(1, capacity);
resrel_history = zeros(1, capacity);
regime_history = cell(1, capacity);
s_alg_history  = [];
s_exp_history  = [];

h_history(1)      = hIn;
res_history(1)    = resnorm;
resrel_history(1) = resrelnorm;
regime_history{1} = 'initial';
nIter = 1;

r_best          = r;
resnorm_best    = resnorm;
resrelnorm_best = resrelnorm;
resPhasor_best  = resPhasor;
hIn_best        = hIn;

stagnationWindow     = options.stagnationWindow;
stagnationRatio      = options.stagnationRatio;
spectral_width       = A.h + B.h;
algebraic_hit_count  = 0;

if options.verbose
    fprintf('\nPhasorArray.mlHmcDivide — adaptive hIn\n')
    fprintf('%5s | %11s | %12s | %-12s | %s\n', 'hIn', 'Res norm', 'Rel res norm', 'Regime', 'Note')
    fprintf('------|-------------|--------------|-------------|------\n')
    fprintf('%5d | %11.4e | %12.4e | %-12s |\n', hIn, resnorm, resrelnorm, 'initial')
end

status    = -1;
statusMsg = '';

while resrelnorm_best > thresholdResidual && hIn < maxhIn
    %% --- Step size selection ---
    regime = 'initial';

    if strcmp(options.updateMethod, 'incremental') || hIn < spectral_width * 1.1 || nIter <= 1
        hIn = hIn + 1;
    elseif strcmp(options.updateMethod, 'adaptive')
        idx_start = find(h_history(1:nIter) >= spectral_width * 1.1, 1);
        if isempty(idx_start) || idx_start >= nIter
            % Not enough asymptotic history yet — fall back to +1
            hIn = hIn + 1;
        else
            idx1 = idx_start;
            idx2 = nIter;
            e1   = resrel_history(idx1);  h1 = h_history(idx1);
            e2   = resrel_history(idx2);  h2 = h_history(idx2);

            s_exp = (log(e2+eps) - log(e1+eps)) / (h2 - h1 + eps);
            s_alg = (log(e2+eps) - log(e1+eps)) / (log(h2+eps) - log(h1+eps));

            target = thresholdResidual;
            h_exp  = h2 + ceil((log(target+eps) - log(e2+eps)) / (s_exp - eps));
            h_alg  = ceil(h2 * (target / (e2+eps))^(1 / (s_alg - eps)));

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

            s_alg_history(end+1) = s_alg; %#ok<AGROW>
            s_exp_history(end+1) = s_exp; %#ok<AGROW>

            deltah = ceil(deltah * 0.8);
            deltah = max(1, deltah);
            deltah = min(deltah, 50);
            deltah = min(deltah, ceil(hIn * 0.5));
            hIn    = h2 + deltah;

            % Algebraic early exit
            if strcmp(regime, 'algebraic') && h_alg > maxhIn
                algebraic_hit_count = algebraic_hit_count + 1;
                if algebraic_hit_count >= 2
                    status    = 4;
                    statusMsg = sprintf( ...
                        'Algebraic convergence too slow (slope=%.2f). Target h=%d unreachable (maxhIn=%d). Best: hIn=%d, resrel=%.2e.', ...
                        s_alg, h_alg, maxhIn, hIn_best, resrelnorm_best);
                    r          = r_best;
                    resnorm    = resnorm_best;
                    resrelnorm = resrelnorm_best;
                    resPhasor  = resPhasor_best;
                    hIn        = hIn_best;
                    if options.verbose
                        fprintf('%5d | %11.4e | %12.4e | %-12s | unreachable (slope=%.2f, target h=%d > maxhIn=%d)\n', ...
                            hIn, resnorm, resrelnorm, regime, s_alg, h_alg, maxhIn)
                    end
                    break
                end
            else
                algebraic_hit_count = 0;
            end
        end
    end

    if status == 4, break, end
    hIn = min(hIn, maxhIn);

    %% --- Solve at new hIn ---
    nIter = nIter + 1;
    [r, resnorm, resrelnorm, resPhasor] = solveAtH(A, B, hIn, Bnorm);

    h_history(nIter)      = hIn;
    res_history(nIter)    = resnorm;
    resrel_history(nIter) = resrelnorm;
    regime_history{nIter} = regime;

    if resrelnorm < resrelnorm_best
        r_best          = r;
        resnorm_best    = resnorm;
        resrelnorm_best = resrelnorm;
        resPhasor_best  = resPhasor;
        hIn_best        = hIn;
    end

    note = '';

    % Convergence check (inside loop for correct table alignment)
    if resrelnorm <= thresholdResidual
        status    = 0;
        statusMsg = sprintf('Converged at hIn=%d (resrel=%.2e).', hIn, resrelnorm);
        note = 'converged';
        if options.verbose
            fprintf('%5d | %11.4e | %12.4e | %-12s | %s\n', hIn, resnorm, resrelnorm, regime, note)
        end
        break
    end

    % Stagnation check
    if nIter >= stagnationWindow
        window     = resrel_history(nIter - stagnationWindow + 1 : nIter);
        rel_improv = (window(1) - min(window)) / (window(1) + eps);
        if rel_improv < stagnationRatio
            status    = 1;
            statusMsg = sprintf('Stagnated at hIn=%d (%.1f%% improvement over %d steps). Best: hIn=%d, resrel=%.2e.', ...
                hIn, rel_improv*100, stagnationWindow, hIn_best, resrelnorm_best);
            note = 'stagnated';
        end
    end

    if options.verbose
        fprintf('%5d | %11.4e | %12.4e | %-12s | %s\n', hIn, resnorm, resrelnorm, regime, note)
    end

    if status == 1
        r          = r_best;
        resnorm    = resnorm_best;
        resrelnorm = resrelnorm_best;
        resPhasor  = resPhasor_best;
        hIn        = hIn_best;
        break
    end
end

%% --- Finalise status ---

if status == -1
    if resrelnorm_best <= thresholdResidual
        status    = 0;
        statusMsg = sprintf('Converged at hIn=%d (resrel=%.2e).', hIn_best, resrelnorm_best);
        r          = r_best;
        resnorm    = resnorm_best;
        resrelnorm = resrelnorm_best;
        resPhasor  = resPhasor_best;
        hIn        = hIn_best;
    else
        status    = 2;
        statusMsg = sprintf('Reached maxhIn=%d without convergence. Best: hIn=%d, resrel=%.2e.', ...
            maxhIn, hIn_best, resrelnorm_best);
        r          = r_best;
        resnorm    = resnorm_best;
        resrelnorm = resrelnorm_best;
        resPhasor  = resPhasor_best;
        hIn        = hIn_best;
        if options.verbose
            fprintf('  → maxhIn reached. Returning best solution (hIn=%d).\n', hIn_best)
        end
    end
end

h_history      = h_history(1:nIter);
res_history    = res_history(1:nIter);
resrel_history = resrel_history(1:nIter);
regime_history = regime_history(1:nIter);

info = packInfo(status, statusMsg, resrelnorm, resnorm, hIn, ...
    h_history, res_history, resrel_history, regime_history, ...
    s_alg_history, s_exp_history, r, resPhasor, options);

end % mlHmcDivide

%% =========================================================================
function [r, resnorm, resrelnorm, resPhasor] = solveAtH(A, B, hIn, Bnorm)
%SOLVEATH  Solve A(t)*X(t) = B(t) at a given hIn via Toeplitz least-squares.
hOut    = hIn + A.h;
res_tb  = A.spTB(hOut, hIn) \ B.spF_tb(hOut);
r       = PhasorArray(TFTB_2_array(res_tb, size(A,2), size(B,2)));
resPhasor  = A*r - B;
resnorm    = norm(resPhasor.value, 'fro');
resrelnorm = resnorm / (Bnorm + eps);
end

%% =========================================================================
function info = packInfo(status, statusMsg, resrelnorm, resnorm, h, ...
        h_history, res_history, resrel_history, regime_history, ...
        s_alg_history, s_exp_history, r, resPhasor, options)
%PACKINFO  Build the info struct with all fields always present.
info.status         = status;
info.statusMsg      = statusMsg;
info.resrelnorm     = resrelnorm;
info.resnorm        = resnorm;
info.h              = h;
info.h_history      = h_history;
info.res_history    = res_history;
info.resrel_history = resrel_history;
info.regime_history = regime_history;
info.s_alg_history  = s_alg_history;
info.s_exp_history  = s_exp_history;

if size(r,1) == size(r,2)
    info.ressym  = norm(value(r - r')/2, 'fro');
    info.resasym = norm(value(r + r')/2, 'fro');
else
    info.ressym  = NaN;
    info.resasym = NaN;
end

if options.storeResidualPhasor
    info.residualPhasor = resPhasor;
else
    info.residualPhasor = [];
end
end
