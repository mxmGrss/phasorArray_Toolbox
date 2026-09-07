function [r, info] = mlHmcDivide(A, B, nvp)
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
%     thresholdResidual   Convergence threshold on the relative residual
%                         (default: 5e-4). Relative to the right-hand side:
%                         norm(A*X-B,'fro') / norm(B,'fro'), never to X.
%                         Tightening it raises the order of X: cheap when the
%                         quotient is smooth, not when it is not. See
%                         ADAPTIVEHSOLVE, and info.resrel_history for the
%                         curve on your own problem.
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
%       .time_history    Step solve time in seconds  ([] if autoUpdateh=false)
%       .regime_history  'initial'/'exponential'/'algebraic'/'stagnated' per iter
%       .s_alg_history   Algebraic slope estimates   ([] if initial regime only)
%       .s_exp_history   Exponential slope estimates ([] if initial regime only)
%       .solskewnorm     ||(X - X')/2||_F, how far the solution is from
%                        Hermitian (NaN if it is not square). Replaces .ressym;
%                        the Hermitian part, formerly .resasym, is recoverable
%                        as sqrt(energy(X) - solskewnorm^2) since the two parts
%                        are orthogonal.
%       .residualPhasor  Residual PhasorArray ([] unless storeResidualPhasor=true)
%
%   See also: lyap, mrHmcDivide, plotHmcConvergence

arguments
    A PhasorArray
    B PhasorArray
    nvp.h                    = []
    nvp.thresholdResidual    (1,1) double  = 5e-4
    nvp.autoUpdateh          (1,1) logical = false
    nvp.maxh                 = []
    nvp.stagnationWindow     (1,1) {mustBeInteger, mustBePositive} = 15
    nvp.stagnationRatio      (1,1) double  = 0.02
    nvp.updateMethod         {mustBeMember(nvp.updateMethod, ...
                                  {'adaptive','incremental'})} = 'adaptive'
    nvp.verbose              (1,1) logical = false
    nvp.storeResidualPhasor  (1,1) logical = false
end

thresholdResidual = nvp.thresholdResidual;
C                 = namedargs2cell(nvp);

%% --- Guard G1: A is a constant scalar ---

if isscalar(A) && A.h == 0
    scalarVal = A{1,1,0};
    if scalarVal == 0
        error('PhasorArray:mlHmcDivide:singularScalar', 'Division by zero scalar.')
    end
    r    = PhasorArray(B.value / scalarVal);
    info = packSolverInfo( ...
        struct('h', 0, 'resnorm', 0, 'resrelnorm', 0), ...
        struct('status', 3, 'statusMsg', 'Scalar constant, direct division.'), ...
        solution            = r, ...
        storeResidualPhasor = nvp.storeResidualPhasor);
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
    info = packSolverInfo( ...
        struct('h', hFinal, 'resnorm', resnorm, 'resrelnorm', resrelnorm), ...
        struct('status', 3, 'statusMsg', 'Scalar-periodic × matrix, element-wise recursion.'), ...
        solution            = r, ...
        storeResidualPhasor = nvp.storeResidualPhasor, ...
        residualPhasor      = resPhasor);
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

if isempty(nvp.h)
    hIn = max(A.h, B.h) + A.h;
    nvp.autoUpdateh = true;
else
    hIn = nvp.h;
end

Bnorm = norm(B.value, 'fro');

%% --- Single-order solve callback ---

% Closes over the operands so adaptiveHSolve only has to vary hIn.
solve = @(hh) solveAtH(A, B, hh, Bnorm);

%% --- Fixed-h: early return ---

if ~nvp.autoUpdateh
    [r, resnorm, resrelnorm, resPhasor] = solve(hIn);
    % No refinement ran, so the histories and the extrapolated order stay at
    % their documented defaults -- packSolverInfo fills them in.
    info = packSolverInfo( ...
        struct('h', hIn, 'resnorm', resnorm, 'resrelnorm', resrelnorm), ...
        struct('status', 3, 'statusMsg', sprintf('Fixed hIn=%d.', hIn)), ...
        solution            = r, ...
        storeResidualPhasor = nvp.storeResidualPhasor, ...
        residualPhasor      = resPhasor);
    return
end

%% --- Adaptive-hIn refinement ---

% Spectral width of the operator: the Cauchy product A*X spreads B's harmonics
% by A.h, so the residual only enters its asymptotic regime past A.h + B.h.
spectral_width = A.h + B.h;

cfg = struct( ...
    'thresholdResidual', thresholdResidual, ...
    'maxh',              nvp.maxh, ...
    'stagnationWindow',  nvp.stagnationWindow, ...
    'stagnationRatio',   nvp.stagnationRatio, ...
    'updateMethod',      nvp.updateMethod, ...
    'verbose',           nvp.verbose, ...
    'hOp',               spectral_width, ...
    'hOutFcn',           @(hh) hh + A.h, ...
    'preamble',          sprintf('\nPhasorArray.mlHmcDivide — adaptive hIn\n  Equation:  A(t)·X(t) = B(t)   [%dx%d \\ %dx%d]\n  Each step solves:  T_tb(A) · F_tb(X) = F_tb(B)  (least squares)\n\n', ...
                                 size(A,1), size(A,2), size(B,1), size(B,2)), ...
    'label',             'hIn');

[best, trace] = adaptiveHSolve(solve, hIn, cfg);

r    = best.sol;
info = packSolverInfo(best, trace, ...
    solution            = best.sol, ...
    storeResidualPhasor = nvp.storeResidualPhasor, ...
    residualPhasor      = best.resPhasor);

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
