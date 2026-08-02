function [res, info] = lyap(pA1, pA2, pA3, nvp)
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
    pA1                                                        % A (both modes)
    pA2                                                        % Q (Lyapunov) or B (Sylvester)
    pA3                               = []                     % [] → Lyapunov; provided → Sylvester
    nvp.T                   (1,1) double  = 2*pi
    nvp.h                                 = []            % [] triggers autoUpdateh
    nvp.thresholdResidual   (1,1) double  = 1e-6
    nvp.autoUpdateh         (1,1) logical = false
    nvp.maxh                              = []            % default: h0 * 20
    nvp.stagnationWindow    (1,1) {mustBeInteger, mustBePositive} = 5
    nvp.stagnationRatio     (1,1) double  = 0.05
    nvp.verbose             (1,1) logical = false
    nvp.storeResidualPhasor (1,1) logical = false
    nvp.systemType          {mustBeMember(nvp.systemType,  {'rectangle','square'})}   = 'rectangle'
    nvp.updateMethod        {mustBeMember(nvp.updateMethod,{'adaptive','incremental'})} = 'adaptive'
    nvp.direction           {mustBeMember(nvp.direction,   {'backward','forward'})}      = 'backward'
    nvp.E                                 = []            % descriptor mass (Lyapunov); [] → standard
    nvp.Ea                                = []            % descriptor left mass (Sylvester)
    nvp.Eb                                = []            % descriptor right mass (Sylvester)
    nvp.derivativeForm      {mustBeMember(nvp.derivativeForm,{'product','sandwich'})}  = 'product'
end

% --- Descriptor (generalized) mode: delegate to the lyapG engine -----------
% E / Ea / Eb supplied via name-value → solve the descriptor equation
%   Lyapunov : d/dt(E'PE) + A'PE + E'PA + Q = 0
%   Sylvester: d/dt(Ea X Eb) + A X Eb + Ea X B + C = 0
% The standard (no E) path falls through to the base solver below.
if ~isempty(nvp.E) || ~isempty(nvp.Ea) || ~isempty(nvp.Eb)
    fwd = {'T',nvp.T, 'h',nvp.h, 'thresholdResidual',nvp.thresholdResidual, ...
           'autoUpdateh',nvp.autoUpdateh, 'maxh',nvp.maxh, ...
           'stagnationWindow',nvp.stagnationWindow, 'stagnationRatio',nvp.stagnationRatio, ...
           'verbose',nvp.verbose, 'storeResidualPhasor',nvp.storeResidualPhasor, ...
           'systemType',nvp.systemType, 'updateMethod',nvp.updateMethod, ...
           'derivativeForm',nvp.derivativeForm, 'direction',nvp.direction};
    if isempty(pA3)
        [res, info] = lyapG(pA1, pA2, nvp.E, [], [], fwd{:});   % Lyapunov descriptor
    else
        [res, info] = lyapG(pA1, pA2, pA3, nvp.Ea, nvp.Eb, fwd{:}); % Sylvester descriptor
    end
    return
end

T                 = nvp.T;
omega             = 2*pi / T;
thresholdResidual = nvp.thresholdResidual;
autoUpdateh       = nvp.autoUpdateh;

%% --- Mode dispatch, validation, and h initialisation ---

isLyapunov = isempty(pA3);

if isLyapunov
    % Lyapunov mode: dP/dt + A'P + PA + Q = 0
    if ~isa(pA1,'PhasorArray'), pA1 = PhasorArray(pA1); end
    if ~isa(pA2,'PhasorArray'), pA2 = PhasorArray(pA2); end
    if size(pA1,1) ~= size(pA1,2)
        error('PhasorArray:lyap:nonSquareA', 'A must be square (got %dx%d).', size(pA1,1), size(pA1,2))
    end
    if size(pA1,1)~=size(pA2,1) || size(pA1,2)~=size(pA2,2)
        error('PhasorArray:lyap:dimensionMismatch', 'A (%dx%d) and Q (%dx%d) must be the same size.', size(pA1,1), size(pA1,2), size(pA2,1), size(pA2,2))
    end
    h = nvp.h;
    if isempty(h)
        h = max(pA1.h, pA2.h);
        autoUpdateh = true;
    end
    % Remap to Sylvester: dP/dt + (A').P + P.(A) + Q = 0  →  A'→pA1, A→pA2, Q→pA3
    pA3 = pA2; pA2 = pA1; pA1 = pA1.';

else
    % Sylvester mode: dM/dt + AM + MB + C = 0
    if ~isa(pA1,'PhasorArray'), pA1 = PhasorArray(pA1); end
    if ~isa(pA2,'PhasorArray'), pA2 = PhasorArray(pA2); end
    if ~isa(pA3,'PhasorArray'), pA3 = PhasorArray(pA3); end
    if size(pA1,1) ~= size(pA1,2)
        error('PhasorArray:lyap:nonSquareA', 'Sylvester: A must be square (got %dx%d).', size(pA1,1), size(pA1,2))
    end
    if size(pA2,1) ~= size(pA2,2)
        error('PhasorArray:lyap:nonSquareB', 'Sylvester: B must be square (got %dx%d).', size(pA2,1), size(pA2,2))
    end
    if size(pA1,1) ~= size(pA3,1)
        error('PhasorArray:lyap:dimensionMismatch', 'Sylvester: rows(A)=%d must equal rows(C)=%d.', size(pA1,1), size(pA3,1))
    end
    if size(pA2,2) ~= size(pA3,2)
        error('PhasorArray:lyap:dimensionMismatch', 'Sylvester: cols(B)=%d must equal cols(C)=%d.', size(pA2,2), size(pA3,2))
    end
    h = nvp.h;
    if isempty(h)
        h = max([pA1.h, pA2.h, pA3.h]);
        autoUpdateh = true;
    end
end

%% --- Single-order solve callback ---

% Closes over the operands so adaptiveHSolve only has to vary h.
solveAtH = @(hh) solveOnce(pA1, pA2, pA3, hh, omega, T, nvp);

%% --- Fixed-h: early return ---

if ~autoUpdateh
    [res, resnorm, resrelnorm, resPhasor] = solveAtH(h);
    info = packInfo(3, sprintf('Fixed h=%d.', h), ...
        resrelnorm, resnorm, h, [], [], [], [], {}, [], [], res, resPhasor, nvp);
    return
end

%% --- Adaptive-h refinement ---

n_sys    = size(pA1, 1) * size(pA2, 1);   % n1*n2: size of the vectorised solution
h_op     = max([pA1.h, pA2.h, pA3.h]);     % spectral width of the operators (fixed)
isSquare = strcmp(nvp.systemType, 'square');

cfg = struct( ...
    'thresholdResidual', thresholdResidual, ...
    'maxh',              nvp.maxh, ...
    'stagnationWindow',  nvp.stagnationWindow, ...
    'stagnationRatio',   nvp.stagnationRatio, ...
    'updateMethod',      nvp.updateMethod, ...
    'verbose',           nvp.verbose, ...
    'hOp',               h_op, ...
    'hOutFcn',           @(hh) hh*isSquare + (hh + h_op)*~isSquare, ...
    'preamble',          buildPreamble(pA3, n_sys, h_op, isLyapunov, isSquare), ...
    'label',             'h');

[best, trace] = adaptiveHSolve(solveAtH, h, cfg);

res  = best.sol;
info = packInfo(trace.status, trace.statusMsg, best.resrelnorm, best.resnorm, best.h, ...
    trace.h_history, trace.resrel_history, trace.res_history, trace.time_history, ...
    trace.regime_history, trace.s_alg_history, trace.s_exp_history, ...
    best.sol, best.resPhasor, nvp);

end % lyap

%% =========================================================================
function [res, resnorm, resrelnorm, resPhasor] = solveOnce(pA1, pA2, pA3, h, omega, T, nvp)
%SOLVEONCE  Solve the Sylvester/Lyapunov equation at a single harmonic order.
res = PhasorArray(SylvHarmonic(pA1, pA2, pA3, h, omega, nvp.systemType, nvp.direction));
[resnorm, resrelnorm, resPhasor] = computeResidual(res, pA1, pA2, pA3, T, nvp.direction);
end

%% =========================================================================
function s = buildPreamble(pA3, n_sys, h_op, isLyapunov, isSquare)
%BUILDPREAMBLE  Header text printed above the verbose refinement table.
n1 = size(pA3, 1);  n2 = size(pA3, 2);
s = sprintf('\nPhasorArray.lyap — adaptive h refinement\n');
if isLyapunov
    s = [s sprintf('  Equation:  dP/dt + A''·P + P·A + Q = 0   [Lyapunov, %dx%d]\n', n1, n2)];
    s = [s sprintf('  Each step solves:  E · vec(X) = vec(Q)\n')];
    hop_label = 'max(hA,hQ)';
else
    s = [s sprintf('  Equation:  dM/dt + A·M + M·B + C = 0   [Sylvester, %dx%d]\n', n1, n2)];
    s = [s sprintf('  Each step solves:  E · vec(X) = vec(C)\n')];
    hop_label = 'max(hA,hB,hC)';
end
if isSquare
    s = [s sprintf('  where  n = n1·n2 = %d·%d = %d,  h_op = %s = %d,  hOut = h  [square mode]\n', n1, n2, n_sys, hop_label, h_op)];
    s = [s sprintf('    E  :  [n·(2·h+1)] × [n·(2·h+1)]   (square — border harmonics discarded)\n')];
else
    s = [s sprintf('  where  n = n1·n2 = %d·%d = %d,  h_op = %s = %d,  hOut = h + h_op\n', n1, n2, n_sys, hop_label, h_op)];
    s = [s sprintf('    E  :  [n·(2·hOut+1)] × [n·(2·h+1)]   (rectangle — full residual L2 coverage)\n')];
end
s = [s sprintf('    X  :  solution truncated to h harmonics\n\n')];
end

%% =========================================================================
function [resnorm, resrelnorm, resPhasor] = computeResidual(res, pA1, pA2, pA3, T, direction)
%COMPUTERESIDUAL  Evaluate Sylvester residual and its norms.
% 'backward': dX/dt + AX + XB + C = 0 ; 'forward': -dX/dt + AX + XB + C = 0.
if strcmp(direction, 'forward')
    signD = -1;
else
    signD = +1;
end
resPhasor  = signD*res.d(T) + pA1*res + res*pA2 + pA3;
resnorm    = norm(resPhasor.value, 'fro');
Cnorm      = norm(pA3.value, 'fro');
resrelnorm = resnorm / (Cnorm + eps);
end

%% =========================================================================
function info = packInfo(status, statusMsg, resrelnorm, resnorm, h, ...
        h_history, resrel_history, res_history, time_history, ...
        regime_history, s_alg_history, s_exp_history, res, resPhasor, nvp)
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

if nvp.storeResidualPhasor
    info.residualPhasor = resPhasor;
else
    info.residualPhasor = [];
end
end
