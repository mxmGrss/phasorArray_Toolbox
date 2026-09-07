function info = packSolverInfo(best, trace, nvp)
%PACKSOLVERINFO  Build the diagnostics struct every harmonic solver returns.
%
%   SYNTAX
%     info = packSolverInfo(best, trace)
%     info = packSolverInfo(best, trace, solution=X, storeResidualPhasor=tf, ...
%                           residualPhasor=R, extra=s)
%
%   DESCRIPTION
%   One definition of the `info` contract, for every solver that reports one:
%   PhasorArray/lyap, /lyapG, /mlHmcDivide and /place. Before this existed the
%   struct was assembled at five sites -- three copies of a local `packInfo`
%   plus two inline `struct(...)` in place -- which had already drifted: the
%   two `packInfo` signatures passed res_history and resrel_history in opposite
%   positions, and the symmetry diagnostic came under three names with a factor
%   two between them.
%
%   BEST and TRACE are what ADAPTIVEHSOLVE returns, passed straight through. A
%   solver that did not run the adaptive loop (a fixed-h branch) builds them
%   itself with the values it has; any field left out lands on the documented
%   default below rather than being absent, so the field set never varies.
%
%   GUARANTEED FIELDS
%     status         0 converged, 1 stagnated, 2 maxh reached, 3 fixed h,
%                    4 convergence judged unreachable
%     statusMsg      the same, in words
%     h              order the accepted solution was computed at
%     resnorm        absolute residual norm
%     resrelnorm     residual norm relative to the right-hand side
%     solskewnorm    norm of the skew-Hermitian part of the SOLUTION, i.e. how
%                    far a solution that should be Hermitian actually is -- the
%                    Lyapunov and Riccati solutions are Hermitian in exact
%                    arithmetic, so this is a drift measure, not a residual.
%                    NaN when the solution is not square, or still carries an
%                    unsolved decision variable. A NORM, homogeneous with
%                    resnorm. The Hermitian part needs no field of its own: the
%                    two parts are orthogonal, so its norm is
%                    sqrt(energy(X) - solskewnorm^2) for a caller holding X.
%                    (This is what lyap called resPsym -- "residual of P's
%                    symmetry", which is not a residual -- and what mlHmcDivide
%                    called ressym, at half the magnitude.)
%     residualPhasor the residual itself, or [] unless storeResidualPhasor
%     h_history, res_history, resrel_history, time_history, regime_history,
%     s_alg_history, s_exp_history
%                    one entry per refinement step; empty ([] or {}) when the
%                    solver ran at a fixed order
%     hForTargetResidual, targetResidual
%                    the order a near-zero residual would need, extrapolated at
%                    exit, and the residual that extrapolation aimed at. NaN
%                    when no refinement ran.
%
%   INPUTS
%     best   - struct from adaptiveHSolve: h, resnorm, resrelnorm.
%     trace  - struct from adaptiveHSolve: status, statusMsg, the histories,
%              hForTargetResidual, targetResidual.
%     solution            - The returned solution, as a PhasorArray, for
%                           solskewnorm. Omit to get NaN.
%     storeResidualPhasor - Publish the residual in info.residualPhasor.
%     residualPhasor      - The value to publish when the flag is set.
%     extra               - Struct of solver-specific fields, merged last. Use
%                           it only for what a single solver reports; anything
%                           every solver reports belongs in the list above.
%
%   See also adaptiveHSolve, PhasorArray/lyap, PhasorArray/mlHmcDivide, mherm

arguments
    best  (1,1) struct
    trace (1,1) struct
    nvp.solution = []
    nvp.storeResidualPhasor (1,1) logical = false
    nvp.residualPhasor = []
    nvp.extra (1,1) struct = struct()
end

info.status         = getOr(trace, 'status', NaN);
info.statusMsg      = getOr(trace, 'statusMsg', '');
info.h              = getOr(best,  'h', NaN);
info.resnorm        = getOr(best,  'resnorm', NaN);
info.resrelnorm     = getOr(best,  'resrelnorm', NaN);

% Symmetry drift of the solution. NaN rather than an error on the two cases
% where it has no meaning: a rectangular solution has no Hermitian part, and one
% still carrying an unsolved decision variable has no numeric norm.
info.solskewnorm = NaN;
if ~isempty(nvp.solution)
    X = nvp.solution;
    if ~isa(X, 'PhasorArray'), X = PhasorArray(X); end
    if size(X,1) == size(X,2) && isnumeric(value(X))
        info.solskewnorm = sqrt(hermEnergy(X, false, skewOption='skew'));
    end
end

info.h_history      = getOr(trace, 'h_history', []);
info.res_history    = getOr(trace, 'res_history', []);
info.resrel_history = getOr(trace, 'resrel_history', []);
info.time_history   = getOr(trace, 'time_history', []);
info.regime_history = getOr(trace, 'regime_history', {});
info.s_alg_history  = getOr(trace, 's_alg_history', []);
info.s_exp_history  = getOr(trace, 's_exp_history', []);

info.hForTargetResidual = getOr(trace, 'hForTargetResidual', NaN);
info.targetResidual     = getOr(trace, 'targetResidual', NaN);

if nvp.storeResidualPhasor
    info.residualPhasor = nvp.residualPhasor;
else
    info.residualPhasor = [];
end

% Solver-specific fields last, so a solver can add to the contract but never
% silently redefine one of its guaranteed fields.
for f = string(fieldnames(nvp.extra)).'
    if isfield(info, f)
        error('packSolverInfo:reservedField', ...
            'extra.%s would overwrite a guaranteed info field.', f);
    end
    info.(f) = nvp.extra.(f);
end
end

% =========================================================================
function v = getOr(s, name, default)
%GETOR  Field of s, or default when the caller did not supply it.
if isfield(s, name)
    v = s.(name);
else
    v = default;
end
end
