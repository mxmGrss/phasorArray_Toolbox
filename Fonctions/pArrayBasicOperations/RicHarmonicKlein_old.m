function [K_final, S_final, htrunc, H] = RicHarmonicKlein(A, B, Q, R, K0, T, options)
%RIC_HARMONIC_KLEIN Iterative Riccati algorithm (Kleinman)
%   A: System matrix (PhasorArray)
%   B: Input matrix (PhasorArray)
%   Q: State weighting matrix
%   R: Control weighting matrix
%   K0: Initial feedback gain
%   T: period of periodic matrices
%   options.max_iter:           Maximum Klein iterations (default: 100)
%   options.htrunc:             Initial harmonic truncation order (auto if empty)
%   options.autoUpdateh:        Auto-update h in Lyapunov solve (default: false)
%   options.residualThreshold:  Convergence threshold (default: 1e-6)
%   options.relChangeThreshold: Relative change threshold on S (default: 1e-3)
%   options.hmax:               Maximum harmonic order (default: inf)
%   options.stagnationWindow:   Sliding window for stagnation detection (default: 5)
%   options.stagnationRatio:    Min relative improvement to avoid stagnation (default: 0.05)

arguments (Input)
    A PhasorArray
    B PhasorArray
    Q
    R
    K0 =ones(size(A,1),size(B,2))'*10
    T = 2*pi
    options.max_iter (1,1) {mustBeInteger, mustBePositive} = 100
    options.htrunc = []
    options.autoUpdateh logical = false
    options.residualThreshold = 1e-3
    options.relChangeThreshold = 1e-3
    options.hmax = min(2* sum([size(A,3), size(B,3), size(Q,3), size(R,3), size(K0,3)]),4* max([size(A,3), size(B,3), size(Q,3), size(R,3), size(K0,3)]))
    options.stagnationWindow (1,1) {mustBeInteger, mustBePositive} = 5
    options.stagnationRatio  (1,1) double = 0.05
end
max_iter          = options.max_iter;
htrunc            = options.htrunc;
autoUpdateh       = options.autoUpdateh;
residualThreshold = options.residualThreshold;

if isempty(htrunc)
    htrunc      = max([size(A,3), size(B,3), size(Q,3), size(R,3), size(K0,3)]);
    autoUpdateh = true;
end

Ak0 = A - B * K0;
LL  = HmqNEig(Ak0, htrunc, T, "fundamental");
if any(real(LL) > 0)
    % warning('RicHarmonicKlein:unstableK0', ...
        % 'K0 does not stabilize A-B*K0. Attempting LQR fallback on DC component...');
    Q_tmp = PhasorArray(Q);
    R_tmp = PhasorArray(R);
    A_dc = A.phas(0);
    B_dc = B.phas(0);
    Q_dc = Q_tmp.phas(0);
    R_dc = R_tmp.phas(0);
    [~, K0_lqr] = icare(A_dc, B_dc, Q_dc, R_dc);
    if isempty(K0_lqr) || any(~isfinite(K0_lqr(:)))
        error('RicHarmonicKlein:noStabilizingK0', ...
            'LQR on DC component also failed. Provide a stabilizing K0 manually.');
    end
    K0  = PhasorArray(K0_lqr);
    Ak0 = A - B * K0;
    LL  = HmqNEig(Ak0, htrunc, T, "fundamental");
    if any(real(LL) > 0)
        error('RicHarmonicKlein:noStabilizingK0', ...
            'LQR DC fallback does not stabilize the periodic system. Provide a stabilizing K0 manually.');
    end
    % warning('RicHarmonicKlein:K0Replaced', 'Using LQR DC fallback as K0.');
end

Q = PhasorArray(Q);
R = PhasorArray(R);

outIsReal = all([isreal(A), isreal(B), isreal(Q), isreal(R)]);

Rinv  = inv(R);
Kk{1} = K0;
S     = cell(1, max_iter);

% Outer-loop stagnation tracking
resRicnorm     = NaN;
resRic_history = [];
S_best         = [];
resRic_best    = Inf;
kk_best        = 1;

for kk = 1:max_iter
    if kk > 1
        Kk{kk} = Rinv * B.' * S{kk-1};
        if outIsReal
            Kk{kk} = mreal(Kk{kk});
        end
    end

    Ak{kk} = A - B * Kk{kk};
    Yk{kk} = Kk{kk}.' * R * Kk{kk};
    if outIsReal
        Ak{kk} = mreal(Ak{kk});
        Yk{kk} = mreal(Yk{kk});
    end

    % Solve Lyapunov: Ak'*S + S*Ak + d(S) + (Yk+Q) = 0
    % lyap handles h-update and stagnation detection internally when autoUpdateh=true
    [S{kk}, res] = lyap( ...
        reduce(Ak{kk}, 'reduceMethod','relative','reduceThreshold',residualThreshold/100,'exclude0Phasor',false), ...
        reduce(Yk{kk}+Q, 'reduceMethod','relative','reduceThreshold',residualThreshold/100,'exclude0Phasor',false), ...
        "T", T, "h", htrunc, ...
        "autoUpdateh",       autoUpdateh, ...
        "maxh",              options.hmax, ...
        "thresholdResidual", residualThreshold, ...
        "stagnationWindow",  options.stagnationWindow, ...
        "stagnationRatio",   options.stagnationRatio); %#ok<NASGU>

    % Recover h actually used by lyap (via h_history, populated when autoUpdateh=true)
    if autoUpdateh && isfield(res, 'h_history') && ~isempty(res.h_history)
        htrunc = res.h_history(end);
    end
    H(kk) = htrunc;

    if outIsReal
        S{kk} = mreal(S{kk});
    end

    % Check convergence (outer Klein loop)
    if kk > 1
        riccati_residual = d(S{kk}, T) + A.' * S{kk} + S{kk} * A ...
                           - S{kk} * B * Rinv * B.' * S{kk} + Q;
        Qnorm      = norm(value(Q), 'fro');
        resRicnorm = norm(value(riccati_residual), 'fro') / (Qnorm + eps);
        rel_change = norm(value(S{kk} - S{kk-1}), 'fro') / (norm(value(S{kk}), 'fro') + eps);

        % Best-solution tracking
        if resRicnorm < resRic_best
            resRic_best = resRicnorm;
            S_best      = S{kk};
            kk_best     = kk;
        end
        resRic_history(end+1) = resRicnorm; %#ok<AGROW>

        % Stagnation detection over sliding window
        if numel(resRic_history) >= options.stagnationWindow
            window     = resRic_history(end-options.stagnationWindow+1 : end);
            rel_improv = (window(1) - min(window)) / (window(1) + eps);
            if rel_improv < options.stagnationRatio
                warning('RicHarmonicKlein:stagnation', ...
                    ['Klein: Riccati relative residual stagnated at %.2e after iter %d ' ...
                     '(%.1f%% improvement over %d steps).\n' ...
                     '  Returning best solution found (iter %d, resrel=%.2e).'], ...
                    resRicnorm, kk, rel_improv*100, options.stagnationWindow, ...
                    kk_best, resRic_best);
                S{kk} = S_best;
                break
            end
        end

        % Normal convergence
        if rel_change < options.relChangeThreshold || resRicnorm < residualThreshold
            fprintf('Converged at iteration %d\n', kk);
            break;
        end
    end

end

S_final = S{kk};
K_final = Kk{kk};

fprintf('Riccati relative residual norm: %.2e\n', resRicnorm);
end
