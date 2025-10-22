function [K_final, S_final,htrunc,H] = RicHarmonicKlein(A, B, Q, R, K0,T, options)
%RIC_HARMONIC_KLEIN Iterative Riccati algorithm
%   A: System matrix
%   B: Input matrix
%   Q: State weighting matrix
%   R: Control weighting matrix
%   K0: Initial feedback gain
%   T: period of periodic matrices
%   max_iter: Maximum number of iterations (optional)
%   htrunc: truncature ordrer for Lyapunov Solving

arguments (Input)
    A PhasorArray
    B PhasorArray
    Q
    R
    K0
    T = 2*pi
    options.max_iter = 100; % Default value for max_iter
    options.htrunc = []
    options.autoUpdateh logical = false
    options.residualThreshold = 1e-6
    options.hmax = inf
end

max_iter = options.max_iter;
htrunc = options.htrunc;
autoUpdateh = options.autoUpdateh;
residualThreshold = options.residualThreshold;

if isempty(htrunc)
    htrunc = max([size(A, 3), size(B, 3), size(Q, 3), size(R, 3), size(K0, 3)]);
    autoUpdateh = true;
end


Ak0 = A - B * K0;
LL = HmqNEig(Ak0,htrunc,T,"fundamental");
if any(real(LL)>0)
    error("Unstable System, please pick K0 such that A(t)-B(t)K0(t) is stable")
end


Rinv = inv(R);
Kk{1} = K0;
S = cell(1, max_iter); % Preallocate S for max_iter iterations

for kk = 1:max_iter
    if kk > 1
        Kk{kk} = Rinv * B.' * S{kk-1};
    end
    Ak{kk} = A - B * Kk{kk};
    Yk{kk} = Kk{kk}.' * R * Kk{kk};
    
    % Solve Lyapunov: Ak’*S + S*Ak + (Yk+Q) + dot(S) = 0
    
    warning('off','phasorArray:lyap:residual')
    [S{kk},res] = lyap(Ak{kk},Yk{kk}+Q,"T",T,"h",htrunc);
    while res.resnorm >residualThreshold && autoUpdateh && htrunc<options.hmax
        htrunc = htrunc+1;
        [S{kk},res] = lyap(Ak{kk},Yk{kk}+Q,"T",T,"h",htrunc);
    end
    % M = kron(eye(2), Ak{kk}.' ) + kron(Ak{kk}.' , eye(2));
    % Ss = (-M * T_tb(htrunc) - N_tb(M, htrunc, T)) \ (F_tb(vec(Yk{kk} + Q), htrunc));
    % SS = F_tb_2_PhasorArray(Ss, 2, 2);
    % S{kk} = SS;
    H(kk) = htrunc;

    % Check convergence
    if kk > 1
        rel_change = norm(value(S{kk} - S{kk-1}), 'fro') / norm(value(S{kk}), 'fro');
        if rel_change < 1e-8
            fprintf('Converged at iteration %d\n', kk);
            break;
        end
    end

end

% Store the final feedback gain and Lyapunov matrix S
S_final = S{kk};
K_final = Kk{kk};


% Verify Riccati equation
riccati_residual = d(S{kk}, T) + A.' * S{kk} + S{kk} * A - S{kk} * B * Rinv * B.' * S{kk} + Q;
fprintf('Riccati residual norm: %.2e\n', norm(value(riccati_residual), 'fro'));
end