% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% ── Polytopic Stability Test ────────────────────────────────────────────
% Is there a common periodic Lyapunov function P(θ) > 0 for all ω ∈ [ω₁, ω₂]?

%% 0. Problem data (replace with your system)
nx = 3;  ha = 3;
rng(42);
A0 = PhasorArray(-5*eye(nx)) + PhasorArray.random(nx, nx, ha) * 0.1;
A1 = PhasorArray.random(nx, nx, ha) * 0.1;
A1 = A1 - A1.';  % Skew-symmetric (like Coriolis) so it stays stable at high omega
omega_bornes = [20, 40] * 2*pi;       % [rad/s]

hP   = 10;     % Lyapunov harmonic order
hlmi = 20;     % LMI truncation order

%% 1. Decision variable
P   = PhasorArray.ndsdpvar(nx, nx, hP);
PTB = P.T_tb(hlmi);

%% 2. Build LMI at each vertex
constraints = [PTB >= 1e-7 * eye(size(PTB))];

for i = 1:numel(omega_bornes)
    omega_i = omega_bornes(i);
    T_i     = 2*pi / omega_i;
    Ai      = A0 + omega_i * A1;

    % Lyapunov inequality:  Ai'P + PAi - ω(NP - PN) < 0
    AiP     = Ai.' * P;
    AiP_TB  = AiP.T_tb(hlmi) - N_tb(nx, hlmi, T_i)' * PTB;

    constraints = [constraints; ...
        AiP_TB + AiP_TB' <= -1e-4 * eye(size(PTB))];
end

%% 3. Solve (feasibility)
sol = optimize(constraints, [], sdpsettings('solver', 'mosek', 'verbose', 0));

if sol.problem == 0
    fprintf('System is ROBUSTLY STABLE over [%.1f, %.1f] rad/s\n', omega_bornes);
    P_opt = sdpval(P);
    figure; stem(P_opt); title('Common Lyapunov P(\theta)');    
    figure; plot(HmqNEig(A0 + omega_bornes(1)*A1, hlmi, 2*pi/omega_bornes(1)), 'o'); legend('Eigenvalues at \omega_1');
    hold on;plot(HmqNEig(A0 + omega_bornes(2)*A1, hlmi, 2*pi/omega_bornes(2)), 'o'); legend('Eigenvalues at \omega_2');
    hold off;
    figure; lsim(A0 + omega_bornes(1)*A1, 2, [], 2*pi/omega_bornes(1)); title('Response at \omega_1');
    hold on;lsim(A0 + omega_bornes(2)*A1, 2, [], 2*pi/omega_bornes(2)); title('Response at \omega_2');
    hold off;
else
    fprintf('Stability LMI INFEASIBLE: %s\n', sol.info);
end
