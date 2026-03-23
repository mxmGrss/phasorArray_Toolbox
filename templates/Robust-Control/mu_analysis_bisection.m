% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% ── μ-Analysis via Bisection ────────────────────────────────────────────
%  Multiplicative uncertainty: (1+δ)(A₀(θ) + ω A₁(θ)),  ω ∈ [ωmin, ωmax]
rng(42)
%% 0. Problem data (replace with your system)
nx = 3;  ha = 4;
A0 = (PhasorArray.random(nx, nx, ha))-eye(nx)*nx*nx;
A1 = (PhasorArray.random(nx, nx, ha))*0.001;
omega_bornes = [20, 40] * 2*pi;

hP = 10;  hlmi = 10;

%% 1. Decision variables
P = PhasorArray.ndsdpvar(nx, nx, hP);
D = PhasorArray.ndsdpvar(nx, nx, hP);   % D free (δ scalar ⇒ DΔ=ΔD for any D)

PTB = P.T_tb(hlmi);
DTB = D.T_tb(hlmi);

%% 2. Build parametric LMI blocks (function of β)
LMI_nom  = cell(2,1);   % nominal Lyapunov + uncertainty coupling
D_beta   = cell(2,1);   % β·D contribution

for i = 1:2
    omega_i = omega_bornes(i);
    T_i     = 2*pi / omega_i;
    Ai      = A0 + omega_i * A1;       % full plant at vertex
    Ni      = N_tb(nx, hlmi, T_i);

    % Nominal closed-loop: A_i P and derivative NP-PN
    AiP  = T_tb((Ai * P),hlmi) - Ni * PTB;   % = (A_i,T − ωN)·P
    % Uncertainty channel: plant Toeplitz (without N)
    AiP_plant = T_tb((Ai * P),hlmi);          % = T_tb(Ai·P)

    LMI_nom{i}  = [AiP + AiP',  AiP_plant'; AiP_plant, -DTB];
    D_beta{i}   = blkdiag(DTB, zeros(size(DTB)));
end

%% 3. Bisection on β  (max |δ|² ≤ β)
tol = 0.01;
beta_lo = 0;  beta_hi = 1;
ops = sdpsettings('solver', 'mosek', 'verbose', 0);

fprintf('ready for bissection\n')
% Expand upper bound until infeasible
beta_cur = beta_hi;
while true
    F = [PTB >= 1e-12*eye(size(PTB)); DTB >= 1e-12*eye(size(DTB)); ...
         LMI_nom{1} + beta_cur * D_beta{1} <= 0; ...
         LMI_nom{2} + beta_cur * D_beta{2} <= 0;...
         trace(DTB)==hlmi];
    s = optimize(F, [], ops);
    if s.problem ~= 0 || beta_hi > 1e6, break; end
    beta_hi = beta_hi * 2;  beta_cur = beta_hi;
end

fprintf('Robustness margin upper bound: (βupper = %.4f)\n', beta_hi);

% Bisect
beta_star = 0;
iterCount =  0;
while (beta_hi - beta_lo) > tol
    iterCount = iterCount+1;
    beta_mid = (beta_hi + beta_lo) / 2;
    F = [PTB >= 1e-12*eye(size(PTB)); DTB >= 1e-12*eye(size(DTB)); ...
         LMI_nom{1} + beta_mid * D_beta{1} <= 0; ...
         LMI_nom{2} + beta_mid * D_beta{2} <= 0;...
         trace(DTB)==hlmi];
    s = optimize(F, [], ops);
    if s.problem == 0
        beta_lo   = beta_mid;
        beta_star = beta_mid;
    else
        beta_hi = beta_mid;
    end
end

fprintf('Robustness margin: |δ|_max = √β* = %.4f  (β* = %.4f)\n', sqrt(beta_star), beta_star);
