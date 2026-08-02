% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% ── H2 Synthesis (state-feedback form) ──────────────────────────────────
% dx/dt = A(t)x + B_w(t)w + B_u(t)u
% z     = C_z(t)x + D_z(t)u            (performance output)

%% 0. Problem data
nx = 3;  nw = 2; nu = 2;  nz = 3;  ha = 3;
T = 0.1;

rng(42);
A  = PhasorArray(-5*eye(nx)) + PhasorArray.random(nx, nx, ha)*0.1; % strictly stable open-loop for illustration
Bw = PhasorArray.random(nx, nw, ha);
Bu = PhasorArray.random(nx, nu, ha);
Cz = PhasorArray.random(nz, nx, ha);
Dz = PhasorArray.random(nz, nu, ha);

hP = 6;  h = 6;

%% 1. Decision variables
Q = PhasorArray.ndsdpvar(nx, nx, hP);
Y = PhasorArray.ndsdpvar(nu, nx, hP, "symmetry","real");
Z = PhasorArray.ndsdpvar(nz, nz, hP);

QT = Q.T_tb(h);
YT = Y.T_tb(h);
ZT = Z.T_tb(h);
N  = N_tb(nx, h, T);

%% 2. Build H₂ LMI
% LMI 1: A Q + Q A' + Bu Y + Y' Bu' - dQ/dt + Bw Bw' < 0
LMI1 = A*Q + Q*A' + Bu*Y + Y'*Bu' + Bw*Bw';
LMI1_tb = LMI1.T_tb(h) - (N*QT - QT*N);
LMI1_tb = 0.5 * (LMI1_tb + LMI1_tb'); % Enforce exact numerical symmetry

% LMI 2: [Z, Cz Q + Dz Y; *, Q] > 0
LMI2_12 = Cz*Q + Dz*Y;
LMI2_12_tb = LMI2_12.T_tb(h);

F = [ QT >= 0 ];
F = [ F, LMI1_tb <= -1e-6 * eye(size(LMI1_tb)) ]; 
F = [ F, [ ZT, LMI2_12_tb ; LMI2_12_tb', QT ] >= 0 ];

%% 4. Solve
% Minimize Trace(Z_0)
obj = trace(Z{:,:,0});
sol = optimize(F, obj, sdpsettings('solver', 'mosek', 'verbose', 1));

%% 5. Extract state-feedback K = Y Q⁻¹
if sol.problem == 0
    Q_val = PhasorArray(value(Q.value));
    Y_val = PhasorArray(value(Y.value));
    K = Y_val / Q_val;
    K = K.reduce('reduceMethod', 'relative', 'reduceThreshold', 1e-4);

    fprintf('Optimal H2 norm squared: %.6f\n', value(obj));
    figure('Name','h2_lmi — Controller gain'); stem(K); sgtitle('Controller gain K(\theta)');
    figure('Name','h2_lmi — Closed-loop eigenvalues'); plot(HmqNEig(A + Bu*K, h, T), '*'); title('Closed-loop eig');
    figure('Name','h2_lmi — Closed-loop response'); lsim(A + Bu*K, 2, [], T); title('Closed-loop response');
else
    warning('LMI infeasible: %s', sol.info);
end
