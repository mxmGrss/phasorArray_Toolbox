% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% ── Polytopic H2 Synthesis ──────────────────────────────────────────────
%  ẋ = [A]X + [B]U + Bw·w,    z = [C]X + [D]U
%  [·]_i = ·₀ + ωᵢ·₁        ω ∈ [ωmin, ωmax]

%% 0. Problem data
nx = 3;  nu = 2;  ny = 4;  ha = 3;
omega_bornes = [100, 1500] * 2*pi / 60;

A0 = (PhasorArray.random(nx, nx, ha));
A1 = (PhasorArray.random(nx, nx, ha));
B0 = (PhasorArray.random(nx, nu, ha));
B1 = (PhasorArray.random(nx, nu, ha));
C0 = (PhasorArray.random(ny, nx, ha));
C1 = (PhasorArray.random(ny, nx, ha));
D0 = (PhasorArray.random(ny, nu, ha));
D1 = (PhasorArray.random(ny, nu, ha));
Bw = PhasorArray(eye(nx));

hP = 10;  hlmi = 10;

%% 1. Decision variables
Y = PhasorArray.ndsdpvar(nu, nx, hP, 'PhasorType', 'full');
P = PhasorArray.ndsdpvar(nx, nx, hP);
W = PhasorArray.ndsdpvar(ny, ny, hP);

PTB  = P.T_tb(hlmi);
WTB  = W.T_tb(hlmi);
BwTB = T_tb(Bw * Bw',hlmi);

%% 2. LMIs at each vertex (products computed before lifting)
F = [PTB >= 1e-6 * eye(size(PTB))];

for ii = 1:numel(omega_bornes)
    omega_i = omega_bornes(ii);
    T_i = 2*pi / omega_i;
    Ni  = N_tb(nx, hlmi, T_i);

    Ai = A0 + omega_i * A1;
    Bi = B0 + omega_i * B1;
    Ci = C0 + omega_i * C1;
    Di = D0 + omega_i * D1;

    AiP  = T_tb((Ai * P),hlmi);
    BiY  = T_tb((Bi * Y),hlmi);
    Pdot = Ni * PTB - PTB * Ni;      % lifted d/dt P(t)

    % Stability: A P + P A^H - dot(P) - B Y - Y^H B^H + BwBw^H < 0
    F = [F; -(AiP + AiP') + Pdot + (BiY + BiY') - BwTB >= 1e-6*eye(size(PTB))];

    % H₂ performance: [W, CP-DY; (CP-DY)^H, P] >= 0
    G21 = T_tb((Ci * P - Di * Y),hlmi);
    F = [F; [WTB, G21; G21', PTB] >= 1e-6*eye(size(WTB)+size(PTB))];
end

%% 4. Solve
gamma = phas(trace(W), 0);
sol   = optimize(F, gamma, sdpsettings('solver', 'mosek', 'verbose', 1));

%% 5. Extract
if sol.problem == 0
    K = sdpval(Y) / sdpval(P);
    fprintf('Optimal H2 norm: %.4e\n', value(gamma));

    figure; stem(K); sgtitle('Gain-scheduled K(\theta)');
    T1 = 2*pi/omega_bornes(1);
    A_cl = A0 + omega_bornes(1)*A1 + (B0 + omega_bornes(1)*B1)*K;
    figure; plot(HmqNEig(A_cl, hlmi, T1), 'o'); title('Closed-loop eig at \omega_1');
    figure; lsim(A_cl, 2, [], T1); title('Closed-loop response at \omega_1');
end
