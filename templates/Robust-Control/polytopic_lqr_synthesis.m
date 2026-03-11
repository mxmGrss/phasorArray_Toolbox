% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% ── Polytopic LQR Synthesis (frequency-affine) ─────────────────────────
%  ẋ = (A₀ + ω A₁ − ωN)X + (B₀ + ω B₁)U     ω ∈ [ωmin, ωmax]

%% 0. Problem data
nx = 3;  nu = 2;  ha = 3;
omega_bornes = [100, 1500] * 2*pi / 60;    % RPM → rad/s

A0 = (PhasorArray.random(nx, nx, ha));
A1 = (PhasorArray.random(nx, nx, ha));
B0 = (PhasorArray.random(nx, nu, ha));
B1 = (PhasorArray.random(nx, nu, ha));

Q = PhasorArray(diag([1 1 10]));    % state weight
R = PhasorArray(eye(nu));           % control weight
Dw = PhasorArray(eye(nx));          % disturbance input

hP = 10;  hlmi = 10;

%% 1. Decision variables
Y = PhasorArray.ndsdpvar(nu, nx, hP, 'PhasorType', 'full');
P = PhasorArray.ndsdpvar(nx, nx, hP);
X = PhasorArray.ndsdpvar(nu, nu, hP);

sdpvar gam;

PTB  = P.T_tb(hlmi);
DwTB = T_tb((Dw * Dw'),hlmi);

%% 2. Stability LMI at each vertex (products computed before lifting)
F = [PTB >= 1e-6 * eye(size(PTB))];

for ii = 1:numel(omega_bornes)
    omega_i = omega_bornes(ii);
    T_i = 2*pi / omega_i;
    Ni  = N_tb(nx, hlmi, T_i);

    Ai = A0 + omega_i * A1;
    Bi = B0 + omega_i * B1;

    AiP  = T_tb((Ai * P),hlmi);     % T_tb(Ai·P)
    BiY  = T_tb((Bi * Y),hlmi);     % T_tb(Bi·Y)
    Pdot = Ni * PTB - PTB * Ni;      % lifted d/dt P(t)

    G1 = -(AiP + AiP') + Pdot + (BiY + BiY') - DwTB;

    F = [F; G1 >= 1e-6 * eye(size(G1))];
end

%% 4. Performance LMI: [X R^{1/2}Y; Y'R^{1/2} P] >= 0
Rsqrt = reduce(R^(1/2));
G2 = [X, Rsqrt*Y; (Y.')*Rsqrt, P];
F  = [F; G2.T_tb(hlmi) >= 1e-6 * eye(size(G2.T_tb(hlmi)))];

% Trace constraint
F = [F; phas(-trace(Q*P) - trace(X) + gam, 0) >= 0];

%% 5. Solve
obj = gam;
sol = optimize(F, obj, sdpsettings('solver', 'mosek', 'verbose', 1));

%% 6. Extract controller
if sol.problem == 0
    K = sdpval(Y) / sdpval(P);
    fprintf('Optimal γ = %.4e\n', value(gam));

    figure; stem(K); sgtitle('K(\theta)');
    T1 = 2*pi/omega_bornes(1);
    A_cl = A0 + omega_bornes(1)*A1 + (B0 + omega_bornes(1)*B1)*K;
    figure; plot(HmqNEig(A_cl, hlmi, T1), 'o'); title('Closed-loop eig at \omega_1');
    figure; lsim(A_cl, 2, [], T1); title('Closed-loop response at \omega_1');
end
