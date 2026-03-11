% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% ── Riccati LMI (Schur Complement) ─────────────────────────────────────
nx = 2;  nu = 2;  ha = 3;
hP = 10;  h = 10;
T = 3;

rng(42);
A = PhasorArray(-2*eye(nx)) + PhasorArray.random(nx, nx, ha)*0.5; % strictly stable open-loop for illustration
B = PhasorArray.random(nx, nu, ha);

% Riccati weights
Q_ric = PhasorArray(eye(nx));
R_ric = PhasorArray(eye(nu));
N_ric = PhasorArray(zeros(nu, nx));     % cross-term (often zero)

%% 1. Decision variable
P  = PhasorArray.ndsdpvar(nx, nx, hP, 'PhasorType', 'symmetric', 'real', true);
PT = P.T_tb(h);
N  = N_tb(nx, h, T);

%% 2. Schur form of the Riccati inequality (products before lifting)
AP    = T_tb((A' * P + P * A),h);    % T_tb(A^H P + P A)
Pdot  = N * PT - PT * N;              % lifted d/dt P(t)
PB_Nc = T_tb((P * B),h) + N_ric.T_tb(h)';   % T_tb(PB) + N_c^H

F11 = AP + Pdot + Q_ric.T_tb(h);
F12 = PB_Nc;
F22 = R_ric.T_tb(h);

F = [PT >= 0, ...
     [F11, F12; F12', F22] >= 0];

%% 4. Solve
sol = optimize(F, -trace(PT), sdpsettings('solver', 'mosek', 'verbose', 1));

%% 5. Extract gain
if sol.problem == 0
    PP = PhasorArray(value(P.value));
    K  = R_ric \ (B.' * PP + N_ric);

    figure; stem(K); sgtitle('Controller gain K(\theta)');
    figure; plot(HmqNEig(A - B*K, h, T), '*'); title('Closed-loop eig');
    figure; lsim(A - B*K, 2, [], T);           title('Closed-loop response');
end
