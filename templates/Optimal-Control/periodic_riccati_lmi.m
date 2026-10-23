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
P = PhasorArray.ndsdpvar(nx, nx, hP);

%% 2. Riccati inequality in Schur form, on the symbols (products are exact)
%   M(t) = [dP/dt + A'P + PA + Q,  PB + N';  (PB + N')',  R] >= 0  for all t.
%   gramConstraint imposes it exactly, whatever the degree of M; a finite
%   section T_tb(M,h) >= 0 would only be a necessary condition.
PB_N = P * B + N_ric';
M = [P.d(T) + A' * P + P * A + Q_ric, PB_N; PB_N', R_ric];

[F1, G1, c1] = M.gramConstraint(0);
[F2, G2, c2] = P.gramConstraint(0);

%% 3. Solve: maximise the average trace of P (lower bound of the Riccati solution)
sol = optimize([F1, F2], -trace(phas(P, 0)), sdpsettings('solver', 'mosek', 'verbose', 1));

%% 4. Validate the numerical certificate (Proposition 4 of the TBLMI paper)
if sol.problem == 0
    report = sdpval(M).checkGramCertificate(0, value(G1), c1);
    fprintf('min_t lambda_min(M(t)) >= %.2e\n', report.lower);
end

%% 5. Extract gain
if sol.problem == 0
    PP = PhasorArray(value(P.value));
    K  = R_ric \ (B.' * PP + N_ric);

    figure('Name','riccati_lmi — Controller gain'); stem(K); sgtitle('Controller gain K(\theta)');
    figure('Name','riccati_lmi — Closed-loop eigenvalues'); plot(HmqNEig(A - B*K, h, T), '*'); title('Closed-loop eig');
    figure('Name','riccati_lmi — Closed-loop response'); lsim(A - B*K, 2, [], T); title('Closed-loop response');
end
