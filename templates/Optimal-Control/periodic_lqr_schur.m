% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% ── Periodic LQR via Schur-Complement LMI ──────────────────────────────
% dx/dt = A(t)x + B(t)u    with A(t), B(t) periodic of period T

%% 0. Problem data (replace with your system)
nx = 3;  nu = 2;  ha = 3;
T = 0.1;

rng(42);
A = PhasorArray(-2*eye(nx)) + PhasorArray.random(nx, nx, ha)*0.5; % strictly stable open-loop for illustration
B = (PhasorArray.random(nx, nu, ha));

% LQ weights (periodic or constant)
Ql = PhasorArray(diag(10*ones(nx,1)));   % state weight
Rl = PhasorArray(diag( 1*ones(nu,1)));   % control weight

%% 1. Decision variables
hP  = 10;           % harmonic order of the Lyapunov variable
h   = 10;           % LMI truncation order

Q_var = PhasorArray.ndsdpvar(nx, nx, hP, 'PhasorType', 'symmetric', 'real', true);
Y     = PhasorArray.ndsdpvar(nu, nx, hP, 'PhasorType', 'full', 'real', true);

QT = Q_var.T_tb(h);
N  = N_tb(A, h, T);

%% 2. Schur-complement LMI
LMI_top = (A*Q_var + B*Y) + (A*Q_var + B*Y).';
LMI_top_tb = LMI_top.T_tb(h) - (N*QT - QT*N);   % -dQ/dt lifting
LMI_top_tb = 0.5 * (LMI_top_tb + LMI_top_tb'); % Enforce exact numerical symmetry

F = [QT >= 1e-4 * eye(size(QT))];     % Q > 0
F = [F, [...
    LMI_top_tb,                       Q_var.T_tb(h),         Y.T_tb(h)';
    Q_var.T_tb(h),                    -T_tb(Ql^-1, h),        zeros((2*h+1)*nx, (2*h+1)*nu);
    Y.T_tb(h),                         zeros((2*h+1)*nu, (2*h+1)*nx), -T_tb(Rl^-1, h)] <= 0];

%% 3. Solve
obj = -trace(Q_var).phas(0);
sol = optimize(F, obj, sdpsettings('solver', 'mosek', 'verbose', 1));

%% 4. Extract controller K(t) = Y(t) / Q(t)
if sol.problem == 0
    K = sdpval(Y) / sdpval(Q_var);
    K = K.reduce('reduceMethod', 'relative', 'reduceThreshold', 1e-4);

    % Verify closed-loop stability
    figure; stem(K); sgtitle('Controller gain K(\theta)');
    figure; plot(HmqNEig(A + B*K, h, T), 'o'); title('Closed-loop eigenvalues');
    figure; lsim(A + B*K, 2, [], T);           title('Closed-loop response');
else
    warning('LMI infeasible: %s', sol.info);
end
