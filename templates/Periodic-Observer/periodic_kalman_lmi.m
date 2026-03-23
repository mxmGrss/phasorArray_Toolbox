% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% ── Periodic Kalman Filter (Schur LMI — dual of LQR) ───────────────────
%  ẋ = A(t)x + B(t)u + w       (process noise w, covariance W)
%  y = C(t)x + v               (measurement noise v, covariance V)

%% 0. Problem data (replace with your system)
nx = 3;  ny = 2;  ha = 3;
T = 0.1;

rng(42);
A = PhasorArray(-2*eye(nx)) + PhasorArray.random(nx, nx, ha)*0.5; % strictly stable open-loop for illustration
C = (PhasorArray.random(ny, nx, ha));

% Noise covariances (periodic or constant)
W = PhasorArray(0.1 * eye(nx));    % process noise
V = PhasorArray(eye(ny));          % measurement noise

%% 1. Decision variables
hP  = 10;
h   = 10;

% Q is the state covariance matrix (symmetric, nx x nx)
Q_var = PhasorArray.ndsdpvar(nx, nx, hP, 'PhasorType', 'symmetric', 'real', true);
% Z = Q L — note: Z is nx × ny because L is nx × ny
Z     = PhasorArray.ndsdpvar(nx, ny, hP, 'PhasorType', 'full', 'real', true);

QT = Q_var.T_tb(h);
N  = N_tb(A, h, T);

%% 2. Dual Schur-complement LMI
%  Dual of the control problem: A Q + Q A' - Z C - C' Z' - dQ/dt + W
LMI_top = (A * Q_var + Q_var * A' - Z * C - C' * Z');
LMI_top_tb = LMI_top.T_tb(h) - (N*QT - QT*N);   % -dQ/dt lifting rules
LMI_top_tb = 0.5 * (LMI_top_tb + LMI_top_tb'); % Enforce numerical symmetry

ZT = Z.T_tb(h);

F = [QT >= 0];
F = [F, [...
    LMI_top_tb,  Q_var.T_tb(h),  ZT;
    Q_var.T_tb(h), -T_tb(W^-1, h), zeros((2*h+1)*nx, (2*h+1)*ny);
    ZT',            zeros((2*h+1)*ny, (2*h+1)*nx), -T_tb(V^-1, h)] <= -1e-6*eye(size(LMI_top_tb,1) + size(QT,1) + size(ZT,2))];

%% 3. Solve
obj = -trace(Q_var{:,:,0});
sol = optimize(F, obj, sdpsettings('solver', 'mosek', 'verbose', 1));

%% 4. Extract observer gain L(t) = Q⁻¹ Z
if sol.problem == 0
    Q_opt = PhasorArray(value(Q_var.value));
    Z_opt = PhasorArray(value(Z.value));
    
    L = Q_opt \ Z_opt;
    L = L.reduce('reduceMethod', 'relative', 'reduceThreshold', 1e-4);

    % Check observer poles
    figure('Name','kalman_lmi — Observer eigenvalues'); plot(HmqNEig(A - L*C, h, T), 'o'); title('Observer eigenvalues');
    figure('Name','kalman_lmi — Observer gain'); stem(L); sgtitle('Observer gain L(\theta)');
    figure('Name','kalman_lmi — Error dynamics'); lsim(A - L*C, 2, [], T); title('Observer error dynamics');
else
    warning('Observer LMI infeasible: %s', sol.info);
end
