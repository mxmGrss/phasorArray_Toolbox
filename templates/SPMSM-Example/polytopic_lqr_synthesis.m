% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END

%% 3a. Frequency-affine decomposition
A0 = A_aug;
A1 = PhasorArray(zeros(nx_aug));     % no explicit ω-affine term for SPMSM with constant L

B0 = B_aug;
B1 = PhasorArray(zeros(nx_aug, nu));

plot(A0)   % visualise the periodic augmented dynamics

%% 3b. Operating range
omega_bornes = [100, 3000] * 2*pi / 60;   % [RPM] → [rad/s electrical]
omega_bornes = omega_bornes * param.pp;     % mechanical → electrical

%% 3c. LQR weights
Q_lqr = PhasorArray(diag([10, 10, 100, 50000]));   % penalise speed error and integrator
R_lqr = PhasorArray(eye(nu));                       % control effort

Dw = PhasorArray(eye(nx_aug));                      % disturbance input matrix

%% 3d. LMI synthesis
hP   = 10;    % Lyapunov harmonic order
hlmi = 10;    % LMI truncation order

Y = PhasorArray.ndsdpvar(nu, nx_aug, hP, 'PhasorType', 'full');
P = PhasorArray.ndsdpvar(nx_aug, nx_aug, hP);
X = PhasorArray.ndsdpvar(nu, nu, hP);
sdpvar gam;

PTB  = P.T_tb(hlmi);
DwTB = T_tb((Dw * Dw'), hlmi);

F = [PTB >= 1e-6 * eye(size(PTB))];

for ii = 1:numel(omega_bornes)
    omega_i = omega_bornes(ii);
    T_i = 2*pi / omega_i;
    Ni  = N_tb(nx_aug, hlmi, T_i);

    Ai = A0 + omega_i * A1;
    Bi = B0 + omega_i * B1;

    AiP  = T_tb((Ai * P), hlmi);
    BiY  = T_tb((Bi * Y), hlmi);
    Pdot = Ni * PTB - PTB * Ni;

    G1 = -(AiP + AiP') + Pdot + (BiY + BiY') - DwTB;
    F = [F; G1 >= 1e-6 * eye(size(G1))];
end

% Performance LMI
Rsqrt = reduce(R_lqr^(1/2));
G2 = [X, Rsqrt*Y; (Y.')*Rsqrt, P];
F  = [F; G2.T_tb(hlmi) >= 1e-6 * eye(size(G2.T_tb(hlmi)))];
F  = [F; phas(-trace(Q_lqr*P) - trace(X) + gam, 0) >= 0];

%% 3e. Solve
sol = optimize(F, gam, sdpsettings('solver', 'mosek', 'verbose', 1));

if sol.problem ~= 0
    error('LMI infeasible: %s', sol.info);
end

K = sdpval(Y) / sdpval(P);
K = K.reduce('reduceMethod', 'relative', 'reduceThreshold', 1e-4);

stem(K);   sgtitle('Gain-scheduled K(\theta)');
plot(K)    % time-domain reconstruction
plot(P)    % Lyapunov function
