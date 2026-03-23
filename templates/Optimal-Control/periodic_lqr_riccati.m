% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% ── Periodic LQR via Iterative Riccati ──────────────────────────────────
% dx/dt = A(t)x + B(t)u    with period T

%% 0. Problem data (replace with your system)
nx = 3;  nu = 2;  ha = 3;
T = 0.1;

rng(42);
A = PhasorArray(5*eye(nx)) + PhasorArray.random(nx, nx, ha)*0.5; % open-loop unstable
B = PhasorArray.random(nx, nu, ha)+rand(nx,nu)*5;

Q = PhasorArray(10*eye(nx));     % state weight
R = PhasorArray(eye(nu));        % control weight
K0 = [];                         % [] → DC LQR fallback (required when A is open-loop unstable)

%% 1. Solve — automatic harmonic truncation growth
[K, S, info_ric] = RicHarmonicKlein(A, B, Q, R, K0, T, ...
    'maxIter', 100, 'thresholdResidual', 1e-5, 'autoUpdateh', true,verbose=2,warmStartFraction=0.7);

%% 2. Inspect result
fprintf('Final truncation order: h = %d\n', info_ric.h);
fprintf('Gain K has %d harmonics\n', (size(K,3)-1)/2);

figure('Name','lqr_riccati — Controller gain'); stem(K); sgtitle('Controller gain K(\theta)');
figure('Name','lqr_riccati — Closed-loop eigenvalues'); plot(HmqNEig(A - B*K, info_ric.h, T), 'o'); title('Closed-loop eigenvalues');
figure('Name','lqr_riccati — Closed-loop response'); lsim(A - B*K, 2, [], T); title('Closed-loop response');
