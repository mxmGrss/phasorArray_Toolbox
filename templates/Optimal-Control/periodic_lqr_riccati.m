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
A = PhasorArray(-5*eye(nx)) + PhasorArray.random(nx, nx, ha)*0.5; % strictly stable open-loop for illustration
B = PhasorArray.random(nx, nu, ha);

Q = PhasorArray(10*eye(nx));     % state weight
R = PhasorArray(eye(nu));        % control weight
K0 = PhasorArray(zeros(nu,nx));  % initial stabilising gain (zero if A is open-loop stable)

%% 1. Solve — automatic harmonic truncation growth
[K, S, htrunc, H] = RicHarmonicKlein(A, B, Q, R, K0, T, ...
    'max_iter', 100, 'residualThreshold', 1e-5, 'autoUpdateh', true);

%% 2. Inspect result
fprintf('Final truncation order: h = %d\n', htrunc);
fprintf('Gain K has %d harmonics\n', (size(K,3)-1)/2);

figure; stem(K); sgtitle('Controller gain K(\theta)');
figure; plot(HmqNEig(A - B*K, htrunc, T), 'o'); title('Closed-loop eigenvalues');
figure; lsim(A - B*K, 2, [], T); title('Closed-loop response');
