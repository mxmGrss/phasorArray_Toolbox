% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% ── Periodic Observer via Iterative Riccati ─────────────────────────────
%  ẋ = A(t)x + w,    y = C(t)x + v

%% 0. Problem data
nx = 3;  ny = 2;  ha = 3;
T = 0.1;

rng(42);
A = PhasorArray(-2*eye(nx)) + PhasorArray.random(nx, nx, ha)*0.5; % strictly stable open-loop for illustration
C = (PhasorArray.random(ny, nx, ha));

W = PhasorArray(0.1 * eye(nx));    % process noise covariance
V = PhasorArray(eye(ny));          % measurement noise covariance

%% 1. Dual Riccati: solve on (A', C', W, V)
%  A_dual = A'   (transpose of A)
%  B_dual = C'   (transpose of C)
%  Q_dual = W    (process noise = "state weight" in dual)
%  R_dual = V    (measurement noise = "control weight" in dual)
L0 = PhasorArray(zeros(ny, nx));   % initial dual gain (ny × nx)

[L_dual, Sigma, info_ric] = RicHarmonicKlein(A.', C.', W, V, L0, T, ...
    'maxIter', 100, 'thresholdResidual', 1e-3, 'autoUpdateh', true);

%% 2. Observer gain: L = Σ C' V⁻¹
%  Note: RicHarmonicKlein returns K = R⁻¹B'S in the control sense
%        which corresponds to L' = V⁻¹ C Σ in the observer sense.
%  So L = L_dual' if viewing through duality.
L = L_dual.';

%% 3. Verify observer poles
fprintf('Observer truncation order: h = %d\n', info_ric.h);
figure('Name','observer_riccati — Observer eigenvalues'); plot(HmqNEig(A - L*C, info_ric.h, T), 'o'); title('Observer eigenvalues');
figure('Name','observer_riccati — Observer gain'); stem(L); sgtitle('Observer gain L(\theta)');

%% 4. Simulate estimation error dynamics
figure('Name','observer_riccati — Error dynamics'); lsim(A - L*C, 2, [], T); title('Observer error dynamics');
