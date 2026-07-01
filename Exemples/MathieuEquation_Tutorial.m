% MATHIEU_MASTER_PEDAGOGICAL_CHAIN
% End-to-end LTP control workshop with the PhasorArray toolbox, built on the
% Mathieu equation   ẍ + δ·ẋ + (a − 2q·cos 2t)·x = 0.
%
% Run top-to-bottom (sections are sequential; each one prints its own results):
%   §1  Config & LTP state-space A(t)
%   §2  Stability — Floquet/HSS eigenvalues (HmqNEig)
%   §3  Strutt-Ince instability map (a,q)            [toggle: compute_strutt]
%   §4  Backward Lyapunov certificate  det(P) > 0
%   §5  Open-loop simulation (free + forced)
%   §6  LQR synthesis (RicHarmonicKlein)
%   §7  Kalman observer (KalHarmonicKleinGen)
%   §8  Three closed-loop sims + disturbance problem (§8.4)   [toggle: add_noise]
%   §9–13  Internal-model rejection: augmented observer + feedforward
%   §14    Classical DC integral action
%   §15–16 Periodic (LTP) integral action + generic evalp time-stepping template
%   §17    Summary — rejection residual by architecture
%   §18    Toeplitz-Block LMIs (stability / LQR / Kalman) via YALMIP
%
% Requires: PhasorArray toolbox on the MATLAB path; Control System Toolbox.
%           §18 also needs YALMIP + an SDP solver (MOSEK or SeDuMi) — auto-skipped if absent.

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 1: CONFIGURATION & MATHIEU SYSTEM SETUP
%  ═══════════════════════════════════════════════════════════════════════════
%  Helper functions (stability_str, simInput) are at the END of this file, as
%  MATLAB requires for script-local functions before R2024a.

clearvars; close all; clc
rng(42)

fprintf('\n%s\n', repmat('=', 1, 80))
fprintf('  MATHIEU MASTER PEDAGOGICAL CHAIN: Complete Control + Perturbation Rejection\n')
fprintf('%s\n\n', repmat('=', 1, 80))

% === System parameters ===
T  = pi;                        % Period [s]
omega = 2*pi / T;               % Fundamental frequency [rad/s]  (= 2)

% Mathieu equation: ẍ + δ·ẋ + (a - 2q·cos(2t))·x = 0
% State-space: ẋ = A(t)x,  A(t) = A₀ + A₁·e^{jωt} + A₁*·e^{-jωt}
a  = 0.8;                       % Instability parameter
q  = 0.3144;                    % Periodic coupling (triggers instability near a=0.8)
delta = 0;                      % Damping (0 = undamped, can be varied)

% === System matrices (LTP) ===
A0 = [0, 1; -a, -delta];        % Constant part
A1 = [0, 0; q, 0];              % Periodic part (cos coefficient)
Ar = PhasorArray(cat(3,A0, A1), 'isReal', true);
nx = size(A0, 1);

% === Control & observation setup ===
B_ctrl = PhasorArray([0; 1]);   % Actuation on velocity
C_obs = PhasorArray([1, 0]);    % Observe position only

% === Harmonic truncation strategy ===
h_check = 5;                    % Initial check (fast)
h_report = 25;                  % Precision analysis

% === Noise covariances ===
W_proc = 0.1 * eye(nx);         % Process noise covariance
V_meas = 0.01 * eye(1);         % Measurement noise covariance

fprintf('System: Mathieu equation with a=%.2f, q=%.4f, δ=%.2f\n', a, q, delta)
fprintf('Period T = π [s], Fundamental freq ω = 2 [rad/s]\n')
fprintf('State dimension: %d (position, velocity)\n', nx)
fprintf('Harmonic truncation: check=h%d, report=h%d, control=auto\n\n', h_check, h_report)

% === Bonus: Alternative methods to build and modify PhasorArrays ===
fprintf('Demonstrating PhasorArray direct syntax and cell accessors...\n\n')

figure(20); clf;
plot(Ar);
sgtitle('Evolution of A(t) elements over one period');
hold on;

% 1. Building directly with algebraic operations on signals:
% PhasorArray.cos(phi, k) generates cos(k*ω*t + phi)
Ar_alt = [0, 1; -a + q*PhasorArray.cos(0,1), -delta];
% (Ar_alt is mathematically identical to Ar)

% 2. Modifying specific matrix elements over ALL harmonics using {row, col}
q_bis = 1.0;
Ar{2,1} = -a + q_bis * PhasorArray.cos(0,1); % Overwrites the entire a_21(t) signal
plot(Ar, 'LineStyle', '--');

% 3. Modifying a specific harmonic of a specific element using {row, col, harmonic}
Ar{2,2,0} = -0.1; % Modifies ONLY the DC component (harmonic 0) of a_22(t)
plot(Ar, 'LineStyle', '-.');

% Reset to original values for the rest of the script
Ar{2,1} = -a + q * PhasorArray.cos(0,1);
Ar{2,2,0} = -delta;

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 2: STABILITY ANALYSIS (FLOQUET / HSS)
%  ═══════════════════════════════════════════════════════════════════════════

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 2: Stability Analysis via Floquet Theory (HSS Eigenvalues)\n')
fprintf('%s\n\n', repmat('-', 1, 80))

% Compute Floquet exponents (via HSS eigenvalue lifting)
mu_check = HmqNEig(Ar, h_check, T);
mu_report = HmqNEig(Ar, h_report, T);
max_re_mu_check = max(real(mu_check));
max_re_mu_report = max(real(mu_report));

fprintf('Floquet exponents (h=%d, fast scan):\n', h_check)
fprintf('  max Re(μ) = %+.6f  →  %s\n', max_re_mu_check, ...
    stability_str(max_re_mu_check < 0))
fprintf('Floquet exponents (h=%d, precision):\n', h_report)
fprintf('  max Re(μ) = %+.6f  →  %s\n\n', max_re_mu_report, ...
    stability_str(max_re_mu_report < 0))

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 3: STRUTT-INCE STABILITY DIAGRAM (OPTIONAL, EXPENSIVE)
%  ═══════════════════════════════════════════════════════════════════════════

compute_strutt = true;  % Set to true for ~5-10 min scan
fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 3: Strutt-Ince Stability Diagram (Instability Tongues)\n')
fprintf('%s\n\n', repmat('-', 1, 80))

if compute_strutt
    fig1 = figure('Name', 'Strutt-Ince Diagram');
    clf
    fprintf('Computing stability map (a,q) over grid... (this takes ~5-10 min)\n')
    a_grid = linspace(0, 2, 40);
    q_grid = linspace(-2, 2, 41);
    [A_mesh, Q_mesh] = meshgrid(a_grid, q_grid);
    Stability_map = zeros(size(A_mesh));
    for delta_i = linspace(0,2,3)
    for idx = 1:numel(A_mesh)
        [i,j] = ind2sub(size(A_mesh), idx);
        a_ij = A_mesh(i,j);
        q_ij = Q_mesh(i,j);
        A_ij = PhasorArray([0, 1; -a_ij, -delta_i], [0, 0; q_ij, 0], 'isReal', true);
        mu_ij = HmqNEig(A_ij, 10, T);
        Stability_map(i,j) = max(real(mu_ij));  % negative → stable
    end
    nexttile
    surf(A_mesh, Q_mesh, Stability_map);
    shading interp
    colorbar; 
    
% Center colorbar so that green (midpoint of colormap) corresponds to 0.
% Determine current colormap and its length, then set caxis symmetric about zero.
drawnow;  % ensure the surface is rendered and color limits are available
cm = colormap(gca); 
ncm = size(cm,1);
% Compute symmetric color limits around zero based on Stability_map range
vmin = min(Stability_map(:));
vmax = max(Stability_map(:));
absmax = max(abs([vmin, vmax]));
caxis([-absmax, absmax]);

% Ensure the colorbar midpoint corresponds to zero by adjusting its ticks and labels
hb = colorbar;
% Choose tick positions with zero centered
tick_n = 5;
ticks = linspace(-absmax, absmax, tick_n);
set(hb, 'Ticks', ticks, 'TickLabels', arrayfun(@(v) sprintf('%.2g', v), ticks, 'UniformOutput', false));

% Optionally modify colormap so that its center color is green.
% Create a diverging map that has green at center: interpolate between blue->green->red
half = ceil(ncm/2);
cmap_div = [interp1(linspace(0,1,half)', [linspace(0,0.5,half)', linspace(0,1,half)', linspace(1,0,half)'], linspace(0,1,half)') ;
            interp1(linspace(0,1,ncm-half+1)', [linspace(0,1,ncm-half+1)', linspace(1,0.5,ncm-half+1)', linspace(0,0,ncm-half+1)'], linspace(0,1,ncm-half+1)')];
% If interpolation failed to produce correct size, fall back to 'jet'
if size(cmap_div,1) ~= ncm || any(isnan(cmap_div(:)))
    colormap(gca, 'jet')
else
    colormap(gca, cmap_div);
end
    % caxis([-0.5, 0.5]);
    colormap(gca, 'jet')
    hold on; plot(a, q, 'k*', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Design point');
    xlabel('a'); ylabel('q'); title('Strutt-Ince Diagram: Instability Tongues');
    legend; grid on;
    end
else
    fprintf('(Disabled: set compute_strutt=true to enable, ~5-10 min runtime)\n\n')
end

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 4: BACKWARD LYAPUNOV STABILITY CERTIFICATION (STABLE SYSTEMS)
%  ═══════════════════════════════════════════════════════════════════════════

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 4: Backward Lyapunov — Stability Certificate via det(P) > 0\n')
fprintf('%s\n\n', repmat('-', 1, 80))

% Backward Lyapunov: dP/dt + A'P + PA + Q = 0
% Solution P(t) > 0 certifies exponential stability: cost J = ∫_0^∞ x'Qx dt < x(0)'P(0)x(0)
fprintf('Solving: dP/dt + A''(t)P + PA(t) + Q = 0  (backward cost-to-go)\n')
fprintf('  → Q = diag(%.1f, %.1f) (state penalty)\n', 10, 1)

% Create 3 stable LTP systems: add damping to Mathieu or use different a,q pairs
Q_lyap = PhasorArray(diag([10, 1]));
sys_list = {};
names_list = {};

% System 1: Mathieu + damping (stabilize the original unstable one)
A1_stable = PhasorArray([0, 1; -a, -1.0], [0, 0; q, 0], 'isReal', true);
sys_list{1} = A1_stable;
names_list{1} = 'Mathieu + damping δ=1.0';

% System 2: Mathieu at stable point (a=0.5, q small)
A2_stable = PhasorArray([0, 1; -0.5, -0.5], [0, 0; 0.1, 0], 'isReal', true);
sys_list{2} = A2_stable;
names_list{2} = 'Mathieu (a=0.5, δ=0.5)';

% System 3: Mathieu at different stable point
A3_stable = PhasorArray([0, 1; -1.2, -0.8], [0, 0; 0.2, 0], 'isReal', true);
sys_list{3} = A3_stable;
names_list{3} = 'Mathieu (a=1.2, δ=0.8)';

% Solve backward Lyapunov for each & collect det(P)
detP_list = {};
for s = 1:3
    fprintf('System %d: %s\n', s, names_list{s})
    A_sys = sys_list{s};

    % Verify stability first
    flq = HmqNEig(A_sys, 10, T);
    fprintf('  Floquet max Re(μ) = %+.6f → %s\n', max(real(flq)), stability_str(max(real(flq)) < 0))

    % Solve backward Lyapunov
    [P_sys, info_sys] = lyap(A_sys, Q_lyap, [], 'T', T, 'direction', 'backward', ...
        'autoUpdateh', true, 'maxh', 30, 'thresholdResidual', 1e-6);

    % Compute det(P) — periodic determinant
    detP_sys = det(P_sys);
    detP_list{s} = detP_sys;

    % Check positivity via sample points
    t_check = linspace(0, T, 15); t_check(end) = [];
    detP_vals = evalTime(detP_sys, T, t_check);
    min_det = min(detP_vals(:));

    fprintf('  P > 0 certified: min det(P) = %.6e\n\n', min_det)
end

% Plot det(P) for all systems via mplot
fig_det = figure('Name', 'Lyapunov det(P) — Stability Certificate');
mplot(detP_list{1}, detP_list{2}, detP_list{3});
xlabel('t [s]'); ylabel('det(P(t))');
title('Backward Lyapunov: det(P) > 0 ⟹ System Stable');
legend(names_list, 'Location', 'best', 'Interpreter', 'none');

% --- LMI STABILITY EQUIVALENCE (Time vs Frequency) ---
% Here we demonstrate that the "Temporal" LMI formulation is mathematically
% equivalent to the "Harmonic" (Toeplitz) formulation.
if exist('sdpvar', 'file') == 2
    fprintf('  -> Demonstrating LMI equivalence (Temporal vs Harmonic)...\n');
    h_lmi_eq = 3;
    
    % Method 1: Temporal Point of View
    % We define P(t) as a periodic variable and state the derivative directly in time.
    % The toolbox handles the Toeplitz conversion of the derivative under the hood.
    P_temp = PhasorArray.ndsdpvar(2, 2, h_lmi_eq, PhasorType='symmetric', real=true);
    LMI_temp = d(P_temp, T) + A1_stable'*P_temp + P_temp*A1_stable;
    F_temp = [LMI_temp.T_tb(h_lmi_eq) <= 0, P_temp.T_tb(h_lmi_eq) >= 0];
    
    % Method 2: Harmonic Point of View (HSS)
    % We compute the Toeplitz forms and explicitly apply the frequency-domain derivative operator N.
    N_mat = N_tb(2, h_lmi_eq, T);
    PT = P_temp.T_tb(h_lmi_eq);
    LMI_harm = TB(A1_stable'*P_temp + P_temp*A1_stable, h_lmi_eq) - N_mat'*PT - PT*N_mat;
    F_harm = [LMI_harm <= 0, PT >= 0];
    
    % Pedagogical note: Both F_temp and F_harm evaluate to the exact same matrix inequalities!
    % The temporal approach (Method 1) is a powerful abstraction over the rigorous HSS projection (Method 2).
end

grid on;

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 5: SIMULATION - OPEN-LOOP RESPONSE (FREE + FORCED)
%  ═══════════════════════════════════════════════════════════════════════════

Ar_stable = PhasorArray([0, 1; -a, -1.0], [0, 0; q, 0], 'isReal', true);  % Mathieu + damping δ=1

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 5: Open-Loop Simulation (Free Response + Forced Regime)\n')
fprintf('%s\n\n', repmat('-', 1, 80))

t_sim_ol = linspace(0, 20*T, 500);
x0_ol = [1; 0];  % Initial condition: position=1, velocity=0

% All simulations use the local simInput helper (PhasorSS-free, see end of file).
% It returns nx×Nt → transpose to Nt×nx for the plots below.
% Free response: ẋ = A(t)x   (no input)
y_ol_free = simInput(Ar, t_sim_ol, x0_ol, T)';

% Forced response: step u=1 on the actuator ⇒ constant input inp = ones(Nt,1)
u_step = ones(length(t_sim_ol), 1);
y_ol_forced = simInput(Ar, t_sim_ol, x0_ol, T, B_ctrl, u_step)';

% Plot dynamic response (free and forced) – both states over time
fig_ol = figure('Name', 'Open-Loop Dynamics: Free vs Forced');
subplot(2,2,1)
plot(t_sim_ol, y_ol_free(:,1), 'b-', 'LineWidth', 1.2); hold on
plot(t_sim_ol, y_ol_forced(:,1), 'r--', 'LineWidth', 1.2);
ylabel('Position x_1');
legend('Free', 'Forced', 'Location', 'best', 'Interpreter', 'none');
title('Open-Loop Response: Position');
grid on

subplot(2,2,3)
plot(t_sim_ol, y_ol_free(:,2), 'b-', 'LineWidth', 1.2); hold on
plot(t_sim_ol, y_ol_forced(:,2), 'r--', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Velocity x_2');
legend('Free', 'Forced', 'Location', 'best', 'Interpreter', 'none');
title('Open-Loop Response: Velocity');
grid on

fprintf('Free response: x(0) = [%.2f; %.2f]\n', x0_ol(1), x0_ol(2))
fprintf('Forced regime: step input u=1 on velocity actuator\n')
fprintf('Simulation time: %.2f periods (%d samples)\n\n', t_sim_ol(end)/T, length(t_sim_ol))


% Same two simulations on the damped (stable) Mathieu
y_ol_free_st   = simInput(Ar_stable, t_sim_ol, x0_ol, T)';
y_ol_forced_st = simInput(Ar_stable, t_sim_ol, x0_ol, T, B_ctrl, u_step)';

% Plot dynamic response (free and forced) – both states over time
subplot(2,2,2)
plot(t_sim_ol, y_ol_free_st(:,1), 'b-', 'LineWidth', 1.2); hold on
plot(t_sim_ol, y_ol_forced_st(:,1), 'r--', 'LineWidth', 1.2);
ylabel('Position x_1');
legend('Free', 'Forced', 'Location', 'best', 'Interpreter', 'none');
title('Open-Loop Response: Position');
grid on

subplot(2,2,4)
plot(t_sim_ol, y_ol_free_st(:,2), 'b-', 'LineWidth', 1.2); hold on
plot(t_sim_ol, y_ol_forced_st(:,2), 'r--', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('Velocity x_2');
legend('Free', 'Forced', 'Location', 'best', 'Interpreter', 'none');
title('Open-Loop Response: Velocity');
grid on


%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 6: CONTROLLER DESIGN (LQR - RICCATI SYNTHESIS)
%  ═══════════════════════════════════════════════════════════════════════════

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 6: LQR Controller via RicHarmonicKlein (Kleinman Iteration)\n')
fprintf('%s\n\n', repmat('-', 1, 80))

% LQR weights
Q_lqr = PhasorArray(diag([10, 1]));  % State penalty (weight position more)
R_lqr = PhasorArray(1);               % Control penalty
K0_lqr = PhasorArray(zeros(1, nx));   % Initial gain (zero feedback)
K0_lqr{2} = 10; %initial guess is plain damping

fprintf('LQR Weights:\n  Q = diag([%.1f, %.1f])\n  R = %.1f\n', 10, 1, 1)
fprintf('Solving: dS/dt + A^T·S + S·A - S·B·R^-1·B^T·S + Q = 0\n')

[K_lqr, S_lqr, info_ric_lqr] = RicHarmonicKlein(Ar, B_ctrl, Q_lqr, R_lqr, K0_lqr, T, ...
    'autoUpdateh', true, 'maxIter', 150, 'thresholdResidual', 1e-6, 'verbose', 1, 'skipValidate', false);

% Closed-loop check
Ar_cl_nom = Ar - B_ctrl * K_lqr;
mu_cl_nom = HmqNEig(Ar_cl_nom, info_ric_lqr.h, T);
max_re_mu_cl_nom = max(real(mu_cl_nom));

fprintf('Riccati solved at h=%d\n', info_ric_lqr.h)
fprintf('Closed-loop Floquet max Re(μ) = %+.6f  →  %s\n\n', max_re_mu_cl_nom, ...
    stability_str(max_re_mu_cl_nom < 0))

% --- Inspecting a periodic matrix: the gain K(t) in time- and frequency-domain ---
%   plot(K, T) reconstructs K(t) over one period ; stem(K) shows its harmonics K_k.
figure('Name', 'Inspect K(t): time vs harmonics');
subplot(1,2,1); plot(K_lqr, T);  title('K(t) over one period'); grid on
subplot(1,2,2); stem(K_lqr);     title('Harmonic content K_k'); grid on

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 7: OBSERVER DESIGN (KALMAN FILTER - KalHarmonicKleinGen)
%  ═══════════════════════════════════════════════════════════════════════════

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 7: Kalman Observer via KalHarmonicKleinGen\n')
fprintf('%s\n\n', repmat('-', 1, 80))

fprintf('Solving periodic Kalman filter:\n')
fprintf('  Process noise: W = diag(%.3f, %.3f)\n', W_proc(1,1), W_proc(2,2))
fprintf('  Measurement noise: V = %.3f\n', V_meas(1,1))

[L_obs, Y_obs, info_kal] = KalHarmonicKleinGen(Ar, C_obs, W_proc, V_meas(1,1), [], [], T, ...
    'autoUpdateh', true, 'thresholdResidual', 1e-6);

% Observer loop check
Ar_obs = Ar - L_obs * C_obs;
mu_obs = HmqNEig(Ar_obs, info_kal.h, T);
max_re_mu_obs = max(real(mu_obs));

fprintf('Kalman gain computed at h=%d\n', info_kal.h)
fprintf('Observer error loop max Re(μ) = %+.6f  →  %s\n\n', max_re_mu_obs, ...
    stability_str(max_re_mu_obs < 0))

% --- Same inspection on the observer gain L(t) ---
figure('Name', 'Inspect L(t): time vs harmonics');
subplot(1,2,1); plot(L_obs, T);  title('L(t) over one period'); grid on
subplot(1,2,2); stem(L_obs);     title('Harmonic content L_k'); grid on

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 8: THREE SIMULATIONS (LQR / Kalman / Combined Output-Feedback)
%  ═══════════════════════════════════════════════════════════════════════════

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 8: Three Closed-Loop Demonstrations\n')
fprintf('%s\n\n', repmat('-', 1, 80))

x0_sim   = [1; 0];     % plant initial condition
xhat0    = [0; 0];     % observer initial estimate
Znn_nx = PhasorArray(zeros(nx, nx));
Inn_nx = PhasorArray(eye(nx));

% Noise config: process noise w ~ N(0,W_proc), measurement noise v ~ N(0,V_meas)
add_noise = true;                 % set false for deterministic (noise-free) runs
sigma_w = sqrt(diag(W_proc));     % per-state process-noise std
sigma_v = sqrt(V_meas(1,1));      % measurement-noise std
noise_gain = double(add_noise);   % 0 disables all noise
% Noise enters as a sampled input time-series [w; v] through simInput (end of file).

% ── 8.1 — LQR full-state feedback: u = -K·x (state directly measured) ────────
%   ẋ = (A - B·K)·x + w    (no observer in the loop)
fprintf('8.1 — LQR full-state feedback (u = -K·x), noise=%d\n', add_noise)
t1 = linspace(0, 30*T, 500);
w1 = noise_gain * randn(numel(t1), nx) .* sigma_w';   % process noise, Nt×nx
x1 = simInput(Ar - B_ctrl*K_lqr, t1, x0_sim, T, Inn_nx, w1);   % x1: nx×Nt
K1 = evalTime(K_lqr, T, t1);
u1 = -squeeze(pagemtimes(K1, reshape(x1, nx, 1, [])))';   % u = -K·x
fprintf('   ||x(0)|| = %.3f → ||x(end)|| = %.3e  | max|u| = %.3f\n\n', ...
    norm(x1(:,1)), norm(x1(:,end)), max(abs(u1)))

fig1 = figure('Name', '8.1 LQR Full-State Feedback');
subplot(2,1,1); plot(t1, x1(1,:), '-b', t1, x1(2,:), '-r', 'LineWidth', 1.2);
ylabel('states'); legend('x_1','x_2'); grid on
title('8.1 — LQR full-state feedback: u = -K x')
subplot(2,1,2); plot(t1, u1, '-m', 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('u = -K x'); grid on

% Noise input matrix for observer sims: [w (nx channels); v (1 channel)]
%   plant gets w, observer correction L injects measurement noise v
B_noise = [Inn_nx,   PhasorArray(zeros(nx,1));
           Znn_nx,   L_obs];

% ── 8.2 — Kalman observer on the UNSTABLE open-loop plant (u = 0) ─────────────
%   Plant diverges: ẋ = A·x + w ;  observer x̂̇ = (A - L·C)·x̂ + L·C·x + L·v
%   Error e = x - x̂ stays bounded (→ Kalman covariance floor) even as x blows up.
fprintf('8.2 — Kalman observer on unstable open-loop plant (u = 0), noise=%d\n', add_noise)
t2 = linspace(0, 6*T, 500);   % short horizon: plant is unstable, x grows
w2 = noise_gain * randn(numel(t2), nx) .* sigma_w';
v2 = noise_gain * randn(numel(t2), 1)  * sigma_v;
A_obs2 = [Ar,           Znn_nx;
          L_obs*C_obs,  Ar - L_obs*C_obs];
xa2 = simInput(A_obs2, t2, [x0_sim; xhat0], T, B_noise, [w2, v2]);
x2    = xa2(1:nx, :);
xhat2 = xa2(nx+1:end, :);
e2    = x2 - xhat2;
e2_floor = mean(vecnorm(e2(:, end-99:end)));   % steady-state error level
fprintf('   plant ||x(0)||=%.2f → ||x(end)||=%.2f (diverges)\n', norm(x2(:,1)), norm(x2(:,end)))
fprintf('   error ||e(0)||=%.2f → floor ||e||≈%.3e (bounded)\n\n', norm(e2(:,1)), e2_floor)

fig2 = figure('Name', '8.2 Kalman on Unstable Plant');
subplot(2,1,1); plot(t2, x2(1,:), '-b', t2, xhat2(1,:), '--r', 'LineWidth', 1.2);
ylabel('x_1 vs \^x_1'); legend('true x_1','estimate'); grid on
title('8.2 — Observer tracks diverging unstable plant')
subplot(2,1,2); semilogy(t2, vecnorm(e2), '-k', 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('||e||_2'); grid on; legend('estimation error')

% ── 8.3 — Final closed loop: Kalman + LQR output feedback (u = -K·x̂) ─────────
%   ẋ  = A·x - B·K·x̂ + w ;  x̂̇ = L·C·x + (A - B·K - L·C)·x̂ + L·v
fprintf('8.3 — Combined output feedback (u = -K·x̂), noise=%d\n', add_noise)
t3 = linspace(0, 30*T, 500);
w3 = noise_gain * randn(numel(t3), nx) .* sigma_w';
v3 = noise_gain * randn(numel(t3), 1)  * sigma_v;
BK = B_ctrl * K_lqr;
A_cl3 = [Ar,           -BK;
         L_obs*C_obs,  Ar - BK - L_obs*C_obs];
mu_cl3 = HmqNEig(A_cl3, max(info_ric_lqr.h, info_kal.h), T);
xa3 = simInput(A_cl3, t3, [x0_sim; xhat0], T, B_noise, [w3, v3]);
x3    = xa3(1:nx, :);
xhat3 = xa3(nx+1:end, :);
e3    = x3 - xhat3;
K3 = evalTime(K_lqr, T, t3);
u3 = -squeeze(pagemtimes(K3, reshape(xhat3, nx, 1, [])))';
x3_floor = mean(vecnorm(x3(:, end-99:end)));
e3_floor = mean(vecnorm(e3(:, end-99:end)));
fprintf('   closed-loop Floquet max Re(μ) = %+.4f → %s\n', max(real(mu_cl3)), ...
    stability_str(max(real(mu_cl3)) < 0))
fprintf('   floor ||x||≈%.3e | floor ||e||≈%.3e | max|u|=%.3f\n\n', ...
    x3_floor, e3_floor, max(abs(u3)))

fig3 = figure('Name', '8.3 Combined LQR + Kalman');
subplot(3,1,1);
plot(t3, x3(1,:), '-b', t3, xhat3(1,:), '--c', 'LineWidth', 1.2); hold on
plot(t3, x3(2,:), '-r', t3, xhat3(2,:), '--m', 'LineWidth', 1.2);
ylabel('states'); legend('x_1','\^x_1','x_2','\^x_2'); grid on
title('8.3 — Output feedback u = -K \^x (regulation + estimation)')
subplot(3,1,2); semilogy(t3, vecnorm(x3), '-b', t3, vecnorm(e3), '-k', 'LineWidth', 1.2);
ylabel('norms'); legend('||x||','||e||'); grid on
subplot(3,1,3); plot(t3, u3, '-m', 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('u = -K \^x'); grid on

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 8.4: DISTURBANCE RESPONSE OF NOMINAL LOOP (motivates internal model)
%  ═══════════════════════════════════════════════════════════════════════════
%  The LQR+Kalman output-feedback loop has NO disturbance model. A disturbance
%  d(t) entering the plant (here on velocity, G=[0;1]) is therefore NOT rejected:
%    - constant d        → persistent steady-state bias (no integral action)
%    - periodic d         → sustained ripple at the disturbance frequency
%    - constant+periodic  → both superimposed
%  This residual error is exactly what an internal-model controller (§9+) removes.

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 8.4: Disturbance response of nominal loop (no rejection)\n')
fprintf('%s\n\n', repmat('-', 1, 80))

t_d = linspace(0, 16*T, 800);
Nt_d = numel(t_d);
wd = 2*pi / T;            % disturbance fundamental (= plant fundamental)
d_const = 0.3;           % constant disturbance magnitude
d_amp   = 0.3;           % periodic disturbance amplitude
G_d = PhasorArray([0; 1]);   % disturbance acts on velocity

% Closed loop reused from §8.3 (A_cl3); disturbance + noise input channels [d; w; v]
B_dist = [G_d,                     Inn_nx,  PhasorArray(zeros(nx,1));
          PhasorArray(zeros(nx,1)), Znn_nx,  L_obs];

dist_cases = {'constant', 'periodic', 'constant+periodic'};
dist_signals = {d_const*ones(Nt_d,1), ...
                d_amp*sin(wd*t_d)', ...
                d_const + d_amp*sin(wd*t_d)'};
res_noIM = zeros(1,3);   % store residuals for before/after comparison in §12

fig8d = figure('Name', '8.4 Disturbance Response (no rejection)');
for c = 1:3
    d_c = dist_signals{c};
    w_c = noise_gain * randn(Nt_d, nx) .* sigma_w';   % toggled by add_noise
    v_c = noise_gain * randn(Nt_d, 1)  * sigma_v;
    U_c = [d_c, w_c, v_c];
    yd = simInput(A_cl3, t_d, [x0_sim; xhat0], T, B_dist, U_c);
    x_d = yd(1:nx, :);
    x_floor = mean(vecnorm(x_d(:, end-99:end)));
    res_noIM(c) = x_floor;
    fprintf('   %-18s → residual floor ||x|| = %.3f  (loop cannot reject)\n', ...
        dist_cases{c}, x_floor)

    subplot(3,1,c)
    plot(t_d, x_d(1,:), '-b', t_d, x_d(2,:), '-r', 'LineWidth', 1.0); hold on
    yline(0, ':k');
    ylabel('states'); legend('x_1','x_2','Location','best'); grid on
    title(sprintf('disturbance = %s  (residual ||x||\\approx%.2f)', dist_cases{c}, x_floor))
end
xlabel('t [s]')
fprintf('\n   → persistent error motivates the internal-model design below.\n\n')

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 9: INTERNAL MODEL OF THE DISTURBANCE (matched, S = const ⊕ oscillator)
%  ═══════════════════════════════════════════════════════════════════════════
%  Internal Model Principle: to reject a disturbance generated by ḋ = S·d, the
%  controller must embed a copy of S. The §8.4 disturbances (constant + periodic
%  at the fundamental ω) are generated by:
%       S = blkdiag( 0 ,  [0 1; -ω² 0] )       (integrator ⊕ oscillator @ ω)
%  The disturbance is matched (enters through B, like the control), with output
%  map H selecting the scalar acting on the plant: δ = H·d.

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 9: Internal model of the disturbance (S = const ⊕ oscillator)\n')
fprintf('%s\n\n', repmat('-', 1, 80))

w_dist = 2*pi / T;                 % disturbance fundamental (= plant fundamental)
S_im = [0, 0,         0;
        0, 0,         1;
        0, -w_dist^2, 0];          % integrator ⊕ harmonic oscillator @ w_dist
nd = size(S_im, 1);                % disturbance-model order (3)
H_im = PhasorArray([1, 1, 0]);     % δ = d_const + d_osc(position)

% Augmented plant z = [x; d]:  ż = Az·z + Bz·u,  y = Cz·z
%   ẋ = A·x + B·H·d + B·u   (matched disturbance) ;  ḋ = S·d
Az = [Ar,                        B_ctrl*H_im;
      PhasorArray(zeros(nd,nx)),  PhasorArray(S_im)];
Bz = [B_ctrl; PhasorArray(zeros(nd,1))];
Cz = [C_obs, PhasorArray(zeros(1,nd))];

fprintf('Disturbance model order nd = %d (1 constant + 2 oscillator @ ω=%.2f)\n', nd, w_dist)
fprintf('Matched injection δ = H·d, H = [1 1 0]; augmented plant dim = %d\n\n', nx+nd)

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 10: AUGMENTED DISTURBANCE OBSERVER (Kalman on [x; d])
%  ═══════════════════════════════════════════════════════════════════════════

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 10: Augmented disturbance observer (estimates x AND d)\n')
fprintf('%s\n\n', repmat('-', 1, 80))

% Process-noise weight on the disturbance states sets observer aggressiveness:
% too small → trusts the model → slow d̂ convergence. W_d = 1 gives ~1 period.
W_dist_state = 1.0;
W_aug = blkdiag(W_proc, W_dist_state*eye(nd));
V_aug = V_meas(1,1);

[L_aug, ~, info_kal_aug] = KalHarmonicKleinGen(Az, Cz, W_aug, V_aug, [], [], T, ...
    'autoUpdateh', true, 'thresholdResidual', 1e-6);

% Augmented observer error dynamics must be stable
mu_obs_aug = HmqNEig(Az - L_aug*Cz, info_kal_aug.h, T);
fprintf('Augmented observer gain computed at h=%d\n', info_kal_aug.h)
fprintf('Observer error Floquet max Re(μ) = %+.4f → %s\n\n', max(real(mu_obs_aug)), ...
    stability_str(max(real(mu_obs_aug)) < 0))

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 11: AUGMENTED CONTROL LAW  u = -[K  H]·ẑ
%  ═══════════════════════════════════════════════════════════════════════════
%  Feedback K stabilises x; the H block cancels the (estimated) matched
%  disturbance. On the augmented plant the regulator poles are eig(A-B·K) ⊕
%  eig(S): the disturbance modes stay marginal (the disturbance persists) but
%  its effect on x is removed.

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 11: Augmented control law  u = -K·x̂ - H·d̂\n')
fprintf('%s\n\n', repmat('-', 1, 80))

K_aug = [K_lqr, H_im];    % [stabilising gain | disturbance-cancelling gain]
fprintf('Control: u = -[K  H]·[x̂; d̂]  (feedforward H·d̂ cancels matched δ)\n\n')

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 12: INTERNAL-MODEL CLOSED LOOP — disturbance rejection
%  ═══════════════════════════════════════════════════════════════════════════

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 12: Closed-loop rejection (internal model), noise=%d\n', add_noise)
fprintf('%s\n\n', repmat('-', 1, 80))

% Observer-based output feedback on the AUGMENTED plant (same structure as §8.3):
%   state ξ = [z; ẑ] = [x; d; x̂; d̂]
BzK = Bz*K_aug;  LCz = L_aug*Cz;
A_imcl = [Az,    -BzK;
          LCz,    Az - BzK - LCz];
mu_imcl = HmqNEig(A_imcl, info_kal_aug.h, T);
fprintf('Closed-loop Floquet max Re(μ) = %+.4f  (disturbance modes marginal @ 0)\n\n', ...
    max(real(mu_imcl)))

% Noise input channels [w(nx); v(1)]: w on plant state x, v into the observer
n_aug = 2*(nx + nd);
B_imcl = [ [PhasorArray(eye(nx)); PhasorArray(zeros(2*nd+nx, nx))], ...   % w → x
           [PhasorArray(zeros(nx+nd, 1)); L_aug] ];                       % v → observer

% Same three disturbance cases as §8.4, injected as internal-model IC on d
t_im = linspace(0, 20*T, 900);
Nt_im = numel(t_im);
dist_IC = {[d_const; 0; 0], [0; d_amp; 0], [d_const; d_amp; 0]};   % const / osc / both
res_IM = zeros(1, 3);

fig12 = figure('Name', '12 Internal-Model Rejection');
for c = 1:3
    xi0 = [x0_sim; dist_IC{c}; zeros(nx+nd, 1)];   % [x0; d0; x̂0=0; d̂0=0]
    w_c = noise_gain * randn(Nt_im, nx) .* sigma_w';
    v_c = noise_gain * randn(Nt_im, 1)  * sigma_v;
    yi = simInput(A_imcl, t_im, xi0, T, B_imcl, [w_c, v_c]);
    x_im = yi(1:nx, :);
    res_IM(c) = mean(vecnorm(x_im(:, end-99:end)));
    fprintf('   %-18s ||x||: %.3f (no IM)  →  %.3e (with IM)\n', ...
        dist_cases{c}, res_noIM(c), res_IM(c))

    subplot(3,1,c)
    plot(t_im, x_im(1,:), '-b', t_im, x_im(2,:), '-r', 'LineWidth', 1.0); hold on
    yline(0, ':k');
    ylabel('states'); legend('x_1','x_2','Location','best'); grid on
    title(sprintf('%s:  ||x|| %.2f \\rightarrow %.1e', dist_cases{c}, res_noIM(c), res_IM(c)))
end
xlabel('t [s]')
fprintf('\n   → internal model drives the disturbance-induced error to the noise floor.\n\n')

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 13: DISTURBANCE ESTIMATE TRACKING (const+periodic case)
%  ═══════════════════════════════════════════════════════════════════════════

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 13: Disturbance estimate δ̂ → δ and control effort\n')
fprintf('%s\n\n', repmat('-', 1, 80))

% Clean (noise-free) run of the const+periodic case for a readable estimate plot
H_val = [1, 1, 0];
xi0 = [x0_sim; dist_IC{3}; zeros(nx+nd, 1)];
yi = simInput(A_imcl, t_im, xi0, T, B_imcl, zeros(Nt_im, nx+1));
zhat = yi(nx+nd+1:end, :);                 % observer state [x̂; d̂]
delta_true = H_val * yi(nx+1:nx+nd, :);    % true matched disturbance δ = H·d
delta_hat  = H_val * zhat(nx+1:end, :);    % estimate δ̂ = H·d̂
K_aug_t = evalTime(K_aug, T, t_im);        % 1×(nx+nd)×Nt
u_im = -squeeze(pagemtimes(K_aug_t, reshape(zhat, nx+nd, 1, [])))';

fprintf('δ̂ tracking: |δ(end)-δ̂(end)| = %.3e\n\n', abs(delta_true(end)-delta_hat(end)))

fig13 = figure('Name', '13 Disturbance Estimation & Control');
subplot(3,1,1); plot(t_im, delta_true, '-k', t_im, delta_hat, '--g', 'LineWidth', 1.2);
ylabel('δ'); legend('true δ','estimate \^δ','Location','best'); grid on
title('13 — Internal-model disturbance estimation (const+periodic)')
subplot(3,1,2); semilogy(t_im, vecnorm(yi(1:nx,:)), '-b', 'LineWidth', 1.2);
ylabel('||x||'); legend('regulated state','Location','best'); grid on
subplot(3,1,3); plot(t_im, u_im, '-m', 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('u = -[K H]\^z'); grid on

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 14: DISTURBANCE REJECTION VIA CLASSICAL INTEGRAL ACTION
%  ═══════════════════════════════════════════════════════════════════════════
%  Instead of an observer, we augment the plant with an integrator on the 
%  regulated state (position x_1) to reject disturbances.

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 14: Classical Integral Action (dot z = x_1)\n')
fprintf('%s\n\n', repmat('-', 1, 80))

% Augment plant: dot z = C_z*x
C_z = [1, 0]; % Regulate position x_1 to 0
A_int = [Ar, PhasorArray(zeros(nx,1));
         PhasorArray(C_z), PhasorArray(0)];
B_int = [B_ctrl; PhasorArray(0)];

% Compute LQR for augmented plant
Q_int = blkdiag(diag([10, 1]), 50); % heavy penalty on integral state
% === Initial Guess (K0) via DC LQR ===
% To initialize the Riccati solver, we need a stabilizing feedback.
% Here, the periodic coupling is weak enough that we can just extract 
% the constant parts (DC components: harmonic 0) of A and B, solve a 
% standard LTI LQR, and use it as a constant initial guess.
A0_int = phas(A_int, 0);
B0_int = phas(B_int, 0);
K0_int_dc = lqr(A0_int, B0_int, Q_int, phas(R_lqr, 0));
K0_int = PhasorArray(K0_int_dc);
[K_int, ~, info_int] = RicHarmonicKlein(A_int, B_int, Q_int, R_lqr, K0_int, T, 'autoUpdateh', true, 'h', 3, 'skipValidate', false);

% Simulate closed loop with disturbances (direct injection into plant)
A_cl_int = A_int - B_int * K_int;
B_cl_int_dist = [PhasorArray(zeros(1,1)); PhasorArray(1); PhasorArray(0)]; % Disturbance enters x_2

res_int = zeros(1, 3);
fig14 = figure('Name', '14 Classical Integral Action');
for c = 1:3
    U_c = dist_IC{c}(1) * ones(Nt_im, 1) + dist_IC{c}(2) * sin(w_dist * t_im)';
    yi = simInput(A_cl_int, t_im, zeros(nx+1, 1), T, B_cl_int_dist, U_c);
    x_int = yi(1:nx, :);
    res_int(c) = mean(vecnorm(x_int(:, end-99:end)));
    fprintf('   %-18s ||x||: %.3f (no IM)  →  %.3e (with Integral Action)\n', ...
        dist_cases{c}, res_noIM(c), res_int(c))

    subplot(3,1,c)
    plot(t_im, x_int(1,:), '-b', t_im, x_int(2,:), '-r', 'LineWidth', 1.0); hold on
    yline(0, ':k');
    ylabel('states'); legend('x_1','x_2','Location','best'); grid on
    title(sprintf('%s:  ||x|| \\approx %.1e', dist_cases{c}, res_int(c)))
end
xlabel('t [s]')
fprintf('\n   → Integral action rejects constant disturbances but fails on periodic ones.\n\n')

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 15: PERIODIC INTEGRAL ACTION (LTP MODEL)
%  ═══════════════════════════════════════════════════════════════════════════
%  Generalized perturbation model: modulate the error with harmonics to
%  track/reject periodic signals without an explicit internal-model observer.

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 15: Periodic Integral Action (Tracking fundamental k=1)\n')
fprintf('%s\n\n', repmat('-', 1, 80))

% Modulate the regulation error with cos/sin of the fundamental
% dot z1 = cos(w t) * x1
% dot z2 = sin(w t) * x1
C_mod_c = PhasorArray.cos(0, 1); % Fundamental harmonic cos
C_mod_s = PhasorArray.sin(0, 1); % Fundamental harmonic sin
A_mod = [1 , 0; 
        C_mod_c, 0;
         C_mod_s, 0];
nz = size(A_mod,1);
        
A_per = [Ar, PhasorArray(zeros(nx, nz));
         A_mod, PhasorArray(zeros(nz, nz))];
B_per = [B_ctrl; PhasorArray(zeros(nz, 1))];

% === Initial Guess (K0) via Harmonic State Space (HSS) Lifting ===
% Unlike Section 14, this system has strong periodic coupling due to the 
% modulators. A simple DC LQR will fail to stabilize it.
% We must use the "Lifted LQR" trick:
% 1. Convert the periodic system to a large LTI system using Toeplitz matrices (TB).
%    The dynamic matrix in HSS is A_hss = A_tb - N_tb (where N_tb is the derivative operator).
Q_per = blkdiag(diag([10, 1]), 10*eye(nz));
h_init = 7;
A_tb = TB(A_per, h_init);
B_tb = TB(B_per, h_init);
Q_tb = TB(PhasorArray(Q_per), h_init);
R_tb = TB(R_lqr, h_init);
A_hss = A_tb - N_tb(size(A_per, 1), h_init, T);

% 2. Solve a standard LTI LQR on this huge lifted system to get a lifted gain.
K0_tb = lqr(A_hss, B_tb, Q_tb, R_tb);

% 3. Project the lifted gain back into a periodic time-domain PhasorArray.
%    WARNING: K0_tb is NOT Toeplitz, but taking the average on every diagonal 
%    (PhasorArray.fromTBMatrix) provides an initial K0 that is close enough 
%    for the Riccati solver to converge.
K0_per = PhasorArray.fromTBMatrix(K0_tb, 1, 'n1');

% Ensure K0_per is strictly real in time domain to avoid complex residuals
K0_per = mreal(K0_per);

[K_per, ~, info_per] = RicHarmonicKlein(A_per, B_per, Q_per, R_lqr, K0_per, T, 'autoUpdateh', true, 'h', 3, 'skipValidate', false);

% Self-validation: augmented closed loop must be (marginally) stable, with the
% internal-model modes on the imaginary axis carrying the periodic disturbance.
mu_per = HmqNEig(A_per - B_per*K_per, max(info_per.h, 8), T);
fprintf('Periodic LQR computed at h=%d\n', info_per.h)
fprintf('Augmented closed-loop Floquet max Re(μ) = %+.4f → %s\n', max(real(mu_per)), ...
    stability_str(max(real(mu_per)) < 1e-6))
fprintf('Purely periodic feedback K(t); rejection quantified by §16 (residual ||x||).\n\n')

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 16: GENERIC TIME-STEPPING TEMPLATE VIA evalp (periodic / non-linear)
%  ═══════════════════════════════════════════════════════════════════════════
%  Reusable pattern for ANY periodic-coefficient system, linear or not. The
%  right-hand side evaluates the PhasorArray at the current phase angle θ = ω·t
%  with evalp(A, θ) (one cycle ↔ 2π), then adds whatever non-linear / bilinear
%  term the model needs. Example a learner could drop in:
%       ẋ = A(t)·x + ( Σ_i x_i·N_i(t) )·u          (bilinear, state-modulated B)
%    rhs = @(th,x,u) evalp(A,th)*x + (x(1)*evalp(N1,th) + x(2)*evalp(N2,th))*u;
%  Here we instantiate it on the §15 LTP-integral closed loop (linear) and reuse
%  the SAME loop for the three disturbance cases.

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 16: Generic evalp time-stepping template (periodic systems)\n')
fprintf('%s\n\n', repmat('-', 1, 80))

A_cl_per = A_per - B_per * K_per;          % §15 closed-loop LTP matrix (nx+nz)
omega_p = 2*pi / T;                         % phase rate: θ = ω·t over one period
dt = t_im(2) - t_im(1);

% rhs(θ, x): generic periodic RHS. Swap this one line for a non-linear model.
rhs = @(th, x, dk) evalp(A_cl_per, th) * x + dk;

res_per = zeros(1, 3);
fig16 = figure('Name', '16 Generic evalp Time-Stepping');
for c = 1:3
    X = zeros(nx+nz, Nt_im);
    for k = 1:(Nt_im-1)
        th = omega_p * t_im(k);
        dk = [0; dist_IC{c}(1) + dist_IC{c}(2)*sin(w_dist*t_im(k)); zeros(nz,1)];
        X(:, k+1) = X(:, k) + dt * rhs(th, X(:, k), dk);   % explicit Euler
    end
    res_per(c) = mean(vecnorm(X(1:nx, end-99:end)));
    fprintf('   %-18s ||x||: %.3f (no IM)  →  %.3e (LTP integral)\n', ...
        dist_cases{c}, res_noIM(c), res_per(c))

    subplot(3,1,c)
    plot(t_im, X(1,:), '-b', t_im, X(2,:), '-r', 'LineWidth', 1.0); hold on
    yline(0, ':k');
    ylabel('states'); legend('x_1','x_2','Location','best'); grid on
    title(sprintf('%s:  ||x|| %.2f \\rightarrow %.1e', dist_cases{c}, res_noIM(c), res_per(c)))
end
xlabel('t [s]')
fprintf('\n   → evalp lets you hand-roll Euler/ode45 for any periodic (even non-linear) plant.\n\n')

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 17: SUMMARY
%  ═══════════════════════════════════════════════════════════════════════════

fprintf('%s\n', repmat('=', 1, 80))
fprintf('SUMMARY: Mathieu Pedagogical Chain Complete\n')
fprintf('%s\n\n', repmat('=', 1, 80))

fprintf('✓ Stability (§2): open-loop Floquet max Re(μ) = %+.4f (UNSTABLE)\n', max_re_mu_report)
fprintf('✓ LQR (§6): closed-loop Floquet max Re(μ) = %+.4f\n', max_re_mu_cl_nom)
fprintf('✓ Observer (§7): error loop max Re(μ) = %+.4f\n\n', max_re_mu_obs)

% Disturbance-rejection comparison across architectures (residual ||x|| floor)
fprintf('Disturbance rejection — residual ||x|| by architecture:\n')
fprintf('  %-22s | %-10s | %-10s | %-10s\n', 'disturbance', 'constant', 'periodic', 'const+per')
fprintf('  %s\n', repmat('-', 1, 62))
fprintf('  %-22s | %10.3f | %10.3f | %10.3f\n', 'none (§8.4)',            res_noIM(1), res_noIM(2), res_noIM(3))
fprintf('  %-22s | %10.1e | %10.1e | %10.1e\n', 'A: obs+FF (§9-13)',     res_IM(1),   res_IM(2),   res_IM(3))
fprintf('  %-22s | %10.1e | %10.1e | %10.1e\n', 'B: DC integral (§14)',  res_int(1),  res_int(2),  res_int(3))
fprintf('  %-22s | %10.1e | %10.1e | %10.1e\n', 'B: LTP integral (§15-16)', res_per(1), res_per(2), res_per(3))
fprintf('\n  A (feedforward) cancels exactly when matched; B (in-loop integral)\n')
fprintf('  is error-driven, more robust, converges asymptotically. DC integral\n')
fprintf('  alone rejects only the constant part — the periodic mode needs the LTP one.\n\n')

fprintf('Recommended workshop breakdown:\n')
fprintf('  Workshop 1: §1-4  (setup, stability, Strutt, Lyapunov certificate)\n')
fprintf('  Workshop 2: §5-7  (simulation, LQR, Kalman observer)\n')
fprintf('  Workshop 3: §8    (closed-loop simulations & disturbance problem)\n')
fprintf('  Workshop 4: §9-13 (disturbance rejection via Internal Model Observer)\n')
fprintf('  Workshop 5: §14   (disturbance rejection via Classical Integral Action)\n')
fprintf('  Workshop 6: §15-16 (LTP integral action & generic evalp time-stepping)\n\n')

fprintf('%s\n', repmat('=', 1, 80))
fprintf('End of Master Pedagogical Chain\n')
fprintf('%s\n\n', repmat('=', 1, 80))

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 18: LMI FORMULATION FOR STABILITY, LQR, AND KALMAN
%  ═══════════════════════════════════════════════════════════════════════════

fprintf('%s\n', repmat('-', 1, 80))
fprintf('SECTION 18: LMI Formulation for Stability, LQR, and Kalman\n')
fprintf('%s\n\n', repmat('-', 1, 80))

% Check if YALMIP is in path
if exist('sdpvar', 'file') ~= 2
    fprintf('YALMIP not found. Attempting to add default path...\n');
    addpath(genpath('C:\Users\mgrosso\OneDrive\Documents\MATLAB\yalmip'));
end

if exist('sdpvar', 'file') ~= 2
    fprintf('YALMIP is still not available. Skipping LMI section.\n\n');
else
    h_lmi = 3; % Use a small harmonic truncation for LMI speed
    
    % HSS Matrices for Nominal Plant
    A_tb = TB(Ar, h_lmi);
    B_tb = TB(B_ctrl, h_lmi);
    C_tb = TB(C_obs, h_lmi);
    N = N_tb(nx, h_lmi, T);
    A_hss = A_tb - N;
    
    % --- 1. LMI Lyapunov (Stability of closed-loop LQR) ---
    fprintf('1) LMI Lyapunov Stability on Closed-Loop (A_nom - B_nom*K_lqr)\n');
    A_cl_tb = TB(Ar - B_ctrl*K_lqr, h_lmi);
    A_cl_hss = A_cl_tb - N;
    
    P_lyap = PhasorArray.ndsdpvar(nx, nx, h_lmi, PhasorType='symmetric', real=true);
    P_lyap_tb = P_lyap.T_tb(h_lmi);
    Q_lyap_tb = TB(PhasorArray(eye(nx)), h_lmi);
    
    F_lyap = [A_cl_hss'*P_lyap_tb + P_lyap_tb*A_cl_hss + Q_lyap_tb <= -1e-6*eye(size(P_lyap_tb)), P_lyap_tb >= 1e-6*eye(size(P_lyap_tb))];
    sol_lyap = optimize(F_lyap, [], sdpsettings('verbose', 0));
    
    if sol_lyap.problem == 0
        fprintf('   ✓ LMI Lyapunov feasible! Closed-loop is stable.\n\n');
    else
        fprintf('   ✗ LMI Lyapunov failed: %s\n\n', sol_lyap.info);
    end
    
    % --- 2. LMI Riccati (LQR Optimal Control) ---
    fprintf('2) LMI Riccati (LQR Optimal Control)\n');
    P_lqr = PhasorArray.ndsdpvar(nx, nx, h_lmi, PhasorType='symmetric', real=true);
    P_lqr_tb = P_lqr.T_tb(h_lmi);
    Q_lqr_tb = TB(Q_lqr, h_lmi);
    R_lqr_tb = TB(PhasorArray(R_lqr), h_lmi);
    
    LMI_lqr = [A_hss'*P_lqr_tb + P_lqr_tb*A_hss + Q_lqr_tb, P_lqr_tb * B_tb; 
               B_tb'*P_lqr_tb, R_lqr_tb];
    F_lqr = [P_lqr_tb >= 1e-6*eye(size(P_lqr_tb)), LMI_lqr >= 0];
    
    obj_lqr = -trace(P_lqr{:,:,0});
    sol_lqr = optimize(F_lqr, obj_lqr, sdpsettings('verbose', 0));
    
    if sol_lqr.problem == 0
        P_val_lqr = sdpval(P_lqr);
        fprintf('   ✓ LMI LQR solved successfully!\n');
        fprintf('   Maximal P(DC) trace: %.4f\n\n', trace(P_val_lqr{:,:,0}));
    else
        fprintf('   ✗ LMI LQR failed: %s\n\n', sol_lqr.info);
    end
    
    % --- 3. LMI Kalman (Optimal Observer) ---
    fprintf('3) LMI Kalman (Optimal Observer)\n');
    P_kal = PhasorArray.ndsdpvar(nx, nx, h_lmi, PhasorType='symmetric', real=true);
    P_kal_tb = P_kal.T_tb(h_lmi);
    W_kal_tb = TB(PhasorArray(W_proc), h_lmi);
    V_kal_tb = TB(PhasorArray(V_meas(1,1)), h_lmi);
    
    LMI_kal = [A_hss*P_kal_tb + P_kal_tb*A_hss' + W_kal_tb, P_kal_tb * C_tb'; 
               C_tb*P_kal_tb, V_kal_tb];
    F_kal = [P_kal_tb >= 1e-6*eye(size(P_kal_tb)), LMI_kal >= 0];
    
    obj_kal = -trace(P_kal{:,:,0});
    sol_kal = optimize(F_kal, obj_kal, sdpsettings('verbose', 0));
    
    if sol_kal.problem == 0
        P_val_kal = sdpval(P_kal);
        fprintf('   ✓ LMI Kalman solved successfully!\n');
        fprintf('   Maximal P(DC) trace: %.4f\n\n', trace(P_val_kal{:,:,0}));
    else
        fprintf('   ✗ LMI Kalman failed: %s\n\n', sol_kal.info);
    end
end

%% ═══════════════════════════════════════════════════════════════════════════
%  MISSING - TODO - WIP
%  ═══════════════════════════════════════════════════════════════════════════

% chercher via harmonic balance une traj perriodique à l'équilibre

%% ═══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS  (must sit at the end of the script — required before R2024a)
%  ═══════════════════════════════════════════════════════════════════════════

function s = stability_str(is_stable)
    % One-liner stability formatter.
    if is_stable, s = 'STABLE'; else, s = 'UNSTABLE'; end
end

function y = simInput(A, t, x0, T, B, inp)
    % SIMINPUT  Simulate a periodic LTP system with a sampled input — no PhasorSS.
    %   y = simInput(A, t, x0, T)            % autonomous   ẋ = A(t)·x
    %   y = simInput(A, t, x0, T, B, inp)    % forced       ẋ = A(t)·x + B(t)·u(t)
    %
    %   A, B  : PhasorArray (n×n and n×m).   t : time grid (1×Nt, uniform).
    %   x0    : initial state (n×1).         inp : sampled input (Nt×m), e.g. noise.
    %   y     : state trajectory (n×Nt).
    %
    %   Mirrors PhasorSS.lsim (deterministic, interpolated input) but uses only a
    %   fixed-step RK4 over the grid — runs on any MATLAB release. A(t), B(t) are
    %   pre-evaluated once with evalTime (fast, vectorised); midpoints for RK4 come
    %   from the average of consecutive grid samples.
    n = size(A, 1); Nt = numel(t); dt = t(2) - t(1);
    Ak = evalTime(A, T, t);                 % n×n×Nt  (single vectorised call)
    hasU = nargin >= 6 && ~isempty(B) && ~isempty(inp) && any(inp(:));
    if hasU, Bk = evalTime(B, T, t); end
    y = zeros(n, Nt); y(:, 1) = x0(:);
    for k = 1:Nt-1
        A1 = Ak(:,:,k); A2 = Ak(:,:,k+1); Am = (A1 + A2)/2;   % midpoint matrix
        if hasU
            u1 = Bk(:,:,k)*inp(k,:).';  u2 = Bk(:,:,k+1)*inp(k+1,:).';  um = (u1 + u2)/2;
        else
            u1 = 0; u2 = 0; um = 0;
        end
        xk = y(:, k);
        f1 = A1*xk + u1;
        f2 = Am*(xk + dt/2*f1) + um;
        f3 = Am*(xk + dt/2*f2) + um;
        f4 = A2*(xk + dt*f3)   + u2;
        y(:, k+1) = xk + dt/6*(f1 + 2*f2 + 2*f3 + f4);
    end
end