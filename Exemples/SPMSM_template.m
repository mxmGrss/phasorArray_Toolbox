%% SPMSM Template — Modelling, Control & LPV Simulation
%  End-to-end: αβ model → integral action → polytopic LQR → LPV sim
%
%  Requires: PhasorArray Toolbox, YALMIP + MOSEK, Control System Toolbox

clear; clc; close all;

% Output directory for wiki figures
assetDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'docs', 'wiki_assets');
if ~isfolder(assetDir), mkdir(assetDir); end

%% ── 1. SPMSM αβ Model ──────────────────────────────────────────────────

%% 1a. Motor parameters
param.r     = 0.5;         % stator resistance [Ω]
param.L     = 1.5e-3;      % stator inductance [H]  (constant for SPMSM)
param.J     = 0.033;        % inertia [kg·m²]
param.Bf    = 0.01;         % viscous friction [N·m·s/rad]
param.pp    = 4;            % pole pairs
param.psi_f = 0.1;          % PM flux linkage [Wb] (rms, Concordia convention)

%% 1b. Multi-harmonic back-EMF  e(θ) = sin(θ) + a₅ sin(5θ) + a₇ sin(7θ)
%   Fourier coefficients: harmonic k → coefficient 1/(2j) for sin(kθ)
a5 = 0.05;   a7 = 0.02;    % amplitude of 5th / 7th harmonics (p.u.)
bemf_phasors = zeros(1, 8);
bemf_phasors(2) = 1/(2i);          % fundamental
bemf_phasors(6) = a5/(2i);         % 5th harmonic
bemf_phasors(8) = a7/(2i);         % 7th harmonic
bemf_a = -ScalarPhasorArray(bemf_phasors, 'isreal', true);

%% 1c. Three-phase back-EMF via phase shift
bemf_abc = bemf_a.PhaseShift([0; -2*pi/3; 2*pi/3]);   % 3×1 PhasorArray
bemf_abc = bemf_abc * param.psi_f * param.pp;           % scale to physical units

fig_bemf_time = figure;
plot(bemf_abc)
sgtitle('Three-phase back-EMF e_{abc}(\theta) — time domain');
exportgraphics(fig_bemf_time, fullfile(assetDir, 'spmsm_bemf_abc_time.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

fig_bemf_stem = figure;
stem(bemf_abc)
sgtitle('Three-phase back-EMF e_{abc} — harmonic spectrum');
exportgraphics(fig_bemf_stem, fullfile(assetDir, 'spmsm_bemf_abc_stem.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

%% 1d. abc state-space  (4 states: i_a, i_b, i_c, ω_m)
A_abc = [-param.r/param.L * PhasorArray(eye(3)),  -1/param.L * bemf_abc;
          bemf_abc'/param.J,                       -param.Bf/param.J];

B_abc = [1/param.L * PhasorArray(eye(3));
         PhasorArray(zeros(1, 3))];

%% 1e. Concordia transform → αβ frame (drop zero-sequence)
C32 = PhasorArray.Concordia();           % 3×3 power-invariant Concordia
T_full = blkdiag(C32, PhasorArray(1));   % 4×4 block-diagonal

A_ab0w = T_full * A_abc * T_full';
B_ab0w = T_full * B_abc * C32';

% Extract αβω (rows/cols [1,2,4] — drop zero-sequence index 3)
A = A_ab0w{[1 2 4], [1 2 4]};    % 3×3 PhasorArray
B = B_ab0w{[1 2 4], [1 2]};      % 3×2 PhasorArray

fprintf('αβ model: %d states, %d inputs, %d harmonics\n', ...
    size(A,1), size(B,2), (size(A,3)-1)/2);

%% ── 2. Augmented System (integral action on speed) ─────────────────────

nx = size(A,1);   nu = size(B,2);

% Speed selection row
C_omega = PhasorArray([0 0 1]);        % 1×3, picks ω_m

% Augmented system: x_tilde = [i_α; i_β; ω; z]
A_aug = [A,                      PhasorArray(zeros(nx,1));
         C_omega,                PhasorArray(0)];

B_aug = [B;
         PhasorArray(zeros(1,nu))];

nx_aug = size(A_aug, 1);   % = 4
fprintf('Augmented system: %d states (3 plant + 1 integrator)\n', nx_aug);

%% ── 3. Polytopic LQR ───────────────────────────────────────────────────

%% 3a. Frequency-affine decomposition
%  The αβ model A(θ) has the property that at electrical frequency ω:
%    d/dt x = A(θ)x  with θ = ωt
%  In the lifted domain: (A_T − ωN)X = 0,  so the "effective" dynamics
%  is A_T at period T = 2π/ω.
%
%  For the frequency-affine form  Ã₀ + ω Ã₁, we set:
%    A0 = A_aug    (the full periodic matrix, evaluated at T = 2π/ω)
%    A1 = 0        (for SPMSM with constant L, no explicit ω-dependence in A)
%  The frequency dependence enters only through N(ω) = diag(jk·ω) in the lifted domain.

A0 = A_aug;
A1 = PhasorArray(zeros(nx_aug));     % no explicit ω-affine term for SPMSM with constant L

B0 = B_aug;
B1 = PhasorArray(zeros(nx_aug, nu));


fig_A0 = figure;
plot(A0)
sgtitle('Augmented plant A_0(\theta) — time domain');
exportgraphics(fig_A0, fullfile(assetDir, 'spmsm_A0_time.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

%% 3b. Operating range
omega_bornes = [100, 3000] * 2*pi / 60;   % [RPM] → [rad/s electrical]
omega_bornes = omega_bornes * param.pp;     % mechanical → electrical

%% 3c. LQR weights
Q_lqr = PhasorArray(diag([10, 10, 100, 50000]));   % penalise speed error and integrator
R_lqr = PhasorArray(eye(nu));                   % control effort

Dw = PhasorArray(eye(nx_aug));                  % disturbance input matrix

%% 3d. LMI synthesis
hP   = 10;    % Lyapunov harmonic order
hlmi = 10;    % LMI truncation order

Y = PhasorArray.ndsdpvar(nu, nx_aug, hP, 'PhasorType', 'full');
P = PhasorArray.ndsdpvar(nx_aug, nx_aug, hP);
X = PhasorArray.ndsdpvar(nu, nu, hP);
sdpvar gam;

PTB  = P.T_tb(hlmi);
DwTB = T_tb((Dw * Dw'),hlmi);

F = [PTB >= 1e-6 * eye(size(PTB))];

for ii = 1:numel(omega_bornes)
    omega_i = omega_bornes(ii);
    T_i = 2*pi / omega_i;
    Ni  = N_tb(nx_aug, hlmi, T_i);

    Ai = A0 + omega_i * A1;
    Bi = B0 + omega_i * B1;

    AiP  = T_tb((Ai * P),hlmi);
    BiY  = T_tb((Bi * Y),hlmi);
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
fprintf('Controller K has %d harmonics\n', (size(K,3)-1)/2);

fig_K_stem = figure; stem(K); sgtitle('Gain-scheduled K(\theta)');
exportgraphics(fig_K_stem, fullfile(assetDir, 'spmsm_K_stem.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

fig_K_time = figure;
plot(K)
sgtitle('Gain K(\theta) — time domain');
exportgraphics(fig_K_time, fullfile(assetDir, 'spmsm_K_time.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

fig_P = figure;
plot(P)
sgtitle('Lyapunov matrix P(\theta) — time domain');
exportgraphics(fig_P, fullfile(assetDir, 'spmsm_P_time.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');
%% ── 4. Closed-Loop & LPV Simulation ────────────────────────────────────

%% 4a. Close the loop: Ã_cl = Ã − B̃K
A_cl = A_aug - B_aug * K;     % 4×4 PhasorArray (closed-loop)

%% 4b. Expand with angle state: dθ/dt = ω
%  New state vector: x_ext = [i_α; i_β; ω_m; z; θ]
%  The θ row is [0 0 0 0 0] (θ does not affect the electrical dynamics)
%  but dθ/dt = ω_m·p (electrical speed = mechanical speed × pole pairs)
%
%  ẋ_ext = [A_cl, 0; C_theta, 0] x_ext + [B_cl; 0] r
%  where C_theta = [0, 0, param.pp, 0] selects ω_mech and converts to ω_elec

nx_ext = nx_aug + 1;

A_ext = [A_cl,                               PhasorArray(zeros(nx_aug, 1));
         PhasorArray([0, 0, param.pp, 0]),   PhasorArray(0)];

% Input: speed reference for the integrator (enters as -ω* in the z-dot equation)
% For simplicity, use no external input during free simulation
B_ext = [[B_aug; PhasorArray(zeros(1, nu))] [zeros(3,1) ; -1  ; 0]];

% Output: all states
C_ext = PhasorArray(eye(nx_ext));
D_ext = PhasorArray(zeros(nx_ext, nu));


fig_A_ext = figure;
plot(A_ext)
sgtitle('Extended closed-loop A_{ext}(\theta) — time domain');
exportgraphics(fig_A_ext, fullfile(assetDir, 'spmsm_A_ext_time.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

%% 4c. Build PhasorSS and set LPV mode
sys_cl = PhasorSS(A_ext, B_ext, C_ext, D_ext, [], ...
    'StateName', {'i_α', 'i_β', 'ω_m', 'z', 'θ'}, ...
    'OutputName', {'i_α', 'i_β', 'ω_m', 'z', 'θ'}, ...
    'InputName', {'v_α', 'v_β','ω_r'});

% Set LPV: θ = x(end), the last state IS the scheduling parameter
sys_cl = sys_cl.setLPV(@(t, x, u) x(end));

fig_ss = figure;
plot(sys_cl,false)
sgtitle('Closed-loop PhasorSS — [A , B ; C , D](\theta)');
exportgraphics(fig_ss, fullfile(assetDir, 'spmsm_phasorss_bode.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

%% 4d. Simulate
omega_0  = 0 * 2*pi/60 * param.pp;    % initial electrical speed [rad/s]
theta_0  = 0;                              % initial angle
x0 = [0; 0; omega_0/param.pp; 0; theta_0];  % [iα; iβ; ωm; z; θ]

t_sim = linspace(0, 3, 5000);           % 0.5 s simulation
u_sim = zeros(length(t_sim), nu+1);          % no external voltage (closed-loop handles it)
u_sim(:,end) = 800/60*2*pi.*((t_sim)/2.*(t_sim<=2) + (t_sim>2)) + omega_0 / param.pp;          % no external voltage (closed-loop handles it)

[y, t_out, x_out] = lsim(sys_cl, t_sim, u_sim, x0);

%% 4e. Plot results
figure;
subplot(3,1,1);
plot(t_out, x_out(:,1), t_out, x_out(:,2));
ylabel('Current [A]'); legend('i_\alpha', 'i_\beta');
title('Closed-loop SPMSM — LPV simulation');

subplot(3,1,2);
plot(t_out, x_out(:,3) * 60/(2*pi));
hold on
plot(t_out, u_sim(:,end) * 60/(2*pi));

ylabel('Speed [RPM]'); legend('\omega_m','\omega_r');

subplot(3,1,3);
plot(t_out, mod(x_out(:,5), 2*pi));
ylabel('\theta [rad]'); xlabel('Time [s]');
legend('\theta_{elec}');
sgtitle('Closed-loop SPMSM — LPV simulation (0 \rightarrow 800 RPM ramp)');
exportgraphics(gcf, fullfile(assetDir, 'spmsm_simulation.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');
