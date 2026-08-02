% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% 1a. Converter parameters
param.rg    = 1.15;          % grid-side resistance [Ω]
param.Lg    = 122e-6;        % grid-side inductance [H]
param.C     = 100e-6;        % DC bus capacitance [F]
param.Rload = 120;           % nominal load resistance [Ω]
param.f     = 50;            % grid frequency [Hz]
param.omega = 2*pi*param.f;  % grid pulsation [rad/s]
param.Erms  = 45;            % AC rms phase voltage [V]
param.Vdc_ref = 150;         % DC bus voltage reference [V]

T = 1/param.f;               % fundamental period [s]

%% 1b. Three-phase grid voltage  e_abc(t) = balanced sinusoidal
%   e_a = -sqrt(2)*Erms*cos(ωt)  →  phasor: E_{a,1} = -Erms*sqrt(2)/2
Ea_phasors = zeros(1, 2);     % harmonics 0 and 1
Ea_phasors(2) = -param.Erms * sqrt(2) / 2;   % 1st harmonic (cosine → real)
Ea = ScalarPhasorArray(Ea_phasors, 'isreal', true);

% Balanced three-phase via phase shift
E_abc = Ea.PhaseShift([0; -2*pi/3; 2*pi/3]);   % 3×1 PhasorArray

plot(E_abc);   title('Three-phase grid voltage e_{abc}(t)');
stem(E_abc);   title('Grid voltage — harmonic spectrum');

%% 1c. System matrices (bilinear form)
C33 = eye(3) - ones(3)/3;     % Laplacian matrix

% State matrix A (constant): dx = Ax + G(x)d + Bv
A = PhasorArray(blkdiag(-param.rg/param.Lg * eye(3), 0));

% Input matrix B (constant):
B = PhasorArray(blkdiag(-eye(3)/param.Lg, -1/param.C));

%% 1d. Bilinear coupling  G(x)*d = A(d)*x
%   A(d) = [0, -(C33*d)/Lg ; d'/C, 0]
%   This will be built symbolically when needed for control design

%% 2a. Multi-frequency dq transforms
%   T_k = PhasorArray.Park(0, k, false)  → 2×3 periodic matrix at order k
%   (amplitude-invariant, 2/3 scaling, zero-sequence removed)

T1 = PhasorArray.Park(0, 1, false);   % fundamental dq transform (2×3)
T2 = PhasorArray.Park(0, 2, false);   % 2nd harmonic rotating frame
T3 = PhasorArray.Park(0, 3, false);   % 3rd harmonic rotating frame
T4 = PhasorArray.Park(0, 4, false);   % 4th harmonic rotating frame
T5 = PhasorArray.Park(0, 5, false);   % 5th harmonic rotating frame

plot(T1);   sgtitle('T_1(\theta) — fundamental dq transform');
plot(T2);   sgtitle('T_2(\theta) — 2nd harmonic rotating frame');
stem(T2);   sgtitle('T_2 — phasor content (only \pm 2 non-zero)');

%% 2b. Define targeted harmonics and build combined output matrix
K_target = [2, 4];   % harmonic orders to reject in i_abc

% q-axis row of T_1 (second row) — selects i_q for power factor control
T1_q = T1{2,:};    % 1×3 PhasorArray

% Stack output matrix: C(t) maps x = [i_abc; v_dc] → y
%   y = [v_dc;  i_q;  i_{dq_2};  i_{dq_4}]
Cout = [PhasorArray(zeros(1,3)),  PhasorArray(1);          % v_dc
        T1_q,                     PhasorArray(0);          % i_q
        T2,                       PhasorArray(zeros(2,1)); % i_{dq_2}
        T4,                       PhasorArray(zeros(2,1))];% i_{dq_4}

nz = size(Cout, 1);   % = 7 integrator states

plot(Cout);   sgtitle('Output matrix C(t) — multi-frequency selection');
stem(Cout);   sgtitle('C(t) — harmonic content');

%% 2c. Harmonic domain output matrix  C_TB = T_tb(Cout, h)
h_out = 10;
C_TB = T_tb(Cout, h_out);   % Toeplitz-block form

% Visualize sparsity — each T_k contributes only 2 non-zero diagonals
figure('Name','acdc — Sparsity C_TB'); spy(C_TB); title('Sparsity of \mathcal{C} = \mathcal{T}(C)');

%% 3a. Compute dq equilibrium from steady-state equations
%   In steady-state dq (Amplitude Invariant Park):
%     0 = -rg*id + Lg*w*iq - dd*Vdc - ed
%     0 = -rg*iq - Lg*w*id - dq*Vdc - eq
%     0 = (3/2)*(dd*id + dq*iq) - iload

Vdc_e = param.Vdc_ref;
iload_n = Vdc_e / param.Rload;    % nominal load current [A]

% Grid voltage in dq (Amplitude Invariant): e_d = -E_m, e_q = 0
%   (depends on convention; e_a = -sqrt(2)*Erms*cos(ωt))
E_m = sqrt(2) * param.Erms;
e_d = -E_m;
e_q = 0;

% Impose i_q = 0 → power factor = 1
i_q = 0;

% Solve for i_d exactly from power balance:
%   rg*id^2 + ed*id + (2/3)*Vdc_e*iload_n = 0
a_c = param.rg;
b_c = e_d;
c_c = (2/3) * Vdc_e * iload_n;
disc = b_c^2 - 4*a_c*c_c;

if disc < 0, error('No physical steady-state solution (load too high).'); end
id_sols = (-b_c + [-1, 1]*sqrt(disc)) / (2*a_c);
[~, idx] = min(abs(id_sols)); % Select normal operation root (lowest current)
i_d = id_sols(idx);

% Duty cycles from current equations
d_d = (-param.rg * i_d + param.Lg * param.omega * i_q - e_d) / Vdc_e;
d_q = (-param.rg * i_q - param.Lg * param.omega * i_d - e_q) / Vdc_e;

fprintf('Equilibrium: i_d=%.3f A, i_q=%.3f A, d_d=%.4f, d_q=%.4f\n', ...
    i_d, i_q, d_d, d_q);

%% 3b. Build harmonic setpoint (PhasorArrays)
%  Relation (Amplitude Invariant):  X_{a,1} = (x_d + j*x_q)/2

% Phase-a current setpoint
Ia_e = ScalarPhasorArray([0, (i_d + 1i*i_q)/2], 'isreal', true);
I_abc_e = Ia_e.PhaseShift([0; -2*pi/3; 2*pi/3]);

% DC voltage setpoint (constant)
Vdc_e_pa = PhasorArray(Vdc_e);   % scalar, only 0th phasor

% Full state setpoint
Xe = [I_abc_e; Vdc_e_pa];       % 4×1 PhasorArray

% Phase-a duty cycle setpoint
Da_e = ScalarPhasorArray([0.5, (d_d + 1i*d_q)/2], 'isreal', true);
D_abc_e = Da_e.PhaseShift([0; -2*pi/3; 2*pi/3]);
De = D_abc_e;                    % 3×1 PhasorArray

plot(Xe);    sgtitle('Nominal periodic setpoint x^e(t)');
plot(De);    sgtitle('Nominal duty cycle d^e(t)');

%% 3c.  Verify equilibrium: 0 = A*x^e - d(x^e)/dt + G(x^e)*d^e + B*v^e
% Build bilinear coupling A(d^e)
% Lg * di/dt = ... - vdc * C33 * d / Lg  (NEGATIVE SIGN)
A_de = [PhasorArray(zeros(3,3)), -C33/param.Lg * De;
        De' / param.C,            PhasorArray(0)];

% Disturbance (grid voltage + load current)
Ve = [E_abc; PhasorArray(iload_n)];

% Check: should be ≈ 0
dXe = d(Xe, T);
residual = A * Xe - dXe + A_de * Xe + B * Ve;
fprintf('Equilibrium residual (should be ~0): %.2e\n', ...
    max(abs(residual.value), [], 'all'));

%% 4a. Build the error dynamics matrix  Â = A + A(D^e)
A_err = A + A_de;          % 4×4 PhasorArray  (constant + periodic)

%% 4b. Lyapunov weight Q
alpha = 1e-4;              % relative weight on v_dc vs currents
Q = PhasorArray(blkdiag(eye(3), alpha));

%% 4c. Solve harmonic Lyapunov equation
h_solve = 10;              % truncation order for solving

P = lyap(A_err, Q, "h", h_solve, "T", T);

% Inspect harmonic content of P

figure('Name','acdc — Lyapunov P harmonics'); stem(P);   sgtitle('Lyapunov matrix P — harmonic content');
figure('Name','acdc — Lyapunov P(t)'); plot(P);   sgtitle('Lyapunov matrix P(t)');

%% 4d. Tune gain H1 (scalar)
%   H1 = (1/50) / sigma_max( G(X^e) ^H P )
%   where G(X^e) is the bilinear coupling evaluated at the setpoint

G_Xe = [-Vdc_e_pa * PhasorArray(C33/param.Lg);
         I_abc_e' / param.C];     % 4×3 PhasorArray

% Compute Toeplitz and extract max singular value
GP_tb = T_tb(G_Xe' * P, h_solve);
sigma_max_GP = max(svd(full(GP_tb)));
H1 = (1/50) / sigma_max_GP;

fprintf('Stabilizing gain H1 = %.4f\n', H1);

%% 5a. Reuse C(t) from §2 (multi-frequency output matrix)
%   Cout already defined in §2:
%     y = [v_dc;  i_q;  i_{dq_2};  i_{dq_4}]   (7×1)
%     Cout = [0, 1;  T1_q, 0;  T2, 0;  T4, 0]   (7×4 PhasorArray)
%
%   No oscillator matrix needed: O = 0.

%% 5b. Integral gains L (block-diagonal, constant)
%   L = blkdiag(ℓ_vdc, ℓ_iq, ℓ_2*I_2, ℓ_4*I_2)
ell_vdc = 0.1;            % v_dc integrator gain
ell_iq  = 0.8165;         % i_q integrator gain (≈ sqrt(2/3))
ell_2   = 0.1143;         % i_{dq_2} integrator gain
ell_4   = 0.1143;         % i_{dq_4} integrator gain

L_full = blkdiag(ell_vdc, ell_iq, ell_2*eye(2), ell_4*eye(2));
L_pa = PhasorArray(L_full);   % 7×7 PhasorArray (constant)

%% 5c. Solve harmonic Sylvester equation  (O = 0)
%   -N M - M (A + A(D^e) - N) + L C = 0
%   ⟹ Sylv_harmonique(O=0, A_err, L*C, h, T)

O_pa = PhasorArray(zeros(nz));   % O = 0: no oscillator!

M = lyap(O_pa, A_err, L_pa * Cout, "h", h_solve, "T", T);

% Inspect M
figure('Name','acdc — Forwarding M harmonics'); stem(M);   sgtitle('Forwarding matrix M — harmonic content (O=0)');
figure('Name','acdc — Forwarding M(t)'); plot(M);   sgtitle('M(t) — time domain');

%% 5d. Tune H2 (forwarding aggressiveness)
%   H2 = β * blkdiag(1, 0.1, I_2, I_2)
GM_tb = T_tb(G_Xe' * (M' * M), h_solve);
sigma_max_GM = max(svd(full(GM_tb)));
beta = (1/H1) * (1/50) / sigma_max_GM;

H2_diag = blkdiag(1, 0.1, eye(2), eye(2));
H2 = PhasorArray(beta * H2_diag);

fprintf('Forwarding gain β = %.4f\n', beta);

%% 6a. Extract periodic gain matrices for time-domain implementation
%   P(t) = P_0 + Σ P_{c,k} cos(kθ) + P_{s,k} sin(kθ)
%   M(t) = M_0 + Σ M_{c,k} cos(kθ) + M_{s,k} sin(kθ)

h_max = 5;  % truncation for implementation (P and M decay rapidly)

% PhasorArray stores harmonics; for implementation extract cos/sin form
P_reduced = reduce(P, h_max);
M_reduced = reduce(M, h_max);

fprintf('Implementation truncation: h_max = %d\n', h_max);

%% 6b. Simulation setup
T_sim = 10*T;                          % simulation time [s]
t_sim = linspace(0, T_sim, 50*4000);    % 4 kHz sampling

% Initial condition: diode rectifier equilibrium (uncontrolled)
% Start from approximate values (can be refined)
x0 = [0; 0; 0; 0];                    % [i_a; i_b; i_c; v_dc] at t=0
z0 = zeros(nz, 1);                    % integrator states (7×1)

%% 6c. ODE integration (time-domain controller)
odefun = @(t, xz) acdc_dynamics(t, xz, param, ...
    P_reduced, M_reduced, PhasorArray(H1), -H2*0.1, L_full, ...
    Xe, De, Cout, C33, E_abc);

[t_out, xz_out] = ode15s(odefun, [0 T_sim], [x0; z0]);

x_out = xz_out(:, 1:4);       % states: [i_a, i_b, i_c, v_dc]
z_out = xz_out(:, 5:4+nz);    % integrator states

%% 6d. Plot results
figure('Name','acdc — Closed-loop simulation');
subplot(3,1,1);
plot(t_out, x_out(:,1), t_out, x_out(:,2), t_out, x_out(:,3));
ylabel('Current [A]'); legend('i_a', 'i_b', 'i_c');
title('Three-phase AC/DC converter — Harmonic forwarding control');

subplot(3,1,2);
plot(t_out, x_out(:,4));
hold on; yline(param.Vdc_ref, '--r');
ylabel('v_{dc} [V]'); legend('v_{dc}', 'v_{dc,ref}');

subplot(3,1,3);
% Compute THD of i_a over sliding window
ylabel('THD(i_a) [%]'); xlabel('Time [s]');

%% 6e. Helper function: closed-loop dynamics
function dxz = acdc_dynamics(t, xz, param, P, M, H1, H2, ...
    L, Xe, De, Cout, C33, E_abc)
    
    nx = 4;  nz = size(L, 1);
    x = xz(1:nx);
    z = xz(nx+1:nx+nz);
    
    theta = param.omega * t;
    
    % Evaluate periodic setpoints at current θ
    xe = real(evalp(Xe, theta));
    de = real(evalp(De, theta));
    xbar = x - xe;
    
    % Evaluate periodic gains at current θ
    Pt = real(evalp(P, theta));
    Mt = real(evalp(M, theta));
    Ct = real(evalp(Cout, theta));
    H1t = real(evalp(H1, theta));
    H2t = real(evalp(H2, theta));
    
    % Grid voltage
    e_abc = real(evalp(E_abc, theta));
    
    % Bilinear coupling G(x)^T
    vdc = x(4);
    i_abc = x(1:3);
    Gx_T = [-C33' * vdc / param.Lg,  i_abc / param.C]';   % 3×4
    
    % Forwarding control law
    w = z - Mt * xbar;
    feedback = Pt * xbar - Mt' * H2t * w;
    d = de - H1t * Gx_T' * feedback;
    
    % Saturate duty cycles to [0, 1]
    d = max(0, min(1, d));
    
    % State dynamics: dx = Ax + A(d)x + Bv
    di_abc = (-param.rg * i_abc - C33 * d * vdc - e_abc) / param.Lg;
    dvdc   = (d' * i_abc - vdc/param.Rload-4*(theta>3*2*pi)) / param.C;
    
    % Integrator dynamics: dz = L * C(t) * xbar   (O = 0)
    dz = L * Ct * xbar;
    
    dxz = [di_abc; dvdc; dz];
end
