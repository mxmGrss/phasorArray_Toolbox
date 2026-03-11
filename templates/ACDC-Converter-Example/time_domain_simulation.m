% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END

%% 6a. Extract periodic gain matrices for time-domain implementation
%   P(t) = P_0 + Σ P_{c,k} cos(kθ) + P_{s,k} sin(kθ)
%   M(t) = M_0 + Σ M_{c,k} cos(kθ) + M_{s,k} sin(kθ)

h_max = 5;  % truncation for implementation (P and M decay rapidly)

% PhasorArray stores harmonics; for implementation extract cos/sin form
P_reduced = reduce(P, h_max);
M_reduced = reduce(M, h_max);

fprintf('Implementation truncation: h_max = %d\n', h_max);

%% 6b. Simulation setup
T_sim = 0.3;                          % simulation time [s]
t_sim = linspace(0, T_sim, 50000);    % 50 kHz sampling

% Initial condition: diode rectifier equilibrium (uncontrolled)
% Start from approximate values (can be refined)
x0 = [0; 0; 0; 0];                    % [i_a; i_b; i_c; v_dc] at t=0
z0 = zeros(nz, 1);                    % integrator states (7×1)

%% 6c. ODE integration (time-domain controller)
odefun = @(t, xz) acdc_dynamics(t, xz, param, ...
    P_reduced, M_reduced, H1, H2, L_full, ...
    Xe, De, Cout, C33, E_abc);

[t_out, xz_out] = ode45(odefun, t_sim, [x0; z0]);

x_out = xz_out(:, 1:4);       % states: [i_a, i_b, i_c, v_dc]
z_out = xz_out(:, 5:4+nz);    % integrator states

%% 6d. Plot results
figure;
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
    xe = evaluate(Xe, theta);
    de = evaluate(De, theta);
    xbar = x - xe;
    
    % Evaluate periodic gains at current θ
    Pt = evaluate(P, theta);
    Mt = evaluate(M, theta);
    Ct = evaluate(Cout, theta);
    
    % Grid voltage
    e_abc = evaluate(E_abc, theta);
    
    % Bilinear coupling G(x)^T
    vdc = x(4);
    i_abc = x(1:3);
    Gx_T = [C33' * vdc / param.Lg,  i_abc / param.C]';   % 3×4
    
    % Forwarding control law
    w = z - Mt * xbar;
    feedback = Pt * xbar - Mt' * H2 * w;
    d = de - H1 * Gx_T' * feedback;
    
    % Saturate duty cycles to [0, 1]
    d = max(0, min(1, d));
    
    % State dynamics: dx = Ax + A(d)x + Bv
    di_abc = (-param.rg * i_abc - C33 * d * vdc - e_abc) / param.Lg;
    dvdc   = (d' * i_abc - vdc/param.Rload) / param.C;
    
    % Integrator dynamics: dz = L * C(t) * xbar   (O = 0)
    dz = L * Ct * xbar;
    
    dxz = [di_abc; dvdc; dz];
end
