% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END

%% 4a. Close the loop: Ã_cl = Ã − B̃K
A_cl = A_aug - B_aug * K;     % 4×4 PhasorArray (closed-loop)

%% 4b. Expand with angle state: dθ/dt = ω_m · pp
nx_ext = nx_aug + 1;

A_ext = [A_cl,                               PhasorArray(zeros(nx_aug, 1));
         PhasorArray([0, 0, param.pp, 0]),   PhasorArray(0)];

% B_ext: voltage inputs (from B_aug) + speed reference (enters z-dot as -ω*)
B_ext = [[B_aug; PhasorArray(zeros(1, nu))] [zeros(3,1) ; -1  ; 0]];

% Output: all states
C_ext = PhasorArray(eye(nx_ext));
D_ext = PhasorArray(zeros(nx_ext, nu));

plot(A_ext)   % visualise extended closed-loop dynamics

%% 4c. Build PhasorSS and set LPV mode
sys_cl = PhasorSS(A_ext, B_ext, C_ext, D_ext, [], ...
    'StateName', {'i_α', 'i_β', 'ω_m', 'z', 'θ'}, ...
    'OutputName', {'i_α', 'i_β', 'ω_m', 'z', 'θ'}, ...
    'InputName', {'v_α', 'v_β', 'ω_r'});

% Set LPV: θ = x(end), the last state IS the scheduling parameter
sys_cl = sys_cl.setLPV(@(t, x, u) x(end));

plot(sys_cl, false)   % PhasorSS - visualisation of A B C D (theta) matrices

%% 4d. Simulate — speed ramp from 0 to 800 RPM over 2 s
omega_0  = 0;                              % start from standstill
theta_0  = 0;
x0 = [0; 0; omega_0/param.pp; 0; theta_0];

t_sim = linspace(0, 3, 5000);
u_sim = zeros(length(t_sim), nu+1);

% Speed reference: ramp from 0 to 800 RPM over [0, 2] s, then constant
u_sim(:,end) = 800/60*2*pi .* ((t_sim)/2.*(t_sim<=2) + (t_sim>2)) ...
               + omega_0 / param.pp;

[y, t_out, x_out] = lsim(sys_cl, t_sim, u_sim, x0);

%% 4e. Plot results
subplot(3,1,1);
plot(t_out, x_out(:,1), t_out, x_out(:,2));
ylabel('Current [A]'); legend('i_\alpha', 'i_\beta');
title('Closed-loop SPMSM — LPV simulation');

subplot(3,1,2);
plot(t_out, x_out(:,3) * 60/(2*pi));
hold on;
plot(t_out, u_sim(:,end) * 60/(2*pi));
ylabel('Speed [RPM]'); legend('\omega_m', '\omega_r');

subplot(3,1,3);
plot(t_out, mod(x_out(:,5), 2*pi));
ylabel('\theta [rad]'); xlabel('Time [s]');
legend('\theta_{elec}');
