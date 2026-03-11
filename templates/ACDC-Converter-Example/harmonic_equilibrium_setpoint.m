% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END

%% 3a. Compute dq equilibrium from steady-state equations
%   In steady-state dq:  0 = -rg*idq - Lg*R*ω*idq - ddq*Vdc - edq
%                         0 = (3/2)*ddq'*idq - iload

Vdc_e = param.Vdc_ref;
iload_n = Vdc_e / param.Rload;    % nominal load current [A]

% Grid voltage in dq: e_d = -sqrt(3)*Erms, e_q = 0
%   (depends on convention; here e_a = -sqrt(2)*Erms*cos(ωt))
e_d = -sqrt(3) * param.Erms;
e_q = 0;

% Impose i_q = 0 → power factor = 1
i_q = 0;

% From DC bus: (3/2)*(d_d*i_d + d_q*i_q) = iload_n
% From current d-axis: 0 = -rg*i_d + Lg*ω*i_q - d_d*Vdc - e_d
% From current q-axis: 0 = -rg*i_q - Lg*ω*i_d - d_q*Vdc - e_q

% Solve for i_d from power balance: P_ac = P_dc
%   (3/2)*e_d*i_d = Vdc * iload_n  (lossless approx, then refine)
i_d = (2/3) * Vdc_e * iload_n / e_d;

% Duty cycles from current equations
d_d = (-param.rg * i_d + param.Lg * param.omega * i_q - e_d) / Vdc_e;
d_q = (-param.rg * i_q - param.Lg * param.omega * i_d - e_q) / Vdc_e;

fprintf('Equilibrium: i_d=%.3f A, i_q=%.3f A, d_d=%.4f, d_q=%.4f\n', ...
    i_d, i_q, d_d, d_q);

%% 3b. Build harmonic setpoint (PhasorArrays)
%  Relation:  I_{a,1} = (1/sqrt(6))*(i_d - j*i_q)
%             D_{a,1} = (1/sqrt(6))*(d_d - j*d_q)

% Phase-a current setpoint
Ia_e = ScalarPhasorArray([0, (i_d - 1i*i_q)/sqrt(6)], 'isreal', true);
I_abc_e = Ia_e.PhaseShift([0; -2*pi/3; 2*pi/3]);

% DC voltage setpoint (constant)
Vdc_e_pa = PhasorArray(Vdc_e);   % scalar, only 0th phasor

% Full state setpoint
Xe = [I_abc_e; Vdc_e_pa];       % 4×1 PhasorArray

% Phase-a duty cycle setpoint
Da_e = ScalarPhasorArray([0.5, (d_d - 1i*d_q)/sqrt(6)], 'isreal', true);
D_abc_e = Da_e.PhaseShift([0; -2*pi/3; 2*pi/3]);
De = D_abc_e;                    % 3×1 PhasorArray

plot(Xe);    sgtitle('Nominal periodic setpoint x^e(t)');
plot(De);    sgtitle('Nominal duty cycle d^e(t)');

%% 3c.  Verify equilibrium: 0 = (A + A(d^e) - N) x^e + B v^e
% Build A(d^e): bilinear coupling evaluated at the setpoint
A_de = [PhasorArray(zeros(3,3)),  C33/param.Lg * De * PhasorArray(1);
        De' / param.C,            PhasorArray(0)];

% Disturbance (grid voltage + load current)
Ve = [E_abc; PhasorArray(iload_n)];

% Check: should be ≈ 0
residual = A * Xe + A_de * Xe + B * Ve;
fprintf('Equilibrium residual (should be ~0): %.2e\n', ...
    max(abs(residual.pvalue()), [], 'all'));
