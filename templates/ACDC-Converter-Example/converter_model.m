% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END

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
