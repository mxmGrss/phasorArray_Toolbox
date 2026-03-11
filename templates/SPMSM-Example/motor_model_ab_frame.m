% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END

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

plot(bemf_abc)     % time-domain view
stem(bemf_abc)     % harmonic content

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
