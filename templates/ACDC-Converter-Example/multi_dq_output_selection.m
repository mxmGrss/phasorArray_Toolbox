% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END

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
figure; spy(C_TB); title('Sparsity of \mathcal{C} = \mathcal{T}(C)');
