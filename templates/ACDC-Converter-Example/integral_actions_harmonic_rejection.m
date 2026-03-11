% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END

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

M = Sylv_harmonique(O_pa, A_err, L_pa * Cout, h_solve, T);

% Inspect M
stem(M);   sgtitle('Forwarding matrix M — harmonic content (O=0)');
plot(M);   sgtitle('M(t) — time domain');

%% 5d. Tune H2 (forwarding aggressiveness)
%   H2 = β * blkdiag(1, 0.1, I_2, I_2)
GM_tb = T_tb(G_Xe' * M' * M, h_solve);
sigma_max_GM = max(svd(full(GM_tb)));
beta = (1/H1) * (1/50) / sigma_max_GM;

H2_diag = blkdiag(1, 0.1, eye(2), eye(2));
H2 = PhasorArray(beta * H2_diag);

fprintf('Forwarding gain β = %.4f\n', beta);
