% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END

%% 4a. Build the error dynamics matrix  Â = A + A(D^e)
A_err = A + A_de;          % 4×4 PhasorArray  (constant + periodic)

%% 4b. Lyapunov weight Q
alpha = 1e-4;              % relative weight on v_dc vs currents
Q = PhasorArray(blkdiag(eye(3), alpha));

%% 4c. Solve harmonic Lyapunov equation
h_solve = 10;              % truncation order for solving

P = Lyap_Harmonique(A_err, Q, h_solve, T);

% Inspect harmonic content of P
stem(P);   sgtitle('Lyapunov matrix P — harmonic content');
plot(P);   sgtitle('Lyapunov matrix P(t)');

%% 4d. Tune gain H1 (scalar)
%   H1 = (1/50) / sigma_max( G(X^e) ^H P )
%   where G(X^e) is the bilinear coupling evaluated at the setpoint

G_Xe = [C33/param.Lg * Vdc_e_pa * PhasorArray(eye(3));
        I_abc_e' / param.C];      % 4×3 PhasorArray

% Compute Toeplitz and extract max singular value
GP_tb = T_tb(G_Xe' * P, h_solve);
sigma_max_GP = max(svd(full(GP_tb)));
H1 = (1/50) / sigma_max_GP;

fprintf('Stabilizing gain H1 = %.4f\n', H1);
