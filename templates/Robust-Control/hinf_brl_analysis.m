% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END

%% ── H-infinity Analysis (Bounded Real Lemma) ───────────────────────────
% Computes the H∞ norm (L2-gain) from disturbance w to output z
% for a periodic system:  dx/dt = A(t)x + B_w(t)w
%                             z = C_z(t)x + D_z(t)w

%% 0. Problem data
nx = 3;  nw = 2;  nz = 2;  ha = 3;
T = 0.1;

rng(42);
A  = PhasorArray(-5*eye(nx)) + PhasorArray.random(nx, nx, ha)*0.1; % strictly stable open-loop
Bw = PhasorArray.random(nx, nw, ha);
Cz = PhasorArray.random(nz, nx, ha);
Dz = PhasorArray.random(nz, nw, ha);

hP = 6;  h = 6;

%% 1. Decision variables
P   = PhasorArray.ndsdpvar(nx, nx, hP);
gam = sdpvar(1, 1);

PT = P.T_tb(h);
N  = N_tb(nx, h, T);

%% 2. Build BRL LMI
% [ A'P + PA - dP/dt + Cz'Cz,  P Bw + Cz'Dz ]
% [          *              ,   -gam^2 I + Dz'Dz ]  < 0
%
% By Schur complement on gam^2 (since gam > 0), a linear form in gam is:
% [ A'P + PA - dP/dt,  P Bw,     Cz' ]
% [      *          , -gam I,    Dz' ] < 0
% [      *          ,    *  , -gam I ]

AP   = T_tb(A' * P + P * A,h);
Pdot = N * PT - PT * N;
PBw  = T_tb(P * Bw,h);

F11 = AP - Pdot;
F12 = PBw;
F13 = Cz.T_tb(h)';
F23 = Dz.T_tb(h)';

% We enforce exact symmetry internally
LMI = [F11, F12, F13; F12', -gam * eye((2*h+1)*nw), F23; F13', F23', -gam * eye((2*h+1)*nz)];
LMI = 0.5 * (LMI + LMI');

F = [PT >= 0]; 
F = [F, LMI <= -1e-6 * eye(size(LMI))];

%% 4. Solve
sol = optimize(F, gam, sdpsettings('solver', 'mosek', 'verbose', 1));

%% 5. Extract result
if sol.problem == 0
    gamma_opt = value(gam);
    fprintf('Optimal H-infinity norm (γ): %.4f\n', gamma_opt);
    
    P_val = PhasorArray(value(P.value));
    figure('Name','hinf_brl — Lyapunov P(t)'); plot(P_val); sgtitle('Lyapunov Matrix P(t)');
    figure('Name','hinf_brl — Open-loop eigenvalues'); plot(HmqNEig(A, h, T), '*'); title('Open-loop eigenvalues');
else
    warning('LMI infeasible: %s', sol.info);
end
