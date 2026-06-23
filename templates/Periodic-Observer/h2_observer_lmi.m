% TEMPLATE-HEADER-BEGIN
% This script is the source of truth for the wiki template.
% Edit this file and run: python scripts/inject_templates.py
% The tag <template/...> in the wiki will be replaced by this content.
% TEMPLATE-HEADER-END HERE. Edit the wiki or regenerate.

%% ── H₂ Observer (Periodic Kalman Filter) ────────────────────────────────
%  dx/dt = A(t)x + w,    y = C(t)x + v,    e = x - x̂
%  w ~ N(0, W(t)),       v ~ N(0, V(t))

%% 0. Problem data
nx = 3;  ny = 2;  ha = 3;
T = 0.5;

rng(42);
A = PhasorArray(-5*eye(nx)) + PhasorArray.random(nx, nx, ha)*0.1; % strictly stable open-loop for illustration
C = PhasorArray.random(ny, nx, ha);

% Process and measurement noise covariances
W = PhasorArray(0.1 * eye(nx));    
V = PhasorArray(eye(ny));          

hP = 6;  h = 6;

%% 1. Decision variables 
% Find P(t) minimizing trace(P) subject to Riccati inequality
P  = PhasorArray.ndsdpvar(nx, nx, hP, 'PhasorType', 'symmetric', 'real', true);

PT = P.T_tb(h);
N  = N_tb(nx, h, T);

%% 2. Periodic Riccati LMI for observer (subsolution from below)
% The Kalman covariance Sigma satisfies dSigma/dt = F(Sigma) where
%   F(P) = A P + P A' + W - P C' V^-1 C P.
%
% Any P satisfying the SUBSOLUTION condition dP/dt <= F(P) is a lower bound
% on the covariance (P(t) <= Sigma(t) by comparison).  The MAXIMUM P in the
% feasible set equals Sigma itself.
%
% Schur form of dP/dt <= F(P):
%   [A P + P A' + W - dP/dt,  P C']  >= 0
%   [         C P           ,   V  ]
%
% Objective: MAXIMIZE trace(P_0), i.e., minimize -trace(P_0).

AP   = T_tb((A * P + P * A'), h);
Pdot = N * PT - PT * N;              % T(dP/dt) = N*T(P) - T(P)*N
PC   = T_tb((P * C'), h);
WT   = W.T_tb(h);
VT   = V.T_tb(h);

LMI_tb = [AP - Pdot + WT,  PC;
          PC',              VT];
LMI_tb = 0.5 * (LMI_tb + LMI_tb'); % Enforce exact symmetry numerically

F = [ PT >= 0 ];
F = [ F, LMI_tb >= 0 ];

%% 4. Solve
% Maximize Trace(P_0) — the maximum subsolution is P = Sigma (Kalman covariance).
obj = -trace(P{:,:,0});
sol = optimize(F, obj, sdpsettings('solver', 'mosek', 'verbose', 1));

%% 5. Extract observer gain: L = P C' V⁻¹
if sol.problem == 0
    PP = PhasorArray(value(P.value));
    
    % Gain calculation (V should be invertible)
    L  = (PP * C.') / V;
    L = L.reduce('reduceMethod', 'relative', 'reduceThreshold', 1e-4);

    fprintf('Observer variance (Trace P): %.6f\n', value(obj));
    figure('Name','h2_observer_lmi — Observer gain'); stem(L); sgtitle('Observer gain L(\theta)');
    figure('Name','h2_observer_lmi — Observer eigenvalues'); plot(HmqNEig(A - L*C, h, T), 'o'); title('H₂-optimal observer eig');
    figure('Name','h2_observer_lmi — Error dynamics'); lsim(A - L*C, 4, [], T); title('Observer error dynamics');
else
    warning('H₂ observer LMI infeasible: %s', sol.info);
end
