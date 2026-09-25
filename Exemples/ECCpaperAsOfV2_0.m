% ECC paper example adapted to PhasorArray Toolbox v2.0.0.
% For the published listings, run ECC_ex.m instead.
% See docs/ECC-reproduction.md for every intentional adaptation.
% Requires YALMIP and an SDP solver for the LMI section.
% Define harmonic coefficients
A_0 = [1.5 1; 1 -0.5];
A_1 = [-4/pi^2 1/2; 0 1/(pi*1i)];
A_2 = [0 0; 1i/2 0];
% Erratum: the triangle in At has a NEGATIVE third Fourier coefficient.
A_3 = [-4/(3*pi)^2 0; 0 1/(pi*3*1i)];

% Create PhasorArray with all harmonics (negative to positive)
A = PhasorArray(cat(3,conj(A_3), conj(A_2), conj(A_1), A_0, A_1, A_2, A_3));
% Equivalently for real signal:
A = PhasorArray(cat(3, A_0, A_1, A_2, A_3),'isreal',true)

T = 1; %period
N = 6; %2^N points used for FFT 
At = @(t) [1+sawtooth(2*pi*t/T,0.5)+0.5, 1+cos(2*pi*t/T); 1-sin(2*2*pi*t/T), -0.5 + square(2*pi*t/T)/2];
A = PhasorArray.funcToPhasorArray(At,T,N)

% Harmonic spectrum visualization (stem plot)
figure
stem(A,'scale','linear','uniformYLim',0)
sgtitle('stemplot of phasor of A(t)')
% Time-domain reconstruction over multiple periods
figure
plot(A)
sgtitle('time evolution of A(t) over one period')

%%
% Neglect small harmonics (threshold-based filtering)
A_neglect = neglect(A{2,2}, 2.e-2, 'reduceMethod', 'absolute');
% Truncate to fixed harmonic order
A_trunc = trunc(A{2,2}, 5);
figure
stem(A{2,2}, 'scale', 'log')
hold on
stem(A_neglect, 'Marker', 'x')
stem(A_trunc, 'Marker', '<')
legend('PhasorArray','Neglected Phasors','TruncatedPhasors')
title('Harmonic content (log scale)')
xlabel('Harmonic'), ylabel('Magnitude')
figure
plot(A{2,2}, T, [0, 2*T])
hold on
plot(A_neglect, T, [0, 2*T])
plot(A_trunc, T, [0, 2*T])
t=0:.01:2*T;
A_eval = arrayfun(@(t) At(t), t, 'UniformOutput', false);
catA_eval = cat(3, A_eval{:});
plot(t, squeeze(catA_eval(2,2,:)), 'k--')
legend('PhasorArray', 'Neglected Phasors', 'TruncatedPhasors', 'Original Signal')
title('Square wave'), xlabel('Time (s)'), ylabel('Amplitude')
%%
% Generate random periodic matrix
rng(0); % reproducible algebra example; the paper does not specify a seed
B = PhasorArray.random(2, 2, 2) % 2x2 with 2 harmonics
% Algebraic operations (computed in harmonic domain)
% Equivalence in the time domain
C = A + B; % Addition: C(t) = A(t) + B(t)
D = A * B; % Multiplication: D(t) = A(t)*B(t)
Ainv = inv(A); % Inversion: Ainv(t) = A(t)^{-1}
% A(t) is singular at t=T/8; its Fourier reconstruction also has singularities.
% This illustrates least squares, NOT a converged regular inverse.
% Use the method underlying A\B to retain its numerical diagnostics.
[E, divisionInfo] = mlHmcDivide(A,B,autoUpdateh=true);
if divisionInfo.status~=0
    warning('PhasorArray:ECC:divisionNotConverged', ...
        'Paper A is singular: division status %d, relative residual %.3g at h=%d.', ...
        divisionInfo.status,divisionInfo.resrelnorm,divisionInfo.h);
end
At = A.'; % Transpose
Ah = A'; % Transpose conjugate
%%
% Create Toeplitz-Block form truncated to h=8 harmonics
h = 8;
A_tb = A.T_tb(h); % Toeplitz-Block form (full matrix)
% Visualize the structure
figure
barsurf(abs(A_tb),1e-3);
title('T(A)_8');
%%
A_four = F_tb(A,h) % Fourier form of PhasorArray A

%%
% For system with state dimension n=2 and period T=1
T = 1;
h = 10; % Truncation order
% Construct N operator in Toeplitz-Block form
Ntb = N_tb(size(A,1), h, T);

%%
lambda = A.HmqNEig(h,T,'fundamental')

%%
Q = PhasorArray.eye(2); % Positive definite Q(t)
[P, lyapunovInfo] = lyap(A,Q,'T',T);
assert(lyapunovInfo.status==0,'ECC:Lyapunov','Lyapunov solve did not converge.');

%%
t=0:0.1:T;
Pt = evalTime(P,T,t);

%%
B = [1;PhasorArray.sin];
Q = PhasorArray.eye(2)*10;
R = PhasorArray.eye(1);
K0 = PhasorArray([10,10]);
htrunc = 6; %initial truncation order for Riccati solution
% v2 hare uses warmStartFraction=1.0; the published listing explicitly uses 0.95.
[K_final, S_final, info] = hare(A, B, Q, R, "K0", K0, "T", T, "autoUpdateh", true, "maxIter", 50, "thresholdResidual", 1e-6, "h", htrunc, "maxh", 500, "verbose", 2)
assert(info.status==0 && info.resRicnorm<=1e-6,'ECC:Riccati','Riccati acceptance failed.');

%%
eig_closedloop=HmqNEig(A-B*K_final,20,T,'fundamental')
paperPoles=sort([-3.4466;-2.3234]);
poleError=max(abs(sort(real(eig_closedloop(:)))-paperPoles));
assert(all(real(eig_closedloop)<0) && poleError<1e-3, ...
    'ECC:ClosedLoop','Closed-loop exponents differ from the paper by >=1e-3.');

%%
% TBLMI LQR Optimization Problem
hlmi = 20; %truncation order for LMI
hP = 10; %truncation order for P in LMI
ht = 10; %truncation order for A in LMI
% Define decision variable (periodic P via harmonics)
P = PhasorArray.ndsdpvar(2, 2, hP);
% Formulate LMI components
LMIvar11 = trunc(A,ht).'*P + P*trunc(A,ht)+ P.d(T) + Q;
LMIvar12 = P*B;
LMIvar22 = R;
% Convert to Toeplitz-Block form
P_tb = P.T_tb(hlmi); 
LMI_tb = [T_tb(LMIvar11, hlmi) , T_tb(LMIvar12, hlmi);
 T_tb(LMIvar12.', hlmi), T_tb(LMIvar22, hlmi)];
% Solve: maximize trace(P_0) subject to constraints
Constraints = [P_tb >= 0, LMI_tb >= 0];
Objective = -trace(P{:,:,0});%use P{i,j,k} to access kth phasor of P_{i,j}
sol = optimize(Constraints, Objective);
assert(sol.problem==0,'ECC:LMI','LMI solve failed: %s',sol.info);

%%
Psol = sdpval(P); % Extract solution as PhasorArray
K=inv(R)*B'*Psol % Optimal feedback
% Compare solutions
figure
plot(Psol, T, 0:T/100:T);
hold on;
plot(S_final, T, 0:T/100:T)

figure
plot(K, T, 0:T/100:T)

%%
% Create LTP state-space system
C = PhasorArray.eye(2);
D = PhasorArray.zeros(2, 2);
sys = PhasorSS(A, B, C, D, T);
% Closed-loop system
sys_cl = feedback(sys, K_final);

figure
plot(sys_cl) %plot of matrix [A-BK,B;C,-DK]

%%
% Simulate initial condition response
x0 = [1; 1];
t_sim = 0:T/100:5*T;

figure
initial(sys_cl, x0, t_sim);
title('Closed-loop response (initial condition)')
figure
u= [1+PhasorArray.cos];
[y,t]=lsim(sys_cl,t_sim,u,x0);
plot(t,y);
title('Closed-loop response')
