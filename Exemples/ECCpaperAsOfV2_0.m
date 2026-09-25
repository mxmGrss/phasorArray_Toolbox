% ECC example for toolbox v2.0.0 (YALMIP + SDP solver).
% Published listings: ECC_ex.m. Differences: docs/ECC-reproduction.md.

%% Model
A_0 = [1.5 1; 1 -0.5];
A_1 = [-4/pi^2 1/2; 0 1/(pi*1i)];
A_2 = [0 0; 1i/2 0];
A_3 = [-4/(3*pi)^2 0; 0 1/(pi*3*1i)];
A = PhasorArray(cat(3,conj(A_3), conj(A_2), conj(A_1), A_0, A_1, A_2, A_3));
A = PhasorArray(cat(3, A_0, A_1, A_2, A_3),'isreal',true)
T = 1; %period
N = 6; %2^N points used for FFT 
At = @(t) [1+sawtooth(2*pi*t/T,0.5)+0.5, 1+cos(2*pi*t/T); 1-sin(2*2*pi*t/T), -0.5 + square(2*pi*t/T)/2];
A = PhasorArray.funcToPhasorArray(At,T,N)
figure
bar(A,scale="linear",uniformYLim=false)
sgtitle('Harmonic coefficients of A(t)')
figure
plot(A)
sgtitle('time evolution of A(t) over one period')

%% Harmonic reduction
A_neglect = neglect(A{2,2}, 2.e-2, 'reduceMethod', 'absolute');
A_trunc = trunc(A{2,2}, 5);
figure
bar(A{2,2}, A_neglect, A_trunc, layout="grouped", scale="log", ...
    labels=["Original", "Neglected", "Truncated"])
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

%% Algebra (A_alg is uniformly positive definite)
A_alg = PhasorArray([2, 0.2; 0.2, 1.5]) ...
    + 0.25*PhasorArray.cos()*PhasorArray.eye(2);
rng(0); % reproducible right-hand side; no seed can repair a singular A
B_alg = PhasorArray.random(2, 2, 2);
C_alg = A_alg + B_alg; % Addition
D_alg = A_alg * B_alg; % Multiplication
Ainv_alg = inv(A_alg); % Sampled inversion, checked in harmonic norm below
[E_alg, divisionInfo] = mlHmcDivide(A_alg,B_alg,autoUpdateh=true, ...
    thresholdResidual=1e-8);
assert(divisionInfo.status==0 && divisionInfo.resrelnorm<=1e-8, ...
    'ECC:Division','The regular algebra example did not converge.');
inverseResidual = norm(reshape(pvalue(A_alg*Ainv_alg-PhasorArray.eye(2)),[],1))/sqrt(2);
assert(inverseResidual<=1e-8,'ECC:Inverse','Inverse residual exceeds 1e-8.');
At_alg = A_alg.'; % Transpose
Ah_alg = A_alg'; % Conjugate transpose

%% Harmonic operators
h = 8;
A_tb = A.T_tb(h); % Toeplitz-Block form (full matrix)
figure
barsurf(abs(A_tb),1e-3);
title('T(A)_8');
A_four = F_tb(A,h) % Fourier form of PhasorArray A
T = 1;
h = 10; % Truncation order
Ntb = N_tb(size(A,1), h, T);
lambda = A.HmqNEig(h,T,'fundamental')

%% Lyapunov
Q = PhasorArray.eye(2); % Positive definite Q(t)
[P, lyapunovInfo] = lyap(A,Q,'T',T);
assert(lyapunovInfo.status==0,'ECC:Lyapunov','Lyapunov solve did not converge.');
t=0:0.1:T;
Pt = evalTime(P,T,t);

%% LQR
B = [1;PhasorArray.sin];
Q = PhasorArray.eye(2)*10;
R = PhasorArray.eye(1);
K0 = PhasorArray([10,10]);
htrunc = 6; %initial truncation order for Riccati solution
[K_final, S_final, info] = hare(A, B, Q, R, "K0", K0, "T", T, "autoUpdateh", true, "maxIter", 50, "thresholdResidual", 1e-6, "h", htrunc, "maxh", 500, "verbose", 2)
assert(info.status==0 && info.resRicnorm<=1e-6,'ECC:Riccati','Riccati acceptance failed.');
eig_closedloop=HmqNEig(A-B*K_final,20,T,'fundamental')
paperPoles=sort([-3.4466;-2.3234]);
poleError=max(abs(sort(real(eig_closedloop(:)))-paperPoles));
assert(all(real(eig_closedloop)<0) && poleError<1e-3, ...
    'ECC:ClosedLoop','Closed-loop exponents differ from the paper by >=1e-3.');

%% LMI formulation
hlmi = 20; %truncation order for LMI
hP = 10; %truncation order for P in LMI
ht = 10; %truncation order for A in LMI
P = PhasorArray.ndsdpvar(2, 2, hP);
LMIvar11 = trunc(A,ht).'*P + P*trunc(A,ht)+ P.d(T) + Q;
LMIvar12 = P*B;
LMIvar22 = R;
P_tb = P.T_tb(hlmi); 
LMI_tb = [T_tb(LMIvar11, hlmi) , T_tb(LMIvar12, hlmi);
 T_tb(LMIvar12.', hlmi), T_tb(LMIvar22, hlmi)];
Constraints = [P_tb >= 0, LMI_tb >= 0];
Objective = -trace(P{:,:,0});%use P{i,j,k} to access kth phasor of P_{i,j}
sol = optimize(Constraints, Objective);
assert(sol.problem==0,'ECC:LMI','LMI solve failed: %s',sol.info);
Psol = sdpval(P); % Extract solution as PhasorArray
K=inv(R)*B'*Psol % Optimal feedback
figure
plot(Psol, T, 0:T/100:T);
hold on;
plot(S_final, T, 0:T/100:T)
figure
plot(K, T, 0:T/100:T)

%% Closed-loop simulation
C = PhasorArray.eye(2);
D = PhasorArray.zeros(2, 2);
sys = PhasorSS(A, B, C, D, T);
sys_cl = feedback(sys, K_final);
figure
plot(sys_cl) %plot of matrix [A-BK,B;C,-DK]
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
