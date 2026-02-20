%% generateWikiFigures.m
% Script to generate illustrative figures for the PhasorArray Toolbox Wiki.
% Figures are saved in docs/wiki_assets/ with transparent backgrounds and SVG format.

clear; clc; close all;

% --- Setup Paths & Folder ---
assetDir = fullfile( 'wiki_assets');
if ~exist(assetDir, 'dir')
    mkdir(assetDir);
end

% Set default figure properties for professional look
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultLineLineWidth', 1.5);

fprintf('Starting Wiki Figure Generation...\n');

%% --- Figure A: Reconstruction (Phasor vs Time) ---
fprintf('Generating Figure A: Reconstruction Comparison (SVG)...\n');
T_sig = 2*pi;
t = linspace(0, T_sig, 1000);
y_true = square(t);

h_list = [1, 3, 10, 40, 120];
colors = lines(length(h_list));

figA = figure('Name', 'reconstruction_comparison', 'Position', [100 100 800 450]);
set(gcf, 'Color', 'none');
plot(t, y_true, 'k--', 'DisplayName', 'Original (Square)', 'LineWidth', 1);
hold on;

for i = 1:length(h_list)
    h = h_list(i);
    % build Ak as positive harmonics (index 1 = harmonic 1)
    Ak = zeros(1, 1, h + 1);
    for k = 1:2:h
        Ak(1, 1, k) = -2j / (pi*k);
    end
    pA = PhasorArray(0, Ak, 'isreal', true);
    y_rec = evalTime(pA, T_sig, t);
    plot(t, squeeze(y_rec), 'Color', colors(i,:), ...
        'DisplayName', sprintf('h = %d', h));
end

xlabel('Time [s]'); ylabel('Amplitude');
sgtitle('Time-Domain Reconstruction vs. Harmonic Order');
legend('Location', 'bestoutside'); grid on;
exportgraphics(figA, fullfile(assetDir, 'reconstruction_h_comparison.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');


%% --- Figure B: Spectral Decay (Energy per harmonic) ---
fprintf('Generating Figure B: Spectral Decay — pageEnergy (SVG)...\n');
rng(42); % Fix seed for reproducibility
h_spec = 15;
A_pos = (rand(1, 1, h_spec+1) + 1i*rand(1, 1, h_spec+1)) .* exp(-0.3*(0:h_spec));
pA_spec = PhasorArray(A_pos(1,1,1), A_pos(1,1,2:end), 'isreal', true);

figB = figure('Name', 'spectral_decay', 'Position', [100 100 600 400]);
set(gcf, 'Color', 'none');
pageEnergy(pA_spec, true, 'none', 'log');
sgtitle('Normalised Energy per Harmonic (log scale)');
exportgraphics(figB, fullfile(assetDir, 'spectral_decay.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');


%% --- Figure C1: Toeplitz Structure (Spy) ---
fprintf('Generating Figure C1: Toeplitz Spy (SVG)...\n');
Ar = PhasorArray.random(3, 3, 5,"time_structure","cmplx"); 
h_tb = 5;
TBM = Ar.T_tb(h_tb);

figC1 = figure('Name', 'toeplitz_spy', 'Position', [100 100 500 500]);
set(gcf, 'Color', 'none');
spy(TBM);
sgtitle(sprintf('Toeplitz Block Matrix Structure (h=%d, n=3)', h_tb));
xlabel('Column Index'); ylabel('Row Index');
exportgraphics(figC1, fullfile(assetDir, 'toeplitz_structure_spy.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');


%% --- Figure C2: Toeplitz BarSurf ---
fprintf('Generating Figure C2: Toeplitz BarSurf (SVG)...\n');
figC2 = figure('Name', 'toeplitz_barsurf', 'Position', [100 100 600 500]);
set(gcf, 'Color', 'none');
barsurf(abs(double(TBM)));
colorbar
sgtitle(sprintf('Toeplitz Magnitude Map (h=%d, n=3)', h_tb));
xlim([-1 (2*h_tb+1)*3+1]);
ylim([-1 (2*h_tb+1)*3+1]);
view([-0 90])
exportgraphics(figC2, fullfile(assetDir, 'toeplitz_structure_barsurf.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

%% --- Figure C2bis: BarSurf ---
fprintf('Generating Figure C2bis: PhasorArray BarSurf (SVG)...\n');
figC2bis = figure('Name', 'phasorarray_barsurf', 'Position', [100 100 600 500]);
set(gcf, 'Color', 'none');
barsurf(Ar);
colorbar
view([-10 50])
title('Harmonic content of A_i(t)')
exportgraphics(figC2bis, fullfile(assetDir, 'phasorArray_barsurf.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');


%% --- Figure C3: Time-Domain Evaluation of Ar ---
fprintf('Generating Figure C3: Time Plot of Ar (SVG)...\n');
figC3 = figure('Name', 'time_plot_Ar', 'Position', [100 100 600 400]);
set(gcf, 'Color', 'none');
plot(Ar,2*pi);
sgtitle('Time-Domain Evaluation of A(t) (T=2*pi)');
grid on;
exportgraphics(figC3, fullfile(assetDir, 'time_domain_evaluation_Ar.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

%% --- Figure C3b: Time-Domain Evaluation of Ar ---
fprintf('Generating Figure C3: Time Plot of Ar (SVG)...\n');
figC3 = figure('Name', 'time_plot3D_Ar', 'Position', [100 100 600 400]);
set(gcf, 'Color', 'none');
plot3D(Ar,2*pi);
sgtitle('Time-Domain Evaluation of A(t) (T=2*pi)');
grid on;
exportgraphics(figC3, fullfile(assetDir, 'time_domain3D_evaluation_Ar.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');


%% --- Figure C4: Harmonic Stem Plot of Ar ---
fprintf('Generating Figure C4: Stem Plot of Ar (SVG)...\n');
figC4 = figure('Name', 'stem_plot_Ar', 'Position', [100 100 600 400]);
set(gcf, 'Color', 'none');
stem(Ar);
sgtitle('Stem Plot of A(t) Harmonics');
grid on;
exportgraphics(figC4, fullfile(assetDir, 'harmonic_spectrum_stem_Ar.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

%% --- Figure C4b: Harmonic Stem3 Plot of Ar ---
fprintf('Generating Figure C4b: Stem3 Plot of Ar (SVG)...\n');
figC4b = figure('Name', 'stem3_plot_Ar', 'Position', [100 100 600 400]);
set(gcf, 'Color', 'none');
stem3(Ar);
sgtitle('3D Stem Plot of A(t) Harmonics');
grid on;
exportgraphics(figC4b, fullfile(assetDir, 'harmonic_spectrum_stem3_Ar.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

%% --- Figure C5: TB vs BT Comparison (BarSurf) ---
fprintf('Generating Figure C5: TB vs BT Comparison (SVG)...\n');
figC5 = figure('Name', 'tb_vs_bt', 'Position', [100 100 1000 450]);
set(gcf, 'Color', 'none');
subplot(1,2,1);
barsurf(abs(double(Ar.T_tb(5))));
title('Toeplitz Block (T\_tb)');
view([0 90]); colorbar;
subplot(1,2,2);
barsurf(abs(double(Ar.T_bt(5))));
title('Block Toeplitz (T\_bt)');
view([0 90]); colorbar;
sgtitle('Comparison: TB vs BT Ordering (h=5)');
exportgraphics(figC5, fullfile(assetDir, 'toeplitz_tb_vs_bt.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');
%% --- Figure D: Harmonic Transfer (DC Gain) ---
fprintf('Generating Figure D: Harmonic Transfer DC Gain (SVG)...\n');
t_sys = 2*pi;
A0 = [-0.1 1; -1 -0.1];
A1 = [0 0.5; 0.5 0];
pA_sys = PhasorArray(A0, A1, 'isreal', true);
B = [0; 1];
C = [1 0];
sys = PhasorSS(pA_sys, B, C, 0, t_sys);

figD = figure('Name', 'harmonic_transfer', 'Position', [100 100 700 500]);
set(gcf, 'Color', 'none');
h_dc = 5;
hmqDcGain(sys, h_dc, t_sys);
colorbar
sgtitle(sprintf('DC Gain (Effective Transfer, h=%d)', h_dc));
view([-0 90])
exportgraphics(figD, fullfile(assetDir, 'harmonic_transfer_dcgain.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');


%% --- Figure E: Gibbs Effect (GIF) ---
fprintf('Generating Figure E: Gibbs Effect Animation (GIF)...\n');
gifFilename = fullfile(assetDir, 'gibbs_effect_animation.gif');
figE = figure('Name', 'gibbs_gif', 'Position', [100 100 600 450]);

h_gif_list = unique(floor(logspace(0, 3, 30)));

for i = 1:length(h_gif_list)
    h = h_gif_list(i);
    Ak = zeros(1, 1, h + 1);
    for k = 1:2:h
        Ak(1, 1, k) = -2j / (pi*k);
    end
    pA = PhasorArray(0, Ak, 'isreal', true);
    y_rec = evalTime(pA, T_sig, t);
    
    plot(t, y_true, 'k--', 'LineWidth', 1); hold on;
    plot(t, squeeze(y_rec), 'r', 'LineWidth', 2);
    sgtitle(sprintf('Gibbs Effect Reconstruction (h = %d)', h));
    xlabel('Time [s]'); ylabel('Amplitude');
    ylim([-1.5 1.5]); grid on;
    hold off;
    
    frame = getframe(figE);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    if i == 1
        imwrite(imind, cm, gifFilename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

%% =====================================================================
%% NEW FIGURES FOR WIKI PAGES THAT WERE MISSING ILLUSTRATIONS
%% =====================================================================

%% --- Figure F: Harmonic Growth (Arithmetic-Operations page) ---
fprintf('Generating Figure F: Harmonic Growth after Multiplication (SVG)...\n');
rng(42);
A_arith = PhasorArray.random(2, 2, 3, "time_structure", "cmplx");
B_arith = PhasorArray.random(2, 2, 5, "time_structure", "cmplx");
C_arith = A_arith * B_arith; % h_C = 3 + 5 = 8

figF = figure('Name', 'harmonic_growth', 'Position', [100 100 900 350]);
set(gcf, 'Color', 'none');
subplot(1,3,1);
pageEnergy(A_arith, true, 'none', 'linear');
title(sprintf('A (h=%d)', A_arith.h)); xlabel('Harmonic'); ylabel('Norm. Energy');
subplot(1,3,2);
pageEnergy(B_arith, true, 'none', 'linear');
title(sprintf('B (h=%d)', B_arith.h)); xlabel('Harmonic'); ylabel('Norm. Energy');
subplot(1,3,3);
pageEnergy(C_arith, true, 'none', 'linear');
title(sprintf('C = A*B (h=%d)', C_arith.h)); xlabel('Harmonic'); ylabel('Norm. Energy');
sgtitle('Harmonic Growth: h_C = h_A + h_B after Cauchy Product');
exportgraphics(figF, fullfile(assetDir, 'harmonic_growth_multiplication.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

%% --- Figure G: Getting-Started illustration ---
fprintf('Generating Figure G: Getting Started Tutorial Output (SVG)...\n');
T_gs = 0.02; w_gs = 2*pi/T_gs;
x_gs = 1 + 0.5 * PhasorArray.cos(0, 1);
A0_gs = [0, w_gs; -w_gs, 0];
A_vary_gs = PhasorArray.eye(2) * PhasorArray.cos(0, 1);
A_gs = PhasorArray(A0_gs) + A_vary_gs;

figG = figure('Name', 'getting_started', 'Position', [100 100 900 400]);
set(gcf, 'Color', 'none');
subplot(1,2,1);
plot(x_gs, T_gs);
title('Periodic Signal x(t) = 1 + 0.5 cos(\omega t)');
xlabel('Time [s]'); ylabel('Amplitude'); grid on;
subplot(1,2,2);
plot(A_gs, T_gs);
sgtitle('Getting Started: Periodic Signal & System Matrix');
exportgraphics(figG, fullfile(assetDir, 'getting_started_output.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

%% --- Figure H: Simulation Guide — LTP system ---
fprintf('Generating Figure H: LTP Simulation (SVG)...\n');
T_sim = 0.3
w_sim = 2*pi/T_sim;

% Build a stable 3-state system with T = 0.1
A0_sim = [-2 1 0; 0 -4 1; -1 0 -7];
A1_sim = [0 0.5 3; 0.5 2 -0.2; 1 0 0.3];
A_sim = PhasorArray(A0_sim, A1_sim, 'isreal', true);
B_sim = [0; 0; 1];
C_sim = eye(3);
sys_ltp = PhasorSS(A_sim, B_sim, C_sim, zeros(3,1), T_sim);

x0_sim = [1; 0.5; -0.3];
t_sim = 0:T_sim/200:(10*T_sim);
u_sig = zeros(1, length(t_sim));
u_sig(t_sim >= 1) = 10; % Step at t=0.5

figH = figure('Name', 'ltp_simulation', 'Position',  [100 100 900 700]);
set(gcf, 'Color', 'none');

% Show the parameter trajectory
subplot(3,1,1);
t_plot = linspace(0, 3, 2000);
p_vals = ones(size(t_plot))*2*pi/T_sim;
plot(t_plot, p_vals, 'b', 'LineWidth', 1.5);
title('Parameter trajectory p(t)');
xlabel('Time [s]'); ylabel('p [rad/s]'); grid on;


subplot(3,1,2);
[y_ini, t_ini] = initial(sys_ltp, x0_sim, 3);
plot(t_ini, squeeze(y_ini));
title('LTP Free Response (initial condition)');
xlabel('Time [s]'); ylabel('States'); grid on;
legend('x_1','x_2','x_3');

subplot(3,1,3);
[y_lsim, t_lsim] = lsim(sys_ltp, t_sim, u_sig, x0_sim);
plot(t_lsim, squeeze(y_lsim));
title('LTP Forced Response (step at t = 1 s)');
xlabel('Time [s]'); ylabel('States'); grid on;
legend('x_1','x_2','x_3');

sgtitle(sprintf('LTP System Simulation (T = %.2f s)', T_sim));
exportgraphics(figH, fullfile(assetDir, 'simulation_ltp.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');


%% --- Figure H2: Simulation Guide — LPV system ---
fprintf('Generating Figure H2: LPV Simulation (SVG)...\n');
% Same system as LTP but with time-varying parameter p(t,x,u)
% p = 2*pi/0.1 * (1 + sigmoid(t, center=1, width=0.5)) + cos(2*pi/0.2 * t).*(t>1.3)
w_fun = @(t) (2*pi/T_sim) * (1 + 2./(1 + exp(-((t - 1.7)/0.25)))) ...
                  + 3*cos(2*2*pi/T_sim * t) .* (t > 1.3);



u_sig(1,t_sim >= 1) = 10; % Step at t=0.5
u_sig(2,:) = w_fun(t_sim); % Step at t=0.5

p_fun = @(t,x,u) x(end);

sys_lpv = PhasorSS(blkdiag(A_sim,0), blkdiag(B_sim,1), [C_sim,zeros(3,1)]);
sys_lpv = sys_lpv.setLPV(p_fun);



figH2 = figure('Name', 'lpv_simulation', 'Position', [100 100 900 700]);
set(gcf, 'Color', 'none');

% Show the parameter trajectory
subplot(3,1,1);
plot(t_sim, u_sig(2,:), 'b', 'LineWidth', 1.5);
title('Trajectory \omega(t)');
xlabel('Time [s]'); ylabel('p [rad/s]'); grid on;

% Free response

u_sig(1,t_sim >= 1) = 0; % Step at t=0.5
u_sig(2,:) = w_fun(t_sim); % Step at t=0.5
subplot(3,1,2);
[y_lpv_ini, t_lpv_ini] = lsim(sys_lpv, t_sim', u_sig', [x0_sim; 0]);
plot(t_lpv_ini, squeeze(y_lpv_ini));
title('LPV Free Response (initial condition)');
xlabel('Time [s]'); ylabel('States'); grid on;
legend('x_1','x_2','x_3');

% Forced response

u_sig(1,t_sim >= 1) = 10; % Step at t=0.5
u_sig(2,:) = w_fun(t_sim); % Step at t=0.5
subplot(3,1,3);
[y_lpv_f, t_lpv_f] = lsim(sys_lpv, t_sim, u_sig, [x0_sim; 0]);
plot(t_lpv_f, squeeze(y_lpv_f));
title('LPV Forced Response (step at t = 0.5 s)');
xlabel('Time [s]'); ylabel('States'); grid on;
legend('x_1','x_2','x_3');

sgtitle('LPV System Simulation — Varying Parameter p(t)');
exportgraphics(figH2, fullfile(assetDir, 'simulation_lpv.svg'), 'BackgroundColor', 'none', 'ContentType', 'vector');

%% --- Figure I: Solvers — Floquet Exponents Convergence (GIF) ---
fprintf('Generating Figure I: Floquet Convergence Animation (GIF)...\n');
T_solv = 0.2;
gifFloquet = fullfile(assetDir, 'floquet_convergence.gif');
figI = figure('Name', 'floquet_convergence', 'Position', [100 100 600 500]);

h_list_floquet = unique(floor(logspace(0, 2.2, 20)));
for idx = 1:length(h_list_floquet)
    h_cur = h_list_floquet(idx);
    ev = A_sim.HmqNEig(h_cur, 2*pi/T_solv);
    plot(real(ev), imag(ev), '*', 'MarkerSize', 8, 'LineWidth', 1.5);
    grid on;
    title(sprintf(['Floquet Exponents (h = %d, $T = %.1f$ s)' ...
        '  (ev of $\\mathcal{A}-\\mathcal{N}(\\omega)$)'], h_cur, T_solv), ...
        'Interpreter', 'latex', 'FontSize', 13);
    xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
    drawnow;

    frame = getframe(figI);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if idx == 1
        imwrite(imind, cm, gifFloquet, 'gif', 'Loopcount', inf, 'DelayTime', 0.4);
    elseif idx == length(h_list_floquet)
        imwrite(imind, cm, gifFloquet, 'gif', 'WriteMode', 'append', 'DelayTime', 2);
    else
        imwrite(imind, cm, gifFloquet, 'gif', 'WriteMode', 'append', 'DelayTime', 0.3);
    end
end

fprintf('All assets generated successfully in %s\n', assetDir);

