%% Quick Test Script for shift2pi
% Run this script to quickly verify that shift2pi works correctly

clear; clc;
fprintf('=== Quick shift2pi Test ===\n');

% Generate test data: 2 seconds of 5 Hz rotation with acceleration
fs = 1000;  % Sample rate
t = (0:1/fs:2)';
omega_base = 2*pi*5;  % 5 Hz base rotation
omega_var = 2*pi*1*sin(2*pi*0.5*t);  % ±1 Hz variation
omega = omega_base + omega_var;
theta = cumtrapz(t, omega);

% Test signal: 50 Hz sinusoid
U = sin(2*pi*50*t);

fprintf('Test data:\n');
fprintf('  Duration: %.1f seconds\n', t(end));
fprintf('  Angular range: %.1f revolutions\n', (theta(end)-theta(1))/(2*pi));
fprintf('  Sample points: %d\n', length(t));

% Test the function
fprintf('\nTesting shift2pi...\n');
tic;
[UHat, istart, i_first_rev, thetaHat, timeHat] = shift2pi(theta, t, U);
elapsed = toc;

fprintf('Results:\n');
fprintf('  ✓ Execution time: %.3f ms\n', elapsed*1000);
fprintf('  ✓ Monotonic region starts at: %d/%d (%.1f%%)\n', ...
        istart, length(theta), 100*istart/length(theta));
if ~isempty(i_first_rev)
    fprintf('  ✓ First revolution at: %d (%.3f seconds)\n', ...
            i_first_rev, timeHat(i_first_rev));
else
    fprintf('  ⚠ No complete revolution detected\n');
end
fprintf('  ✓ Output size: %d points\n', length(UHat));

% Quick validation
if length(UHat) == length(thetaHat) && length(UHat) == length(timeHat)
    fprintf('  ✓ Output dimensions consistent\n');
else
    fprintf('  ✗ Output dimension mismatch!\n');
end

if ~any(isnan(UHat))
    fprintf('  ✓ No NaN values in UHat\n');
else
    fprintf('  ⚠ UHat contains %d NaN values\n', sum(isnan(UHat)));
end

% Test with isOmega=true
fprintf('\nTesting with isOmega=true...\n');
try
    [UHat2, istart2, i_first_rev2] = shift2pi(omega, t, U, isOmega=true);
    fprintf('  ✓ isOmega=true works correctly\n');
    
    % Compare results (should be very similar)
    max_diff = max(abs(UHat - UHat2));
    if max_diff < 1e-10
        fprintf('  ✓ Results identical to position input (max diff: %.2e)\n', max_diff);
    else
        fprintf('  ⚠ Small difference from position input (max diff: %.2e)\n', max_diff);
    end
catch ME
    fprintf('  ✗ isOmega=true failed: %s\n', ME.message);
end

% Test empty signal
fprintf('\nTesting empty signal...\n');
try
    [UHat_empty, ~, ~, thetaHat_empty, timeHat_empty] = shift2pi(theta, t, []);
    if isempty(UHat_empty) && ~isempty(thetaHat_empty)
        fprintf('  ✓ Empty signal handled correctly\n');
    else
        fprintf('  ✗ Empty signal handling issue\n');
    end
catch ME
    fprintf('  ✗ Empty signal test failed: %s\n', ME.message);
end

fprintf('\n=== Test Complete ===\n');
fprintf('If you see mostly ✓ symbols above, the function is working correctly!\n');

% Optional: Create a simple plot to visualize results
try
    figure('Name', 'shift2pi Results', 'Position', [100 100 800 600]);
    
    % Plot original and shifted signals
    subplot(3,1,1);
    plot(t, U, 'b-', 'LineWidth', 1);
    hold on;
    plot(t, UHat, 'r--', 'LineWidth', 1);
    xlabel('Time [s]');
    ylabel('Signal');
    title('Original vs Shifted Signal');
    legend('Original U', 'Shifted UHat', 'Location', 'best');
    grid on;
    
    % Plot angular positions
    subplot(3,1,2);
    plot(t, theta, 'b-', 'LineWidth', 1);
    hold on;
    plot(t, thetaHat, 'r--', 'LineWidth', 1);
    xlabel('Time [s]');
    ylabel('Angular Position [rad]');
    title('Original vs Shifted Angular Position');
    legend('Original θ', 'Shifted θHat', 'Location', 'best');
    grid on;
    
    % Plot angular velocity
    subplot(3,1,3);
    omega_actual = [diff(theta)./diff(t); nan];  % Approximate derivative
    plot(t, omega_actual, 'g-', 'LineWidth', 1);
    xlabel('Time [s]');
    ylabel('Angular Velocity [rad/s]');
    title('Angular Velocity');
    grid on;
    
    fprintf('\nPlot created for visual verification.\n');
    
catch
    fprintf('Note: Could not create plot (figure window may not be available).\n');
end
