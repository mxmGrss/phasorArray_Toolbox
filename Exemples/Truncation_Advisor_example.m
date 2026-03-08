%% Truncation Advisor - Quick Example
% This script demonstrates how to get practical truncation settings
% for Lyapunov, Riccati and LMI workflows.

clear; clc;

fprintf('=== Truncation Advisor Example ===\n');

% 1) Default recommendation for Lyapunov
advLyap = truncationAdvisor("lyap", "profile", "balanced");

% 2) Higher-accuracy recommendation for Riccati
advRic = truncationAdvisor("riccati", "profile", "high-accuracy", "verbose", false);

% 3) LMI-oriented split orders (hLMI, hP, hA)
advLmi = truncationAdvisor("lmi", "profile", "balanced", "verbose", false);

fprintf('Lyap  : h0=%d, hmax=%d, target=%.1e\n', ...
    advLyap.recommended.h0, advLyap.recommended.hmax, advLyap.recommended.targetResidual);
fprintf('Riccati: h0=%d, hmax=%d, target=%.1e\n', ...
    advRic.recommended.h0, advRic.recommended.hmax, advRic.recommended.targetResidual);
fprintf('LMI   : hLMI=%d, hP=%d, hA=%d\n', ...
    advLmi.recommended.hLMI, advLmi.recommended.hP, advLmi.recommended.hA);

% Example usage with existing solver API (illustrative):
% [X,res] = lyap(A,Q,'h',advLyap.solverOptions.h,...
%                   'autoUpdateh',advLyap.solverOptions.autoUpdateh,...
%                   'thresholdResidual',advLyap.solverOptions.thresholdResidual);

fprintf('Done.\n');
