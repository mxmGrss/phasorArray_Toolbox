function phasor = computeAngularPhasors(theta, time, omega, signal, harmonics, istart, k, f, method)
%COMPUTEANGULARPHASORS Core SFT computation via cumulative integration.
%
%   phasor = computeAngularPhasors(theta, time, omega, signal, harmonics,
%                                   istart, k, f, method)
%
%   Implements the running-integral form of the angular SFT:
%     phasor(t) = [I(t) - I_shifted(t)] / 2π
%   where I(t) = ∫ base(τ)·x(τ) dτ  (cumtrapz) and I_shifted uses the
%   fractional interpolation k, f from find2piAntecedant.
%
%   ⚠ PHASE REFERENCE: theta is used as-is. If theta(1) ≠ 0, harmonic k
%   carries an absolute phase offset exp(−j·k·theta(1)).  Do NOT shift
%   theta to start at zero before calling this function.
%
%   Inputs:
%     theta     : column [N×1], unwrapped phase [rad]
%     time      : row [1×N], time [s]
%     omega     : column [N×1], pulsation [rad/s]
%     signal    : row [1×N], pre-sanitized double (no Inf/NaN)
%     harmonics : row [1×nH], harmonic orders to compute
%     istart    : scalar, first index of valid one-revolution window
%     k         : row [1×N], integer antecedent indices (find2piAntecedant)
%     f         : row [1×N], fractional antecedent parts
%     method    : 'angle' — integrate w.r.t. phase dφ
%                 'mixed' — integrate w.r.t. time, weighted by ω
%
%   Output:
%     phasor : [nH×N] complex matrix
%              Region [1 : istart-1] contains I(t)/(2π)  (no valid window).
%              Region [istart : N]   contains the full sliding window phasor.
%
%   See also: angularsft, validateAngularInputs, find2piAntecedant

nH       = numel(harmonics);
nSamples = numel(theta);

% Analysis basis [nH × nSamples]: exp(−j k θ)
baseExp = exp(-1j * harmonics(:) * theta(:)');

switch method
    case 'angle'
        base  = baseExp;
        x_int = theta;     % integrate w.r.t. phase
    case 'mixed'
        base  = baseExp .* omega(:)';
        x_int = time;      % integrate w.r.t. time
end

% Running integral: I(t) = ∫_0^t base(τ)·x(τ) dτ   [nH × nSamples]
Integral_cumul = cumtrapz(x_int, base .* signal, 2);

% Shifted integral via fractional interpolation over the stable region
% [1 : istart-1] stays zero → phasor = I(t)/2π there (partial window, v1 behaviour)
Integral_shifted = zeros(nH, nSamples);

k_safe   = min(k(istart:end), nSamples - 1);  % guard: k(i)<i always, but safe at last sample
f_region = f(istart:end);                      % [1 × M], M = nSamples - istart + 1

% Vectorized over all harmonics at once: [nH × M] matrix indexing
Integral_shifted(:, istart:end) = ...
    Integral_cumul(:, k_safe) + ...
    f_region .* (Integral_cumul(:, k_safe + 1) - Integral_cumul(:, k_safe));

phasor = (Integral_cumul - Integral_shifted) / (2*pi);

end
