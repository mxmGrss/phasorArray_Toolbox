function [phasor_cell, theta, omega, meta, phasorStruct] = angularsft_v2(theta, time, omega, signals, harmonics, NameSignals, PlotTAPRI, nvp)
% ANGULARSFT_V2 Compute Sliding Fourier Transform for time-varying frequency signals.
%
% Refactored version of angularsft with explicit handling of fractional
% interpolation using find2piAntecedant.
%
%   [phasor_cell, theta, omega, meta, phasorStruct] = ANGULARSFT_V2(theta, time, omega, signals, harmonics, NameSignals, PlotTAPRI, options)
%
%   MATHEMATICAL FORMULATION:
%   ------------------------
%   For a signal x with time-varying fundamental frequency ω(t), the k-th harmonic phasor is:
%
%   Method 'angle' (x is a function of phase φ):
%       h_k(t) = (1/2π) ∫_{θ(t)-2π}^{θ(t)} x(φ) exp(-jkφ) dφ
%
%   Method 'mixed' (x is a function of time τ, equivalent via dφ = ω(τ)dτ):
%       h_k(t) = (1/2π) ∫_{t-T(t)}^{t} x(τ) ω(τ) exp(-jkθ(τ)) dτ
%
%   where T(t) is the time needed to complete one revolution: θ(t) - θ(t-T(t)) = 2π
%
%   INPUTS:
%       theta       - Angular phase vector [rad] (monotonically increasing after transient)
%       time        - Time vector [s] (same length as theta)
%       omega       - Instantaneous pulsation [rad/s] (empty → computed from theta/time)
%       signals     - Cell array of signals to analyze, or single vector
%       harmonics   - Harmonics to compute (default: 0:5), cell array or single vector
%       NameSignals - (Optional) Cell array of signal names for plotting
%       PlotTAPRI   - [T A P R I] logical vector for plots:
%                     T: Time-domain signal, A: Abs(phasor), P: Phase, R: Real, I: Imag
%                     Default: [true true false false false]
%       options     - Name-value arguments:
%           .xAxes       - 'time' | 'phase' | 'revolution' (default: 'time')
%           .method      - 'angle' | 'mixed' (default: 'angle')
%           .orientation - 'ver' | 'hor' (default: 'hor')
%           .plotDebut   - Plot starting transient region (default: true)
%           .plotOmega   - Plot omega profile (default: false)
%           .plotlang    - 'fr' | 'en' (default: 'fr')
%           .interpMethod- 'find2pi' | 'interp1' (default: 'find2pi')
%
%   OUTPUTS:
%       phasor_cell  - Cell array of complex phasor matrices [nHarmonics × nSamples]
%       theta        - Processed phase vector (unwrapped, normalized)
%       omega        - Processed pulsation vector
%       meta         - Structure with diagnostic info:
%                      .istart      - Index where monotonic region starts
%                      .k           - Integer indices for θ-2π interpolation
%                      .f           - Fractional parts for interpolation
%                      .nRevolutions- Total number of revolutions
%       phasorStruct - Structure array for plotting (legacy format)
%
%   DIFFERENCES FROM angularsft (v1):
%       - Uses find2piAntecedant by default (faster, explicit fractional interpolation)
%       - Cleaner separation: validation → computation → plotting
%       - Explicit meta output with diagnostic info
%       - Option to fallback to interp1 for compatibility testing
%       - Improved documentation with mathematical formulas
%
%   EXAMPLE:
%       % Generate test signal with varying frequency
%       t = linspace(0, 2, 1000);
%       omega = 2*pi * (5 + 3*sin(2*pi*t));
%       theta = cumtrapz(t, omega);
%       x = cos(theta) + 0.5*cos(3*theta);
%       
%       % Compute phasors
%       [phasors, ~, ~, meta] = angularsft_v2(theta, t, omega, x, 0:5);
%       
%       % Plot results
%       figure; plot(t, abs(phasors{1}));
%
%   See also: angularsft, find2piAntecedant, shift2pi, plotAngularSFT

arguments
    theta
    time
    omega
    signals
    harmonics = 0:5
    NameSignals = {}
    PlotTAPRI {mustBeNumericOrLogical} = [true true false false false]
    nvp.xAxes {mustBeMember(nvp.xAxes, {'time', 'phase', 'revolution'})} = 'time'
    nvp.method {mustBeMember(nvp.method, {'angle', 'mixed'})} = 'angle'
    nvp.orientation {mustBeMember(nvp.orientation, {'ver', 'hor'})} = 'hor'
    nvp.plotDebut logical = true
    nvp.plotOmega logical = false
    nvp.plotlang = 'fr'
    nvp.interpMethod {mustBeMember(nvp.interpMethod, {'find2pi', 'interp1'})} = 'find2pi'
end

%% STEP 1: Validate and normalize inputs
[theta, time, omega, signals, harmonics, NameSignals, PlotTAPRI] = ...
    validateInputs_v2(theta, time, omega, signals, harmonics, NameSignals, PlotTAPRI);

%% STEP 2: Find monotonic region and compute interpolation indices
switch nvp.interpMethod
    case 'find2pi'
        % Use find2piAntecedant (optimized, explicit fractional interpolation)
        [k, f] = find2piAntecedant(theta);
        
        % Find where valid interpolation starts (k > 1 means theta(k) exists for theta-2π)
        istart = find(k > 1, 1, 'first');
        if isempty(istart)
            istart = 1;
            warning('angularsft_v2:InsufficientRevolutions', ...
                'Signal does not span 2π. Results may be unreliable.');
        end
        
    case 'interp1'
        % Fallback to shift2pi + interp1 (for compatibility testing)
        [~, istart, ~, ~, ~, k, f, ~] = shift2pi(theta, time);
        
        % k and f from shift2pi are for the entire signal, extract monotonic region
        if istart > 1
            k = k(istart:end);
            f = f(istart:end);
        end
end

%% STEP 3: Compute phasors for each signal
nSignals = numel(signals);
phasor_cell = cell(nSignals, 1);

for ii = 1:nSignals
    % Extract and sanitize signal
    Signal_ii = double(signals{ii});
    if ~isrow(Signal_ii)
        Signal_ii = Signal_ii.';
    end
    Signal_ii(isinf(Signal_ii)) = 0;
    
    % Extract harmonics
    Hmcs_ii = harmonics{ii};
    if ~isrow(Hmcs_ii)
        Hmcs_ii = Hmcs_ii';
    end
    
    % Compute phasors
    phasor_ii = computePhasors_v2(theta, time, omega, Signal_ii, Hmcs_ii, ...
                                   istart, k, f, nvp.method, nvp.interpMethod);
    
    % Store results
    phasor_cell{ii} = phasor_ii;
    
    % Build phasorStruct for plotting (legacy format)
    phasorStruct(ii).name = NameSignals{ii};
    phasorStruct(ii).phasors = phasor_ii;
    phasorStruct(ii).signal = signals{ii};
    phasorStruct(ii).time = time;
    phasorStruct(ii).harmonics = Hmcs_ii;
    phasorStruct(ii).theta = theta;
    phasorStruct(ii).omega = omega;
    phasorStruct(ii).IDX = k;  % For legacy compatibility
end

%% STEP 4: Build metadata output
meta = struct();
meta.istart = istart;
meta.k = k;
meta.f = f;
meta.nRevolutions = theta(end) / (2*pi);
meta.method = nvp.method;
meta.interpMethod = nvp.interpMethod;

%% STEP 5: Plot if requested
if sum(PlotTAPRI) > 0
    plotAngularSFT(phasorStruct, PlotTAPRI, ...
                   "orientation", nvp.orientation, ...
                   "plotDebut", nvp.plotDebut, ...
                   "plotOmega", nvp.plotOmega, ...
                   "xAxes", nvp.xAxes, ...
                   "lang", nvp.plotlang);
end

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [theta, time, omega, signals, harmonics, NameSignals, PlotTAPRI] = ...
    validateInputs_v2(theta, time, omega, signals, harmonics, NameSignals, PlotTAPRI)
%VALIDATEINPUTS_V2 Validate and normalize all inputs
%
% This function performs comprehensive input validation and normalization:
%   - Converts signals to cell array if needed
%   - Computes omega from theta if missing
%   - Computes theta from omega if missing
%   - Unwraps theta
%   - Ensures consistent dimensions (time: row, theta/omega: column)

    % Normalize signals to cell array
    if ~iscell(squeeze(signals))
        signals = squeeze(signals);
        if ~isvector(signals)
            if size(signals, 1) == numel(time)
                signals = num2cell(signals, 1);
            elseif size(signals, 2) == numel(time)
                signals = num2cell(signals, 2);
            else
                error('angularsft_v2:DimensionMismatch', ...
                      'Signal dimensions do not match time vector.');
            end
        else
            signals = {squeeze(signals)};
        end
    end
    
    % Normalize harmonics to cell array
    if ~iscell(squeeze(harmonics))
        harmonics = {squeeze(harmonics)};
    end
    
    % Replicate harmonics if single element provided
    if numel(harmonics) == 1
        harmonics = repmat(harmonics, 1, numel(signals));
    end
    
    % Generate default signal names
    if isempty(NameSignals)
        NameSignals = arrayfun(@(i) sprintf('signal %d', i), 1:numel(signals), 'UniformOutput', false);
    end
    if ~iscell(NameSignals)
        NameSignals = {NameSignals};
    end
    
    % Handle empty theta
    if isempty(theta)
        theta = 0;
    end
    
    % Compute omega if missing
    if isempty(omega) || numel(omega) == 0
        if numel(theta) ~= numel(time)
            error('angularsft_v2:MissingOmega', ...
                  'omega cannot be empty if theta is not a time vector.');
        end
        omega = gradient(theta) ./ gradient(time);
    elseif numel(omega) == 1
        omega = ones(size(time)) * omega;
    end
    
    % Compute theta if scalar (integrate omega)
    if numel(theta) == 1
        theta = cumtrapz(time, omega) + theta;
    end
    
    % Unwrap theta
    theta = unwrap(theta);
    
    % Ensure consistent dimensions
    omega = omega(:);  % Column vector
    time = time(:)';   % Row vector
    theta = theta(:);  % Column vector
    
    % Validate dimensions
    if numel(theta) ~= numel(time) || numel(omega) ~= numel(time)
        error('angularsft_v2:DimensionMismatch', ...
              'theta, time, and omega must have the same length.');
    end
end

function phasor = computePhasors_v2(theta, time, omega, signal, harmonics, istart, k, f, method, interpMethod)
%COMPUTEPHASORS_V2 Core phasor computation with explicit fractional interpolation
%
% Computes the sliding Fourier transform using cumulative integration and
% shifted integral subtraction.
%
% ALGORITHM:
%   1. Compute base: exp(-jkθ) or exp(-jkθ)·ω depending on method
%   2. Cumulative integral: I(t) = ∫_0^t base(τ)·signal(τ) dτ
%   3. Shifted integral: I_shifted(t) = I(t-T(t)) where θ(t-T(t)) = θ(t)-2π
%   4. Phasor: h(t) = (I(t) - I_shifted(t)) / 2π

    nH = numel(harmonics);
    nSamples = numel(theta);
    
    % Build analysis base
    baseExp = exp(-1j * harmonics' * theta');  % [nH × nSamples]
    
    switch method
        case 'angle'
            % Integrate w.r.t. phase: dφ
            base = baseExp;
            x_int = theta;
        case 'mixed'
            % Integrate w.r.t. time: ω(τ) dτ
            base = baseExp .* omega';
            x_int = time;
    end
    
    % Cumulative integral
    Integral_cumul = cumtrapz(x_int, base .* signal, 2);
    
    % Shifted integral via interpolation
    Integral_shifted = zeros(size(Integral_cumul));
    
    switch interpMethod
        case 'find2pi'
            % Explicit fractional interpolation: I(k) + f·(I(k+1) - I(k))
            for h = 1:nH
                Int_h = Integral_cumul(h, :);
                
                % Vectorized interpolation over monotonic region
                k_region = k(max(1, istart):end);
                f_region = f(max(1, istart):end);
                
                % Ensure k+1 does not exceed array bounds
                k_safe = min(k_region, nSamples - 1);
                
                Integral_shifted(h, istart:end) = Int_h(k_safe) + ...
                    f_region .* (Int_h(min(k_safe + 1, nSamples)) - Int_h(k_safe));
            end
            
        case 'interp1'
            % Standard MATLAB interpolation (for validation)
            Integral_shifted(:, istart:end) = interp1(theta(istart:end)', ...
                Integral_cumul(:, istart:end).', theta(istart:end)' - 2*pi, 'linear', 0).';
    end
    
    % Handle NaNs
    Integral_shifted(isnan(Integral_shifted)) = 0;
    
    % Final phasor: difference / 2π
    phasor = (Integral_cumul - Integral_shifted) / (2*pi);
end
