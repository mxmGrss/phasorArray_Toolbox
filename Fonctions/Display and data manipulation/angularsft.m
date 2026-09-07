function [phasor_cell, theta, omega, IDX, phasorStruct] = angularsft(theta, time, omega, signals, harmonics, NameSignals, PlotTAPRI, optarg)
%ANGULARSFT Sliding Fourier Transform for signals with time-varying frequency.
%
%   [phasor_cell, theta, omega, IDX, phasorStruct] = angularsft(theta, time,
%       omega, signals, harmonics, NameSignals, PlotTAPRI, options)
%
%   MATHEMATICAL FORMULATION
%   ------------------------
%   For a signal x with time-varying pulsation ω(t), the k-th harmonic phasor is:
%
%     Method 'angle':   h_k(t) = (1/2π) ∫_{θ(t)−2π}^{θ(t)} x(φ) exp(−jkφ) dφ
%     Method 'mixed':   h_k(t) = (1/2π) ∫_{t−T(t)}^{t} x(τ)ω(τ)exp(−jkθ(τ)) dτ
%
%   where T(t) satisfies θ(t) − θ(t−T(t)) = 2π (duration of one revolution).
%
%   ⚠ PHASE REFERENCE: theta is used as-is throughout the computation.
%   If theta(1) ≠ 0, harmonic k is rotated by exp(−j·k·theta(1)) relative
%   to a zero-phase reference. This is a convention, not a bug. Never
%   subtract theta(1) from theta before calling this function.
%
%   INPUTS
%   ------
%   theta       : Phase vector [rad]. Unwrapped, non-decreasing in stable region.
%                 Empty → integrated from omega. Scalar → initial phase offset.
%   time        : Time vector [s], same length as theta.
%   omega       : Pulsation [rad/s]. Empty → gradient(theta)/gradient(time).
%                 Scalar → constant pulsation.
%   signals     : Cell array of signals, or single vector. Each element [1×N].
%   harmonics   : Harmonic orders to track (default: 0:5). Scalar or vector →
%                 applied to all signals. Cell array → one vector per signal.
%   NameSignals : (optional) Cell array of signal names for plots.
%   PlotTAPRI   : (optional) 5-element logical [T A P R I]:
%                   T: time-domain, A: |phasor|, P: phase, R: real, I: imag.
%                   Default: [true true false false false].
%   options     : Name-value arguments:
%     .xAxes       — 'time' | 'phase' | 'revolution'  (default: 'time')
%     .method      — 'angle' | 'mixed'                (default: 'angle')
%     .orientation — 'ver' | 'hor'                    (default: 'hor')
%     .plotDebut   — show transient region             (default: true)
%     .plotOmega   — show ω profile                   (default: false)
%     .plotlang    — 'fr' | 'en'                      (default: 'fr')
%
%   OUTPUTS
%   -------
%   phasor_cell  : Cell {nSignals×1} of [nH×N] complex phasor matrices.
%                  Rows = harmonics, columns = time samples.
%                  Columns [1 : istart-1]: partial window (not a full revolution).
%                  Columns [istart : N]:   full sliding-window phasor.
%   theta        : Processed phase vector (column, unwrapped).
%   omega        : Processed pulsation vector (column).
%   IDX          : Integer antecedent indices k from find2piAntecedant [1×N].
%   phasorStruct : Struct array for plotAngularSFT.
%                  Fields: name, phasors, signal, time, harmonics, theta, omega,
%                          IDX, meta (with .istart, .nRevolutions).
%
%   EXAMPLE
%   -------
%   t     = linspace(0, 4, 2000)';
%   om    = 2*pi*(10 + 2*sin(2*pi*t));
%   theta = cumtrapz(t, om);
%   x     = cos(theta) + 0.4*cos(3*theta);
%   [p, ~, ~, ~, S] = angularsft(theta, t, om, x, [1 3]);
%   figure; plot(t, abs(p{1}))
%
%   See also: computeAngularPhasors, validateAngularInputs, find2piAntecedant,
%             plotAngularSFT

arguments
    theta
    time
    omega
    signals
    harmonics   = 0:5
    NameSignals = {}
    PlotTAPRI {mustBeNumericOrLogical} = [true true false false false]
    optarg.xAxes {mustBeMember(optarg.xAxes, {'time','phase','revolution'})} = 'time'
    optarg.method {mustBeMember(optarg.method, {'angle','mixed'})} = 'angle'
    optarg.orientation {mustBeMember(optarg.orientation, {'ver','hor'})} = 'hor'
    optarg.plotDebut logical = true
    optarg.plotOmega logical = false
    optarg.plotlang = 'fr'
end

% ── 1. Validate and normalise ─────────────────────────────────────────────────
[theta, time, omega, signals, harmonics, NameSignals] = ...
    validateAngularInputs(theta, time, omega, signals, harmonics, NameSignals);

% ── 2. Find θ−2π antecedents (O(n), full vector including transient) ──────────
[k, f] = find2piAntecedant(theta);
IDX = k;

istart = find(k > 1, 1, 'first');
if isempty(istart)
    istart = 1;
    warning('angularsft:InsufficientRevolutions', ...
            'Signal does not span 2π. Phasor output will be unreliable.');
end

% ── 3. Compute phasors, signal by signal ─────────────────────────────────────
nSignals = numel(signals);

% Preallocate struct array to avoid AGROW
phasorTemplate = struct('name','', 'phasors',[], 'signal',[], 'time',[], ...
                        'harmonics',[], 'theta',[], 'omega',[], 'IDX',[], 'meta',[]);
phasorStruct = repmat(phasorTemplate, 1, nSignals);
phasor_cell  = cell(nSignals, 1);

for ii = 1:nSignals
    try
        sig_ii = double(signals{ii});
        if ~isrow(sig_ii)
            sig_ii = sig_ii.';
        end
        sig_ii(isinf(sig_ii)) = 0;

        hm_ii = harmonics{ii};
        if ~isrow(hm_ii)
            hm_ii = hm_ii.';
        end

        phasor_ii = computeAngularPhasors(theta, time, omega, sig_ii, hm_ii, ...
                                           istart, k, f, optarg.method);

        phasor_cell{ii}             = phasor_ii;
        phasorStruct(ii).name       = NameSignals{ii};
        phasorStruct(ii).phasors    = phasor_ii;
        phasorStruct(ii).signal     = signals{ii};
        phasorStruct(ii).time       = time;
        phasorStruct(ii).harmonics  = hm_ii;
        phasorStruct(ii).theta      = theta;
        phasorStruct(ii).omega      = omega;
        phasorStruct(ii).IDX        = k;
        phasorStruct(ii).meta       = struct( ...
            'istart',       istart, ...
            'nRevolutions', (theta(end) - theta(1)) / (2*pi));

    catch ME
        warning('angularsft:signalError', 'Signal "%s" skipped: %s', ...
                NameSignals{ii}, ME.message);
        phasor_cell{ii} = [];
    end
end

% ── 4. Plot if requested ──────────────────────────────────────────────────────
if sum(PlotTAPRI) > 0
    plotAngularSFT(phasorStruct, PlotTAPRI, ...
                   "orientation", optarg.orientation, ...
                   "plotDebut",   optarg.plotDebut, ...
                   "plotOmega",   optarg.plotOmega, ...
                   "xAxes",       optarg.xAxes, ...
                   "lang",        optarg.plotlang);
end

end
