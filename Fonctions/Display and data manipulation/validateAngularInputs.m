function [theta, time, omega, signals, harmonics, NameSignals] = validateAngularInputs(theta, time, omega, signals, harmonics, NameSignals)
%VALIDATEANGULARINPUTS Validate and normalize inputs for angularsft.
%
%   [theta, time, omega, signals, harmonics, NameSignals] = ...
%       validateAngularInputs(theta, time, omega, signals, harmonics, NameSignals)
%
%   After this function:
%     theta      : column [N×1], unwrapped, non-decreasing
%     time       : row [1×N]
%     omega      : column [N×1] [rad/s]
%     signals    : cell {1×p} of row double vectors [1×N]
%     harmonics  : cell {1×p}, one harmonic vector per signal
%     NameSignals: cell {1×p} of char
%
%   See also: angularsft, computeAngularPhasors

arguments
    theta
    time
    omega
    signals
    harmonics = 0:5
    NameSignals = {}
end

validateattributes(time,{'numeric'},{'vector','real','finite'});
if numel(time)<2 || any(diff(time)<=0)
    error('angularsft:InvalidTime','At least two strictly increasing time samples are required.');
end

% --- Normalize signals to cell array ---
if ~iscell(squeeze(signals))
    signals = squeeze(signals);
    if ~isvector(signals)
        if size(signals, 1) == numel(time)
            signals = num2cell(signals, 1);
        elseif size(signals, 2) == numel(time)
            signals = num2cell(signals, 2);
        else
            error('angularsft:invalidSignalDimension', ...
                  'Signal size [%d×%d] does not match time vector length (%d).', ...
                  size(signals,1), size(signals,2), numel(time));
        end
    else
        signals = {squeeze(signals)};
    end
end

% --- Normalize harmonics to cell array and replicate if scalar ---
if ~iscell(squeeze(harmonics))
    harmonics = {squeeze(harmonics)};
end
if isscalar(harmonics)
    harmonics = repmat(harmonics, 1, numel(signals));
end

% --- Default signal names ---
if isempty(NameSignals)
    NameSignals = arrayfun(@(i) sprintf('signal %d', i), 1:numel(signals), 'UniformOutput', false);
end
if ~iscell(NameSignals)
    NameSignals = {NameSignals};
end
if numel(harmonics)~=numel(signals) || numel(NameSignals)~=numel(signals)
    error('angularsft:SignalCount','Provide one harmonic vector and name per signal.');
end
for ii=1:numel(signals)
    validateattributes(signals{ii},{'numeric'},{'vector','numel',numel(time)});
    validateattributes(harmonics{ii},{'numeric'},{'vector','real','finite','integer'});
end

% --- Handle scalar/empty theta ---
if isempty(theta)
    theta = 0;
end

% --- Compute omega from theta/time if missing ---
if isempty(omega)
    if numel(theta) ~= numel(time)
        error('angularsft:missingOmega', ...
              'omega cannot be empty when theta length differs from time.');
    end
    omega = gradient(theta(:)) ./ gradient(time(:));
elseif isscalar(omega)
    omega = ones(numel(time), 1) * double(omega);
end

% --- Integrate theta from omega if theta is scalar ---
if isscalar(theta)
    theta = cumtrapz(time(:), omega(:)) + double(theta);
end

% --- Unwrap theta ---
theta = unwrap(theta);

% --- Enforce consistent dimensions ---
omega = omega(:);   % column
time  = time(:)';   % row
theta = theta(:);   % column

if numel(theta) ~= numel(time) || numel(omega) ~= numel(time)
    error('angularsft:dimensionMismatch', ...
          'theta, time, and omega must all have the same length.');
end

end
