function [Adetph,Adet_t] = PhasorDet(Aph,varg)
    %PHASORDET Computes the phasors of the determinant of a time-varying matrix A(t)
    %
    %   [Adetph, Adet_t] = PhasorDet(Aph, varg) calculates the phasors for the
    %   determinant of a matrix A(t) represented in phasor form, where the input
    %   matrix A(t) is assumed to be periodic in nature. The determinant is computed
    %   for each time slice, and the result is transformed back into phasor form.
    %
    %   Inputs:
    %       - Aph: The 3D array representing the phasors of A(t). The third dimension
    %         corresponds to the harmonic components of A(t).
    %       - varg: Optional arguments, provided as a structure or name-value pairs:
    %           - nT: The number of time steps (default is 1).
    %           - T: The period of the matrix A(t) (default is 1).
    %           - m: The truncation order (default is calculated from the harmonic length).
    %           - plot: Boolean flag to plot results (default is false).
    %           - reduceThreshold: Threshold for reducing the phasors (default is 1e-15).
    %           - reduceMethod: The method to use for reduction (default is 'relative').
    %           - autoTrunc: Flag to automatically truncate based on filtered differences (default is false).
    %
    %   Outputs:
    %       - Adetph: The phasor representation of the determinant of A(t) in 3D array form.
    %       - Adet_t: The time-domain representation of the determinant.
    %
    %   Notes:
    %       - If `Aph` is a `PhasorArray`, it is first converted to its value form.
    %       - The determinant of each time slice of A(t) is computed using the `det` function.
    %       - If `autoTrunc` is set to true, the phasors are truncated based on filtered differences.
    %
    %   Example usage:
    %       [Adetph, Adet_t] = PhasorDet(Aph, 'nT', 10, 'T', 1, 'm', 8, 'plot', true);
    %
    %   See also: PhasorArray, ReduceArray, PhasorArray2time, TimeArray2Phasors

arguments
Aph
varg.nT=1
varg.T=1
varg.m=[]
varg.plot=false
varg.reduceThreshold = 1e-15
varg.reduceMethod = 'relative'
varg.autoTrunc = false
end

Aph = pvalue(Aph);

nT=varg.nT;
T=varg.T;
m=varg.m;

hA=(size(Aph,3)-1)/2;

if isempty(m)
m=max(ceil(log2(hA+1)+1),8);
end

if hA==0
    Adetph=det(Aph);
    return
end

n=2^m;
t=0:T/n:nT*T-T/n;
if isa(Aph,"ndsdpvar") || isa(Aph,"sdpvar")
    Aph=value(Aph);
end

Aph=ReduceArray(Aph);

At=PhasorArray2time(Aph,T,t,"plot",false);

Adet_t=arrayfun(@(k) det(At(:,:,k)), 1:size(At,3));
Adet_t=permute(Adet_t,[1 3 2]);



Adetphcomp=TimeArray2Phasors(Adet_t,nT);

% --- Auto-Truncation Logic (Heuristic) ---
% Goal: Detect the noise floor where spectral magnitude stops decreasing
% and becomes dominated by numerical noise.
h = (numel(Adetphcomp)-1)/2;
magnitude = squeeze(abs(Adetphcomp(1,1,(h+1):end)));

% 1. Convert to log-scale (dB-like)
% Protect against log10(0) -> -Inf by clamping to eps
safe_magnitude = max(magnitude, eps); 
log10Ph = log10(safe_magnitude);

% 2. Estimate Noise Floor (from tail of the spectrum)
% Use last 10 samples or all if length < 10
n_tail = min(10, numel(log10Ph));
if n_tail > 0
    noise_floor = mean(log10Ph(end-n_tail+1:end));
else
    noise_floor = -15; % Default fallback
end

% Replace actual zeros (which became eps -> -16) with noise floor for cleaner diff
log10Ph(magnitude < 10*eps) = noise_floor;

% 3. Compute Spectral Slope (Derivative of log-magnitude)
% Positive slope means noise floor reached (or oscillation)
diffPh = diff(log10Ph, 1);

% 4. Smooth the slope to avoid local false positives
% Use lowpass filter if available/robust, else simple moving average
try
    % Try standard MATLAB lowpass with steepness control
    filt_diff = lowpass(diffPh, 0.05); 
catch
    % Fallback: moving average (robust and standard)
    filt_diff = movmean(diffPh, 5);
end

if ~varg.autoTrunc
    % Standard reduction based on threshold
    Adetph = ReduceArray(Adetphcomp, 'reduceMethod', varg.reduceMethod, 'reduceThreshold', varg.reduceThreshold);
else
    % Auto-truncation: Cut where the smoothed slope turns positive (noise floor)
    % Find first positive slope + buffer
    idx_cutoff = find(filt_diff > 0, 1);
    if isempty(idx_cutoff)
        idx_cutoff = numel(filt_diff); % Keep all if no noise floor found
    end
    % Add safety margin (+5 harmonics)
    trunc_h = min(idx_cutoff + 5, h);
    Adetph = ReduceArray(Adetphcomp, trunc_h);
end

if varg.plot
    clf
    tiledlayout("flow")
    nexttile
    h=(numel(Adetph)-1)/2;
    plot(t,squeeze(reshape(Adet_t,[],1,numel(t))))
    hold on
    Ait=PhasorArray2time(Adetph,T,t,"plot",false);
    plot(t,squeeze(reshape(Ait,[],1,numel(t))),'--')
    hold off
    legend('real det','reconstructed det from phasor')
    title('det(A)(t)')
    nexttile
    stem(0:h,squeeze(abs(Adetph(1,1,(h+1):end))))
    set(gca,'YScale','log')
    title('Stem of abs of phasor of det(A)')
    if varg.autoTrunc
    nexttile
    plot(diffPh)
    hold on 
    plot(filt_diff)
    fprintf("recommended truncature order is %d\n",find(filt_diff>0,1))
    title('filtered diff of log10(abs(phasor))')
    end
    
    nexttile
    h=(numel(Adetphcomp)-1)/2;
    stem(0:h,squeeze(abs(Adetphcomp(1,1,(h+1):end))))
    set(gca,'YScale','log')
    title('Stem of abs of phasor of det(A) before trunc')
end

