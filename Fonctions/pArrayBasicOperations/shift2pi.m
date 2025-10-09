function [UHat, istart, i_first_rev, thetaHat, timeHat, L, F, IDX] = shift2pi(theta, time, U, optargin)
    %SHIFT2PI Shifts the signal U by 2*pi in the theta domain.
    %
    % This function performs angular domain signal shifting, which is useful for
    % analyzing periodic signals in rotating machinery or cyclic processes.
    %
    % T(t) is retrieved as time - timeHat where timeHat corresponds to thetaHat = theta - 2*pi.
    %
    % Syntax:
    %   [UHat, istart, i_first_rev, thetaHat, timeHat, L, F, IDX] = ...
    %       SHIFT2PI(theta, time, U, optargin)
    %
    % Inputs:
    %   theta    - Angular position vector [rad] or angular velocity if isOmega=true
    %   time     - Time vector [s]
    %   U        - Signal matrix to shift (optional, can be empty)
    %   optargin - Options structure with fields:
    %              .isOmega - logical, true if theta is angular velocity [rad/s]
    %
    % Outputs:
    %   UHat       - Signal U shifted by 2*pi in angular domain
    %   istart     - Index where monotonic increasing region starts
    %   i_first_rev- Index where shifted theta first becomes >= 0
    %   thetaHat   - Angular position vector shifted by -2*pi [rad]
    %   timeHat    - Time vector corresponding to shifted angular positions
    %   L          - Integer indices for interpolation approximation
    %   F          - Fractional parts for interpolation approximation  
    %   IDX        - Index mapping vector for shift operation
    %
    % See also ANGULARSFT, INTERP1
    
    arguments
        theta (:,1) double
        time (:,1) double
        U (:,:) double = []
        optargin.isOmega (1,1) logical = false
    end
    % Constants
    ANGULAR_SHIFT = 2*pi;
    
    % Step 1: Preprocess and validate inputs
    [theta, time, U] = preprocessInputs(theta, time, U, optargin.isOmega);
    
    % Step 2: Find the monotonic increasing region
    istart = findMonotonicRegionStart(theta);
    
    % Step 3: Compute index mapping for shifted signal
    IDX = computeShiftIndexMapping(theta, istart, ANGULAR_SHIFT);
    
    % Step 4: Generate shifted outputs
    [UHat, thetaHat, timeHat] = generateShiftedOutputs(theta, time, U, istart, ANGULAR_SHIFT);
    
    % Step 5: Find first revolution index
    i_first_rev = findFirstRevolutionIndex(thetaHat);
    
    % Step 6: Compute interpolation coefficients (optional, expensive)
    if nargout > 5
        [L, F] = computeInterpolationCoefficients(thetaHat);
    else
        L = [];
        F = [];
    end
end
%% Helper Functions

function [theta, time, U] = preprocessInputs(theta, time, U, isOmega)
    %PREPROCESSINPUTS Validate and preprocess all input data
    %
    % Handles input validation, shape normalization, and omega integration
    
    % Ensure row vectors for consistent processing
    if iscolumn(theta)
        theta = theta';
    end
    if iscolumn(time)
        time = time';
    end
    
    % Validate input dimensions
    validateInputDimensions(theta, time, U);
    
    % Normalize U to row-major format if provided
    if ~isempty(U) && size(U, 1) == numel(theta)
        U = U';
    end
    
    % Convert angular velocity to position if needed
    if isOmega
        theta = cumtrapz(time, theta);
    end
end

function validateInputDimensions(theta, time, U)
    %VALIDATEINPUTDIMENSIONS Check that input arrays have compatible sizes
    %
    % Provides detailed error messages to help users diagnose dimension problems
    
    nTheta = numel(theta);
    nTime = numel(time);
    
    % Check theta and time compatibility
    if nTheta ~= nTime
        error('shift2pi:DimensionMismatch', ...
              ['theta and time must have the same number of samples.\n' ...
               'Found: theta has %d samples, time has %d samples.\n' ...
               'Hint: Ensure both arrays represent the same time series.'], ...
               nTheta, nTime);
    end
    
    % Validate minimum data requirements
    if nTheta < 2
        error('shift2pi:InsufficientData', ...
              ['At least 2 data points are required for angular shifting.\n' ...
               'Found: %d data points.\n' ...
               'Hint: Angular domain operations require multiple samples to define motion.'], ...
               nTheta);
    end
    
    % Check U compatibility if provided
    if ~isempty(U)
        [nRows, nCols] = size(U);
        if nRows ~= nTheta && nCols ~= nTheta
            error('shift2pi:DimensionMismatch', ...
                  ['Signal matrix U must have the same number of samples as theta.\n' ...
                   'Found: theta has %d samples, U has dimensions [%d × %d].\n' ...
                   'Hint: U should be either [%d × N] or [N × %d] where N is the number of signals.'], ...
                   nTheta, nRows, nCols, nTheta, nTheta);
        end
        
        % Additional check for very large matrices (potential memory issues)
        if nRows * nCols > 1e6
            warning('shift2pi:LargeMatrix', ...
                    'Signal matrix U is very large (%d × %d). This may cause memory issues.', ...
                    nRows, nCols);
        end
    end
end

function istart = findMonotonicRegionStart(theta)
    %FINDMONOTONICREGIONSTART Find index where theta becomes monotonically increasing
    %
    % This function finds the first index from which all subsequent elements
    % of theta are strictly increasing, ensuring uniqueness for interpolation.
    %
    % The monotonic region is essential for angular domain operations because:
    % - Interpolation requires unique mapping: each angle should appear only once
    % - Non-monotonic regions create ambiguity in the shift operation
    % - Starting from a monotonic region ensures stable numerical behavior
    %
    % Algorithm: Uses forward search with minimum length requirement
    
    if numel(theta) <= 1
        istart = 1;
        return;
    end
    
    differences = diff(theta);
    
    % Forward search: find first position where all remaining differences are positive
    % This is more reliable than backward search for ensuring sufficient points
    istart = 1;  % Default to start
    
    for i = 1:length(differences)
        % Check if all remaining differences are positive (strictly increasing)
        if all(differences(i:end) > 0)
            % Ensure we have at least 2 points for interpolation
            remainingPoints = length(theta) - i;
            if remainingPoints >= 1  % Need at least 2 total points (index i and beyond)
                istart = i;
                break;
            end
        end
    end
    
    % Final validation: ensure we have enough points for interpolation
    remainingPoints = length(theta) - istart + 1;
    if remainingPoints < 2
        % Fallback: use as much data as possible, even if not perfectly monotonic
        istart = max(1, length(theta) - 1);
        remainingPoints = length(theta) - istart + 1;
        
        if remainingPoints < 2
            % Last resort: use the entire array
            istart = 1;
            warning('shift2pi:InsufficientData', ...
                    'Unable to find sufficient monotonic region. Using entire array (may be unreliable).');
        else
            warning('shift2pi:NonOptimalMonotonic', ...
                    'Using non-optimal monotonic region with %d points. Results may be less accurate.', ...
                    remainingPoints);
        end
    end
    
    % Additional diagnostic
    if istart > 1
        lostPoints = istart - 1;
        lostPercentage = 100 * lostPoints / length(theta);
        if lostPercentage > 50
            warning('shift2pi:LargeDataLoss', ...
                    'Monotonic region search discarded %.1f%% of data (%d points). Consider data quality.', ...
                    lostPercentage, lostPoints);
        end
    end
end

function IDX = computeShiftIndexMapping(theta, istart, angularShift)
    %COMPUTESHIFTINDEXMAPPING Create index mapping for angular shift operation
    %
    % Computes which original indices correspond to the shifted angular positions.
    % This mapping answers: "For each current angular position θ(t), which past 
    % index had the angular position θ(t) - 2π?"
    %
    % Physical meaning: IDX(k) tells us which past sample to use when we want
    % the signal value from 2π radians ago relative to sample k.
    
    nSamples = numel(theta);
    IDX = nan(size(theta));  % Use NaN for undefined regions (more explicit)
    
    % Extract monotonic region for interpolation
    thetaMonotonic = theta(istart:end);
    indexRange = istart:nSamples;
    
    % Validate that we have enough points for interpolation
    if length(thetaMonotonic) < 2
        warning('shift2pi:InsufficientMonotonicData', ...
                'Monotonic region has only %d point(s). Cannot perform interpolation.', ...
                length(thetaMonotonic));
        return;  % Return all NaN values
    end
    
    shiftedTheta = thetaMonotonic - angularShift;  % θ(t) - 2π
    
    % Find which past indices correspond to the shifted angular positions
    % This is the core operation: for each θ(t), find when we had θ(t) - 2π
    try
        interpolatedIndices = interp1(thetaMonotonic, indexRange, shiftedTheta ...
                                     );
        
        % Round to nearest integer indices and clamp to valid range
        % Note: extrapolation may give indices outside [1, nSamples]
        clampedIndices = clampToValidRange(round(interpolatedIndices), 1, nSamples);
        
        % Only fill the monotonic region where the mapping is well-defined
        IDX(istart:end) = clampedIndices;
        
    catch ME
        warning('shift2pi:InterpolationFailed', ...
                'Index mapping interpolation failed: %s. Using fallback approach.', ME.message);
        
        % Fallback: simple linear mapping (less accurate but functional)
        IDX(istart:end) = max(1, min(nSamples, round(linspace(1, nSamples, length(thetaMonotonic)))));
    end
    
    % The indices before istart remain NaN, indicating no valid shift mapping
end

function clampedValues = clampToValidRange(values, minVal, maxVal)
    %CLAMPTOVALIDRANGE Clamp array values to specified range
    %
    % Efficiently handles NaN values and bounds checking
    
    clampedValues = values;
    clampedValues(isnan(clampedValues)) = minVal;
    clampedValues(clampedValues < minVal) = minVal;
    clampedValues(clampedValues > maxVal) = maxVal;
end

function [UHat, thetaHat, timeHat] = generateShiftedOutputs(theta, time, U, istart, angularShift)
    %GENERATESHIFTEDOUTPUTS Create all shifted output arrays
    %
    % Performs the actual shifting operation using interpolation.
    % 
    % Physical interpretation:
    %   For each time t with angular position θ(t), this function finds:
    %   - The past time t-T(t) when θ(t-T(t)) = θ(t) - 2π  
    %   - Then UHat(t) = U(t-T(t)), where T(t) is the variable time delay
    %     needed to complete exactly one revolution (2π radians)
    %   - T(t) varies with angular velocity: faster rotation → smaller T(t)
    
    % Extract monotonic region
    thetaMonotonic = theta(istart:end);
    timeMonotonic = time(istart:end);
    
    % Validate sufficient data for interpolation
    if length(thetaMonotonic) < 2
        warning('shift2pi:InsufficientMonotonicData', ...
                'Monotonic region has only %d point(s). Cannot generate shifted outputs.', ...
                length(thetaMonotonic));
        
        % Return minimal valid outputs
        UHat = [];
        thetaHat = [];
        timeHat = [];
        return;
    end
    
    % Analyze data quality before processing
    dataQuality = analyzeAngularDataQuality(theta, time, istart);
    if dataQuality.hasWarnings
        % Report quality issues but continue processing
        for i = 1:length(dataQuality.warnings)
            warning('shift2pi:DataQuality', '%s', dataQuality.warnings{i});
        end
    end
    
    shiftedTheta = thetaMonotonic - angularShift;  % θ(t) - 2π
    
    % Shift the signal U if provided
    % This interpolation finds U at the times when θ was 2π radians earlier
    if ~isempty(U)
        try
            UMonotonic = U(:, istart:end);
            UHat = interp1(thetaMonotonic, UMonotonic', shiftedTheta,'linear')';
            % Result: UHat(t) = U(t - T(t)) where T(t) is variable revolution period
            UHat(isnan(UHat)) = 0;  % Handle NaN values gracefully
        catch ME
            warning('shift2pi:SignalInterpolationFailed', ...
                    'Signal interpolation failed: %s. Returning empty UHat.', ME.message);
            UHat = [];
        end
    else
        UHat = [];
    end
    
    % Create shifted angular and time vectors
    thetaHat = shiftedTheta;  % θ(t) - 2π
    
    try
        timeHat = interp1(thetaMonotonic, timeMonotonic, shiftedTheta,'linear');
        % Result: timeHat(t) = t - T(t), the actual past times
        timeHat(isnan(timeHat)) = timeMonotonic(1);  % Handle NaN values gracefully
    catch ME
        warning('shift2pi:TimeInterpolationFailed', ...
                'Time interpolation failed: %s. Using linear time approximation.', ME.message);
        
        % Fallback: linear time approximation
        timeHat = linspace(timeMonotonic(1), timeMonotonic(end), length(shiftedTheta))';
    end
end

function dataQuality = analyzeAngularDataQuality(theta, time, istart)
    %ANALYZEANGULARDATAQUALITY Analyze the quality of angular data for processing
    %
    % This function performs comprehensive quality analysis of the angular data
    % to identify potential issues that could affect the shifting operation.
    %
    % Returns a structure with quality metrics and warnings.
    
    dataQuality = struct();
    dataQuality.warnings = {};
    dataQuality.hasWarnings = false;
    
    % Extract monotonic region for analysis
    thetaMonotonic = theta(istart:end);
    timeMonotonic = time(istart:end);
    
    if length(thetaMonotonic) < 2
        dataQuality.warnings{end+1} = 'Insufficient data in monotonic region for quality analysis';
        dataQuality.hasWarnings = true;
        return;
    end
    
    % 1. Check angular velocity consistency
    angularVelocity = diff(thetaMonotonic) ./ diff(timeMonotonic);
    meanOmega = mean(angularVelocity);
    stdOmega = std(angularVelocity);
    
    if stdOmega / meanOmega > 0.5  % More than 50% variation
        dataQuality.warnings{end+1} = sprintf(...
            'High angular velocity variation detected (CV=%.1f%%). This may affect shift accuracy.', ...
            100 * stdOmega / meanOmega);
        dataQuality.hasWarnings = true;
    end
    
    % 2. Check for very low angular velocities (potential near-standstill)
    minOmega = min(angularVelocity);
    if minOmega < 0.1 * meanOmega
        dataQuality.warnings{end+1} = sprintf(...
            'Very low angular velocity detected (min=%.3f rad/s, mean=%.3f rad/s). Near-standstill conditions may cause interpolation issues.', ...
            minOmega, meanOmega);
        dataQuality.hasWarnings = true;
    end
    
    % 3. Check angular range coverage
    angularRange = thetaMonotonic(end) - thetaMonotonic(1);
    if angularRange < 4*pi  % Less than 2 full revolutions
        dataQuality.warnings{end+1} = sprintf(...
            'Limited angular range (%.2f radians = %.1f revolutions). Shift operation may have limited validity.', ...
            angularRange, angularRange/(2*pi));
        dataQuality.hasWarnings = true;
    end
    
    % 4. Check time sampling uniformity
    dt = diff(timeMonotonic);
    meanDt = mean(dt);
    if max(dt) > 2 * meanDt || min(dt) < 0.5 * meanDt
        dataQuality.warnings{end+1} = 'Non-uniform time sampling detected. This may affect interpolation accuracy.';
        dataQuality.hasWarnings = true;
    end
    
    % Store quality metrics for potential use
    dataQuality.metrics.angularRange = angularRange;
    dataQuality.metrics.meanAngularVelocity = meanOmega;
    dataQuality.metrics.angularVelocityCV = stdOmega / meanOmega;
    dataQuality.metrics.monotonicRegionLength = length(thetaMonotonic);
    dataQuality.metrics.monotonicRegionRatio = length(thetaMonotonic) / length(theta);
end

function [L, F] = computeInterpolationCoefficients(thetaHat)
    %COMPUTEINTERPOLATIONCOEFFICIENTS Compute integer and fractional interpolation parts
    %
    % This function computes the interpolation coefficients for advanced
    % signal processing operations. It's computationally expensive and only
    % computed when explicitly requested.
    %
    % Physical meaning: For each angular position thetaHat(k), find the 
    % interpolation coefficients to reconstruct the signal value at the
    % position thetaHat(k) - 2π using linear interpolation.
    %
    % The algorithm finds, for each point thetaHat(k), the index L(k) and 
    % fractional part F(k) such that:
    %   thetaHat(k) - 2π ≈ thetaHat(L(k)) + F(k) * (thetaHat(L(k)+1) - thetaHat(L(k)))
    %
    % Performance note: This implementation uses sequential search which is
    % efficient for monotonic data but could be optimized with binary search
    % for very large datasets.
    
    nPoints = numel(thetaHat);
    L = zeros(1, nPoints);
    F = zeros(1, nPoints);
    
    if nPoints == 0
        return;  % Handle empty input gracefully
    end
    
    currentIndex = 1;  % Start search from beginning
    angularShift = 2*pi;
    
    for k = 1:nPoints
        targetValue = thetaHat(k) - angularShift;  % θ(k) - 2π
        
        % Handle boundary conditions explicitly
        if targetValue < thetaHat(1)
            % Target is before the start of our data - use first interval
            currentIndex = 1;
            F(k) = 0;  % No interpolation needed, use boundary value
        elseif targetValue >= thetaHat(end)
            % Target is beyond the end of our data - use last interval
            currentIndex = max(1, nPoints - 1);
            F(k) = 1;  % Full interpolation to the end
        else
            % Normal case: find the bracketing interval
            currentIndex = findBracketingInterval(thetaHat, targetValue, currentIndex);
            
            % Compute fractional part for interpolation
            if currentIndex < nPoints
                deltaTheta = thetaHat(currentIndex + 1) - thetaHat(currentIndex);
                if abs(deltaTheta) > eps  % Avoid division by very small numbers
                    F(k) = (targetValue - thetaHat(currentIndex)) / deltaTheta;
                else
                    F(k) = 0;  % Handle nearly identical points
                end
            else
                F(k) = 0;  % Boundary case
            end
        end
        
        % Store the bracketing index
        L(k) = currentIndex;
    end
    
    % Ensure fractional parts are in valid range [0,1]
    F = max(0, min(1, F));
end

function idx = findBracketingInterval(thetaArray, targetValue, startIdx)
    %FINDBRACKETINGINTERVAL Find interval that brackets the target value
    %
    % Efficiently finds the index such that:
    %   thetaArray(idx) <= targetValue < thetaArray(idx+1)
    
    idx = startIdx;
    nPoints = numel(thetaArray);
    
    % Advance index while target is beyond current interval
    while idx < nPoints && thetaArray(idx + 1) <= targetValue
        idx = idx + 1;
    end
end

function i_first_rev = findFirstRevolutionIndex(thetaHat)
    %FINDFIRSTREVOLUTIONINDEX Find index where first complete revolution occurs
    %
    % This function finds when thetaHat has increased by 2*pi from its starting value.
    % This is the fundamental definition of completing one revolution.
    %
    % Input:
    %   thetaHat - Shifted angular position vector (theta - 2*pi)
    %
    % Output:
    %   i_first_rev - Index where first revolution completes (empty if none)
    
    if isempty(thetaHat) || length(thetaHat) <= 1
        i_first_rev = [];
        return;
    end
    
    % The core logic: find where we've rotated 2*pi from the start
    i_first_rev = find(thetaHat >= thetaHat(1) + 2*pi, 1, 'first');
    
    % Optional: warn if no revolution found but we have significant data
    if isempty(i_first_rev) && (thetaHat(end) - thetaHat(1)) > pi
        warning('shift2pi:IncompleteRevolution', ...
                'Angular data spans %.2f radians but no complete 2π revolution detected.', ...
                thetaHat(end) - thetaHat(1));
    end
end

%% Legacy Functions (kept for backward compatibility)

function i_inc = findFirstIncreasingIndex(A)
    %FINDFIRSTINCREASINGINDEX Legacy function - use findMonotonicRegionStart instead
    %
    % This function is kept for backward compatibility but is deprecated.
    % Use findMonotonicRegionStart for new code.
    
    warning('shift2pi:DeprecatedFunction', ...
            'findFirstIncreasingIndex is deprecated. Use findMonotonicRegionStart instead.');
    
    i_inc = findMonotonicRegionStart(A);
end

function [L, F] = findIntandFractionnalApprox(A)
    %FINDINTANDFRACTIONNALPPROX Legacy function - use computeInterpolationCoefficients instead
    %
    % This function is kept for backward compatibility but is deprecated.
    % Use computeInterpolationCoefficients for new code.
    
    warning('shift2pi:DeprecatedFunction', ...
            'findIntandFractionnalApprox is deprecated. Use computeInterpolationCoefficients instead.');
    
    [L, F] = computeInterpolationCoefficients(A);
end