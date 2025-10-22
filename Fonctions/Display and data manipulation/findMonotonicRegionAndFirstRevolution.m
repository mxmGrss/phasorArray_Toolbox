function [ind, indFirstRev,region_starts, region_lengths, first_rev_info] = findMonotonicRegionAndFirstRevolution(theta)
%FINDMONOTONICREGIONANDFIRSTREVOLUTION Find monotonic region and first revolution completion.
%
%   [IND, INDFIRSTREV] = FINDMONOTONICREGIONANDFIRSTREVOLUTION(THETA) finds:
%   1) A monotonic increasing region THETA(IND(1):IND(2))
%   2) The exact location where a complete 2π revolution occurs within that region
%
%   This function is essential for angular domain operations requiring both
%   monotonic data (for stable interpolation) and precise revolution detection.
%
%   Inputs:
%       THETA - Angular position vector [rad] (numeric vector)
%
%   Outputs:
%       IND        - [ind_start, ind_end] indices defining monotonic region where
%                    THETA(IND(1):IND(2)) is strictly increasing
%       INDFIRSTREV - Structure with exact 2π revolution location:
%                     .k - Integer index where THETA(k) <= start_angle + 2π < THETA(k+1)
%                     .f - Fractional part [0,1] for interpolation to exact 2π
%                     Empty [] if no complete revolution found in monotonic region
%
%   Mathematical relationship:
%       start_angle + 2π = THETA(k) + f * (THETA(k+1) - THETA(k))
%       where start_angle = THETA(IND(1))
%
%   Algorithm:
%       1) Find longest monotonic increasing subsequence from any starting point
%       2) Within that region, find exact 2π revolution completion using interpolation
%       3) Optimize for maximum data utilization while ensuring monotonicity
%
%   Examples:
%       % Simple monotonic case
%       theta = linspace(0, 4*pi, 100);
%       [ind, indFirstRev] = findMonotonicRegionAndFirstRevolution(theta);
%       if ~isempty(indFirstRev)
%           exact_2pi_angle = theta(indFirstRev.k) + indFirstRev.f * ...
%                            (theta(indFirstRev.k+1) - theta(indFirstRev.k));
%           fprintf('Exact 2π completion at interpolated position: %.6f rad\n', exact_2pi_angle);
%       end
%       
%       % Non-monotonic start with oscillations
%       t = linspace(0, 6, 200);
%       theta = -0.5*sin(10*t) + 2*t;  % Small oscillations + increasing trend
%       [ind, indFirstRev] = findMonotonicRegionAndFirstRevolution(theta);
%       
%       % Plot results with exact interpolation
%       figure;
%       plot(theta, 'b-', 'DisplayName', 'Original');
%       hold on;
%       plot(ind(1):ind(2), theta(ind(1):ind(2)), 'g-', 'LineWidth', 2, ...
%            'DisplayName', 'Monotonic Region');
%       if ~isempty(indFirstRev)
%           % Plot exact 2π location
%           exact_idx = indFirstRev.k + indFirstRev.f;
%           exact_angle = theta(indFirstRev.k) + indFirstRev.f * ...
%                        (theta(indFirstRev.k+1) - theta(indFirstRev.k));
%           plot(exact_idx, exact_angle, 'ro', 'MarkerSize', 8, ...
%                'DisplayName', 'Exact 2π Revolution');
%       end
%       legend; xlabel('Sample'); ylabel('Theta [rad]');
%       title('Monotonic Region and Exact Revolution Detection');
%
%   Applications:
%       - Preparing data for angular domain interpolation
%       - Precise cycle detection in rotating machinery
%       - Phase unwrapping preprocessing
%       - Periodic signal analysis setup
%
%   See also: FIND2PIANTECEDANT, UNWRAP, DIFF, ISSORTED
%
%   PhasorArray Toolbox
%   Author: Maxime Grosso
%   Version: 1.0

    % Input validation
    arguments
        theta {mustBeNumeric, mustBeVector}
    end
    
    % Ensure column vector for consistent processing
    theta = theta(:);
    n = length(theta);
    
    % Handle edge cases
    if n < 2
        if n == 1
            ind = [1, 1];
            indFirstRev = [];
        else
            ind = [];
            indFirstRev = [];
        end
        return;
    end
    
    % Find all monotonic regions and their first revolutions in one pass
    [region_starts, region_lengths, first_rev_info] = findAllMonotonicRegions(theta);
    % Find the optimal monotonic region (with embedded revolution detection)
    [ind, firstRev] = findOptimalMonotonicRegion(theta, region_starts, region_lengths, first_rev_info);

        indFirstRev = firstRev;
        ind = ind(1:2);  % Keep only the [start, end] indices
end

function [ind, firstRev] = findOptimalMonotonicRegion(theta,region_starts, region_lengths, first_rev_info)
    %FINDOPTIMALMONOTONICREGION Find the longest monotonic increasing region
    %
    % Strategy: Find the region that maximizes data utilization while
    % ensuring strict monotonicity for stable interpolation.
    
    n = length(theta);
    
    region_count = length(region_starts);
    
    % Find the best (longest) region
    if region_count == 0
        % No monotonic regions found - fallback
        best_start = max(1, n-1);
        best_length = 1;
        best_end = n;
        
        warning('findMonotonicRegionAndFirstRevolution:NoMonotonicRegions', ...
                'No monotonic regions found in data. Using fallback region.');
    else
        % Find longest region
        [best_length, best_idx] = max(region_lengths);
        best_start = region_starts(best_idx);
        best_end = min(best_start + best_length, n);
        
        % Store the first revolution info for the selected region
        best_region_revolution = first_rev_info(best_idx);
        
        % Analyze region distribution and issue warnings if needed
        analyzeMonotonicRegions(region_count, region_starts, region_lengths, ...
                               best_start, best_length, n);
    end
    
    % Validation: ensure we have at least 2 points
    if best_end <= best_start
        % Fallback: use the last two points if nothing else works
        best_start = max(1, n-1);
        best_end = n;
        best_region_revolution = [];
        
        warning('findMonotonicRegionAndFirstRevolution:PoorMonotonicity', ...
                'Could only find %d monotonic points out of %d total. Data may be noisy.', ...
                best_length, n-1);
    end
    
    ind = [best_start, best_end];
    
    % Store revolution info for later retrieval (avoid recomputation)
    firstRev = best_region_revolution;
    
    % Quality reporting
    data_usage = (best_end - best_start + 1) / n;
    if data_usage < 0.5
        warning('findMonotonicRegionAndFirstRevolution:LowDataUsage', ...
                'Monotonic region uses only %.1f%% of available data. Consider data preprocessing.', ...
                100 * data_usage);
    end
end

function [region_starts, region_lengths, first_rev_info] = findAllMonotonicRegions(theta)
    %FINDALLMONOTONICREGIONS Find all independent monotonic increasing regions and first revolutions
    %
    % Scans through the theta vector and identifies all maximal monotonic 
    % increasing subsequences. During exploration, also detects the first 2π
    % revolution within each region for efficiency.
    %
    % Inputs:
    %   theta - Angular position vector [rad]
    %
    % Outputs:
    %   region_starts   - Start indices of each monotonic region
    %   region_lengths  - Length of each monotonic region (in samples)
    %   first_rev_info  - Cell array of first revolution info for each region
    %                     Each cell contains struct with .k, .f fields or [] if no revolution
    %
    % Algorithm: Single-pass scan that simultaneously finds regions and revolutions
    
    n = length(theta);
    differences = diff(theta);
    
    % Initialize storage for regions and revolutions
    region_starts = [];
    region_lengths = [];
    first_rev_info = struct('k', [], 'f', []);
    first_rev_info = first_rev_info([]);  % Initialize as empty struct array
    
    % Optimized search - skip overlapping regions while detecting revolutions
    start_idx = 1;
    while start_idx <= n-1
        % Find how far we can go from this starting point while staying monotonic
        % AND look for first 2π revolution during exploration
        current_length = 0;
        check_idx = start_idx;
        start_angle = theta(start_idx);
        target_2pi = start_angle + 2*pi;
        revolution_found = struct('k', [], 'f', []);
        
        while check_idx <= n-1 && differences(check_idx) > 0  % Strictly increasing
            current_length = current_length + 1;
            check_idx = check_idx + 1;
            
            % Check for first 2π revolution during region exploration
            if isempty(revolution_found) && check_idx <= n
                current_angle = theta(check_idx);
                if current_angle >= target_2pi
                    % Found first 2π crossing - compute exact interpolation
                    prev_angle = theta(check_idx - 1);
                    if current_angle > target_2pi
                        % Interpolation needed
                        delta_theta = current_angle - prev_angle;
                        if abs(delta_theta) > eps
                            f = (target_2pi - prev_angle) / delta_theta;
                        else
                            f = 0;
                        end
                        f = max(0, min(1, f));  % Clamp to [0,1]
                        
                        revolution_found = struct();
                        revolution_found.k = check_idx - 1;
                        revolution_found.f = f;
                    else
                        % Exact match
                        revolution_found = struct();
                        revolution_found.k = check_idx;
                        revolution_found.f = 0;
                    end
                end
            end
        end
        
        % Record this region if it has any length
        if current_length > 0
            region_starts(end+1) = start_idx;
            region_lengths(end+1) = current_length;
            first_rev_info(end+1) = revolution_found;  % Empty [] if no revolution found
            
            % Optimization: jump to the end of current monotonic region
            start_idx = start_idx + current_length;
        else
            % Single step if no monotonic region found at this position
            start_idx = start_idx + 1;
        end
    end
    
    % Convert to row vectors for consistency
    region_starts = region_starts(:)';
    region_lengths = region_lengths(:)';
end

function analyzeMonotonicRegions(region_count, region_starts, region_lengths, ...
                                best_start, best_length, total_points)
    %ANALYZEMONOTONICREGIONS Analyze distribution of monotonic regions and warn if ambiguous
    %
    % Provides detailed analysis when there's no clear "best" monotonic region
    
    if region_count == 0
        warning('findMonotonicRegionAndFirstRevolution:NoMonotonicRegions', ...
                'No monotonic regions found in data. Results will be unreliable.');
        return;
    end
    
    % Sort regions by length for analysis
    [sorted_lengths, sort_idx] = sort(region_lengths, 'descend');
    sorted_starts = region_starts(sort_idx);
    
    % Check if there are multiple regions of similar length
    if region_count > 1
        % Find how many regions are close to the best length
        length_threshold = 0.8;  % Regions within 80% of best length are "competitive"
        competitive_mask = sorted_lengths >= (length_threshold * best_length);
        num_competitive = sum(competitive_mask);
        
        % Issue warning and display table if choice is ambiguous
        if num_competitive > 1 && best_length < 0.7 * total_points
            warning('findMonotonicRegionAndFirstRevolution:AmbiguousChoice', ...
                    ['Multiple similar monotonic regions found (%d competitive regions). ' ...
                     'Choice may be sensitive to data ordering. See region analysis below.'], ...
                    num_competitive);
            
            % Display detailed table of regions
            displayMonotonicRegionTable(sorted_starts, sorted_lengths, ...
                                       best_start, competitive_mask, total_points);
        end
    end
    
    % Summary statistics
    total_monotonic_points = sum(region_lengths);
    monotonic_coverage = total_monotonic_points / (total_points - 1) * 100;
    
    if region_count > 5 || monotonic_coverage < 50
        fprintf('\nMonotonic Region Summary:\n');
        fprintf('  Total regions found: %d\n', region_count);
        fprintf('  Total monotonic points: %d/%d (%.1f%%)\n', ...
                total_monotonic_points, total_points-1, monotonic_coverage);
        fprintf('  Largest region: %d points (%.1f%% of data)\n', ...
                best_length, best_length/(total_points-1)*100);
        fprintf('  Selected region: indices [%d:%d]\n', ...
                best_start, best_start + best_length);
    end
end

function displayMonotonicRegionTable(sorted_starts, sorted_lengths, ...
                                    best_start, competitive_mask, total_points)
    %DISPLAYMONOTONICREGIONTABLE Display formatted table of monotonic regions
    %
    % Shows detailed breakdown when region choice is ambiguous
    
    fprintf('\n');
    fprintf('┌─────────────────────────────────────────────────────────────┐\n');
    fprintf('│                    MONOTONIC REGION ANALYSIS                │\n');
    fprintf('├──────┬──────────────┬────────┬──────────┬─────────┬─────────┤\n');
    fprintf('│ Rank │ Start Index  │ Length │ End Index│ Coverage│ Status  │\n');
    fprintf('├──────┼──────────────┼────────┼──────────┼─────────┼─────────┤\n');
    
    num_regions_to_show = min(10, length(sorted_lengths));
    
    for i = 1:num_regions_to_show
        start_idx = sorted_starts(i);
        length_val = sorted_lengths(i);
        end_idx = start_idx + length_val;
        coverage = length_val / (total_points - 1) * 100;
        
        % Determine status
        if start_idx == best_start
            status = 'SELECTED';
        elseif competitive_mask(i)
            status = 'COMPETITIVE';
        else
            status = 'SHORTER';
        end
        
        fprintf('│ %4d │ %12d │ %6d │ %8d │ %6.1f%% │ %-7s │\n', ...
                i, start_idx, length_val, end_idx, coverage, status);
    end
    
    if length(sorted_lengths) > num_regions_to_show
        fprintf('│ ...  │     ...      │  ...   │   ...    │  ...    │   ...   │\n');
    end
    
    fprintf('└──────┴──────────────┴────────┴──────────┴─────────┴─────────┘\n');
    fprintf('Note: COMPETITIVE regions are within 80%% of the longest region.\n\n');
end