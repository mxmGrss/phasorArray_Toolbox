function symlogxy(X_in, Y_in, LineSpec, th_init)
% SYMLOGXY Plot data on a continuous bi-logarithmic scale (Symmetric Log).
%
%   SYMLOGXY(X, Y) plots data using the transformation:
%       y_scaled = sign(x) * log10(1 + |x|/th)
%   This allows visualizing data that spans several orders of magnitude, 
%   including zero and negative values, with a linear region around zero.
%
%   SYNTAX:
%       symlogxy(X, Y)
%       symlogxy(X, Y, LineSpec)
%       symlogxy(X, Y, LineSpec, th_init)
%
%   INPUTS:
%       X_in      - X-axis data (vector or cell array of vectors)
%       Y_in      - Y-axis data (vector or cell array of vectors)
%       LineSpec  - (Optional) MATLAB LineSpec string or cell array of strings
%       th_init   - (Optional) Initial threshold for the linear/log transition
%
%   PhasorArray Toolbox
%   Author: Maxime GROSSO
%   Date: March 2026

    arguments
        X_in
        Y_in
        LineSpec = '.'
        th_init (1,1) double = 1
    end
    
    % 1. Input handling (accepts Vectors OR Cell Arrays)
    if ~iscell(X_in), X = {X_in}; else, X = X_in; end
    if ~iscell(Y_in), Y = {Y_in}; else, Y = Y_in; end
    if ~iscell(LineSpec), LineSpec = {LineSpec}; end

    ax = gca;
    was_held = ishold(ax); % Store initial hold state
    
    % 2. Threshold determination (Auto-compute or retrieve from AppData if held)
    if was_held && isappdata(ax, 'symlog_thx')
        % If adding to existing plot, reuse its threshold to maintain spatial mapping
        thx = getappdata(ax, 'symlog_thx');
        thy = getappdata(ax, 'symlog_thy');
    else
        % New plot: compute optimal global thresholds
        thx = compute_auto_th(X, th_init);
        thy = compute_auto_th(Y, th_init);
        % Save to axis for future "hold on" additions
        setappdata(ax, 'symlog_thx', thx);
        setappdata(ax, 'symlog_thy', thy);
    end

    hold(ax, 'on'); % Force hold during plotting
    
    % Variables to track min/max of the new data batch
    minX_new = inf; maxX_new = -inf;
    minY_new = inf; maxY_new = -inf;

    % 3. Plot curves with transformation
    for i = 1:length(X)
        x_data = X{i}; y_data = Y{i};
        
        x2 = scale(x_data, thx);    
        y2 = scale(y_data, thy);    
        
        plot(ax, x2, y2, LineSpec{min(i, length(LineSpec))}, 'MarkerSize', 12);
        
        % Update bounds for the new data
        minX_new = min(minX_new, min(x_data(:)));
        maxX_new = max(maxX_new, max(x_data(:)));
        minY_new = min(minY_new, min(y_data(:)));
        maxY_new = max(maxY_new, max(y_data(:)));
    end
    
    % Plot center lines (only for the first plot)
    if ~was_held
        xline(0, 'k:', 'HandleVisibility', 'off');
        yline(0, 'k:', 'HandleVisibility', 'off');
    end

    % 4. Compute combined global limits (Old + New)
    if was_held
        % Re-transform current visual limits to find the math values
        old_xlim = inv_scale(xlim(ax), thx);
        old_ylim = inv_scale(ylim(ax), thy);
        
        minminX = min(minX_new, old_xlim(1)); maxmaxX = max(maxX_new, old_xlim(2));
        minminY = min(minY_new, old_ylim(1)); maxmaxY = max(maxY_new, old_ylim(2));
    else
        minminX = minX_new; maxmaxX = maxX_new;
        minminY = minY_new; maxmaxY = maxY_new;
    end

    % 5. Apply Ticks and Limits with formatting
    format_symlog_axis(ax, 'X', minminX, maxmaxX, thx);
    format_symlog_axis(ax, 'Y', minminY, maxmaxY, thy);
    
    grid(ax, 'on');
    
    % Restore hold state
    if ~was_held, hold(ax, 'off'); end
end

% =========================================================================
% INTERNAL HELPER FUNCTIONS
% =========================================================================

function th = compute_auto_th(data_cell, th_init)
    % Flatten data series for global analysis
    all_data =[];
    for i = 1:numel(data_cell)
        all_data =[all_data; data_cell{i}(:)];
    end
    
    % Ignore perfect zeros
    nz_data = all_data(all_data ~= 0);
    th = th_init;
    
    if isempty(nz_data), return; end
    
    % Reduce threshold until 90% of points are in the Log zone
    while nnz(abs(nz_data)/th > 1)/numel(nz_data) < 0.9
        th = th / 10;
    end
end

function new = scale(old, th)
    % Direct Transformation (Real -> Symlog Space)
    new = sign(old) .* log10(1 + abs(old)/th);
end

function old = inv_scale(new, th)
    % Inverse Transformation (Symlog Space -> Real)
    % Used for reading current axis limits
    old = sign(new) .* th .* (10.^abs(new) - 1);
end

function format_symlog_axis(ax, axis_name, min_val, max_val, th)
    th_pwr = floor(log10(th)); % Starting power for graduations
    
    % Targeted graduation generation
    ticks_pos =[];
    if max_val > th/10
        ticks_pos = 10.^(th_pwr : ceil(log10(max_val)));
    end
    
    ticks_neg =[];
    if min_val < -th/10
        ticks_neg = -10.^(th_pwr : ceil(log10(abs(min_val))));
    end
    
    real_ticks = sort([0, ticks_pos, ticks_neg]);
    
    % Apply spatial transformation to ticks
    scaled_ticks = scale(real_ticks, th);
    
    % Generate LaTeX strings for labels
    labels = cell(size(real_ticks));
    for i = 1:length(real_ticks)
        v = real_ticks(i);
        if v == 0
            labels{i} = '0';
        else
            sgn = ''; if v < 0, sgn = '-'; end
            labels{i} = sprintf('%s10^{%d}', sgn, round(log10(abs(v))));
        end
    end
    
    % Apply to axis with adaptive padding
    padding = 0.05 * abs(scale(max_val, th) - scale(min_val, th));
    if padding == 0, padding = 0.2; end % safety for single point
    
    if strcmpi(axis_name, 'X')
        xlim(ax,[scale(min_val, th)-padding, scale(max_val, th)+padding]);
        set(ax, 'XTick', scaled_ticks, 'XTickLabel', labels, 'TickLabelInterpreter', 'tex');
    else
        ylim(ax,[scale(min_val, th)-padding, scale(max_val, th)+padding]);
        set(ax, 'YTick', scaled_ticks, 'YTickLabel', labels, 'TickLabelInterpreter', 'tex');
    end
end