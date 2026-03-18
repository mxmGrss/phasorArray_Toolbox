function sympolar(Z_in, LineSpec, th_init)
% SYMPOLAR Polar plot with logarithmic magnitude scale.
%
%   SYMPOLAR(Z) plots complex numbers Z on a polar plot where the radial
%   coordinate is transformed using rho_scaled = log10(1 + |Z|/th).
%   This allows visualizing poles/multipliers across several orders of 
%   magnitude, extending from the origin.
%
%   SYNTAX:
%       sympolar(Z)
%       sympolar(Z, LineSpec)
%       sympolar(Z, LineSpec, th_init)
%
%   INPUTS:
%       Z_in      - Complex vector or cell array of complex vectors
%       LineSpec  - (Optional) MATLAB LineSpec string or cell array of strings
%       th_init   - (Optional) Initial threshold for the radial transition
%
%   PhasorArray Toolbox
%   Author: Maxime GROSSO
%   Date: March 2026

    arguments
        Z_in
        LineSpec = '.'
        th_init (1,1) double = 0.1
    end

    % 1. Input handling
    if ~iscell(Z_in), Z = {Z_in}; else, Z = Z_in; end
    if ~iscell(LineSpec), LineSpec = {LineSpec}; end

    ax = gca;

    % 1.1 Handle Cartesian to Polar conversion if necessary
    if ~isa(ax, 'matlab.graphics.axis.PolarAxes')
        if ishold(ax)
            warning('sympolar:cartesianAxesWithHold', 'Current axes are Cartesian while hold is on; creating a new polar figure.')
            figure;
            ax = polaraxes;
        else
            % Workaround to replace current axes with polar axes
            dummy = polarplot(0, 0);
            ax = gca;
            delete(dummy);
            ax.ColorOrderIndex = 1;
        end
    end

    was_held = ishold(ax);

    % 2. Radial threshold determination
    if was_held && isappdata(ax, 'symlog_thr')
        th = getappdata(ax, 'symlog_thr');
    else
        % Flatten for global threshold computation
        all_z =[];
        for i = 1:numel(Z)
            all_z =[all_z; Z{i}(:)];
        end
        rho_all = abs(all_z(all_z ~= 0));
        th = th_init;
        if ~isempty(rho_all)
            % Ensure 90% of data is in the log region
            while nnz(rho_all/th > 1)/numel(rho_all) < 0.9
                th = th / 10;
            end
        end
        setappdata(ax, 'symlog_thr', th);
    end

    hold(ax, 'on');

    % 3. Radial scaling and plotting
    max_rho_math = 0;
    for i = 1:length(Z)
        z = Z{i};
        rho = abs(z);
        theta = angle(z);
        
        % Log-Radial Transformation
        rho2 = log10(1 + rho/th);
        
        polarplot(ax, theta, rho2, LineSpec{min(i, length(LineSpec))}, 'MarkerSize', 12);
        max_rho_math = max(max_rho_math, max(rho(:)));
    end

    % 4. Format Radial Ticks
    format_radial_axis(ax, max_rho_math, th);
    
    if ~was_held, hold(ax, 'off'); end
end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function format_radial_axis(pax, max_rho, th)
    % Identify best powers for radial circles
    th_pwr = floor(log10(th));
    ticks = [0, 10.^(th_pwr : ceil(log10(max_rho)))];
    
    % Scale physical ticks to visual coordinates
    scaled_ticks = log10(1 + ticks/th);
    
    % Prepare labels
    labels = cell(size(ticks));
    for i = 1:length(ticks)
        if ticks(i) == 0
            labels{i} = '0';
        else
            labels{i} = sprintf('10^{%d}', round(log10(ticks(i))));
        end
    end
    
    % Apply to polar axes
    set(pax, 'RTick', scaled_ticks, 'RTickLabel', labels, 'TickLabelInterpreter', 'tex');
    rlim(pax, [0, max(scaled_ticks)*1.05]);
end