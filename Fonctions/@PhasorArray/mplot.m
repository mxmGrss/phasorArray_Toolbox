function [r,t] = mplot(o1, nvp)
    % MPLOT Plot multiple `PhasorArray` objects on the same axes.
    %
    %   MPLOT(A1, A2, ..., Name, Value) is the multi-array companion to PLOT.
    %   It shares all the same name-value options but accepts multiple PhasorArray
    %   objects as repeating arguments and supports an inline linestyle syntax.
    %
    %   Inline linestyle syntax (even count):
    %     mplot(A1, '--', A2, '-.')     % explicit per-object linestyles
    %
    %   Normal syntax (odd count — all PhasorArray):
    %     mplot(A1, A2, A3)             % cycles through default linestyles
    %
    %   Name-Value Pair Arguments (all forwarded to PLOT):
    %     'T'           - Period (default: 2*pi).
    %     't'           - Time vector or [] for auto-grid (default: []).
    %     'linestyle'   - Cell array of linestyles cycled across objects.
    %                     Default: {'-','--',':','-.'}. Sentinel {} = not set.
    %     ...all other PLOT name-value arguments are accepted...
    %
    %   Outputs:
    %     r - Time-domain values from the first object (same as PLOT).
    %     t - Time vector used for evaluation.
    %
    %   Example:
    %     A1 = PhasorArray.random(2,2,5);
    %     A2 = PhasorArray.random(2,2,5);
    %     mplot(A1, '--', A2, '-.');       % inline linestyles
    %     mplot(A1, A2, 'T', 2*pi);        % same period, default styles
    %
    %   See also: plot, PhasorArray2time.
    arguments (Repeating)
        o1
    end
    arguments
        nvp.T       = 2*pi
        nvp.t       = []
        nvp.linestyle = {}     % empty = sentinel (not user-provided)
        nvp.plot    logical = true
        nvp.explosed logical = true
        nvp.hold    logical = false
        nvp.DispImag logical = []
        nvp.DispReal logical = []
        nvp.ZeroCentered logical = false
        nvp.title   = []
        nvp.GlobalYLim logical = false
        nvp.linkaxes = 'x'
        nvp.forceReal = false
        nvp.grid    = 'on'
    end

    defaultStyles = {'-'};   % solid line for all objects by default

    % --- Inline linestyle syntax: mplot(A1,'--', A2,'-.', ...) ---
    if mod(numel(o1), 2) == 0
        hasInlineStyles = all(cellfun(@(x) isa(x,'PhasorArray'), o1(1:2:end))) && ...
            all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), o1(2:2:end)));
    else
        if ~all(cellfun(@(x) isa(x,'PhasorArray'), o1))
            error('PhasorArray:mplot:badInput', ...
                ['Invalid input: when an odd number of arguments is provided, ' ...
                'all must be PhasorArray objects.\n' ...
                'For inline linestyles use an even count: mplot(A1,''--'', A2,''-.'', ...).']);
        end
        hasInlineStyles = false;
    end

    if hasInlineStyles
        if ~isempty(nvp.linestyle)
            warning('PhasorArray:mplot:linestyleConflict', ...
                ['Linestyles were specified both inline and via the ''linestyle'' name-value argument. ' ...
                'Inline linestyles take precedence; ''linestyle'' is ignored.']);
        end
        nvp.linestyle = o1(2:2:end);
        o1 = o1(1:2:end);
    end

    % Apply default linestyle list if none provided
    if isempty(nvp.linestyle)
        nvp.linestyle = defaultStyles;
    end

    n        = numel(o1);
    nstyles  = numel(nvp.linestyle);
    varhold  = ishold || nvp.hold;

    for k = 1:n
        oi  = o1{k};
        ls  = nvp.linestyle{mod(k-1, nstyles) + 1};
        isFirst = (k == 1);

        if nargout > 0 && isFirst
            [r, t] = plot(oi, nvp.T, nvp.t, ...
                'plot',         nvp.plot, ...
                'explosed',     nvp.explosed, ...
                'hold',         varhold && ~isFirst, ...
                'DispImag',     nvp.DispImag, ...
                'DispReal',     nvp.DispReal, ...
                'ZeroCentered', nvp.ZeroCentered, ...
                'title',        nvp.title, ...
                'LineStyle',     ls, ...
                'GlobalYLim',   nvp.GlobalYLim, ...
                'linkaxes',     nvp.linkaxes, ...
                'forceReal',    nvp.forceReal, ...
                'grid',         nvp.grid);
        else
            plot(oi, nvp.T, nvp.t, ...
                'plot',         nvp.plot, ...
                'explosed',     nvp.explosed, ...
                'hold',         ~isFirst || (varhold && isFirst), ...
                'DispImag',     nvp.DispImag, ...
                'DispReal',     nvp.DispReal, ...
                'ZeroCentered', nvp.ZeroCentered, ...
                'title',        nvp.title, ...
                'LineStyle',     ls, ...
                'GlobalYLim',   nvp.GlobalYLim, ...
                'linkaxes',     nvp.linkaxes, ...
                'forceReal',    nvp.forceReal, ...
                'grid',         nvp.grid);
        end
    end

    if varhold; hold on; else; hold off; end
end

