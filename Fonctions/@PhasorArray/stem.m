function Tout = stem(o1,varopt)
    % STEM Generate a stem plot for one or more `PhasorArray` objects.
    %
    %   STEM(o1, varopt) generates a stem plot to visualize the phasors of one or more
    %   `PhasorArray` objects. Each input object is plotted with a distinct marker style,
    %   and the function integrates seamlessly with existing figures, axes, or tiled layouts.
    %
    %   Inputs:
    %     o1          - (Repeating, `PhasorArray`) One or more `PhasorArray` objects to be plotted.
    %
    %   Name-Value Pair Arguments:
    %     'scale'     - (char) Y-axis scale. Options:
    %                     - 'log' (default): Logarithmic scale.
    %                     - 'linear': Linear scale.
    %     'explosed'  - (logical) Determines subplot arrangement:
    %                     - true (default): Each matrix component is plotted in its own subplot.
    %                     - false: All components are combined in a single plot.
    %     'display'   - (char) Specifies which part of the phasor to display:
    %                     - 'real': Real part.
    %                     - 'imag': Imaginary part.
    %                     - 'both': Both real and imaginary parts.
    %                     - 'abs' (default): Magnitude of the phasors.
    %     'marker'    - (cell array of chars) Marker styles for each `PhasorArray`.
    %                     - Default: {"o","*","x","square","diamond","^","v",">","<"}.
    %     'side'      - (char) Determines which harmonics are displayed:
    %                     - 'both': Includes both positive and negative harmonics.
    %                     - 'oneSided' (default): Displays only non-negative harmonics.
    %     'parent'    - (graphics handle) Parent figure or axes for the plot.
    %                     - Default: gcf (current figure).
    %
    %   Outputs:
    %     Tout        - (tiledlayout object) Handle to the tiled layout used for the plot.
    %
    %   Behavior:
    %     - Multiple `PhasorArray` Objects: Each object in `o1` is plotted separately using unique markers.
    %       If more objects are provided than marker styles, markers are cycled.
    %     - Real vs. Complex Phasors: Automatically adjusts `side` to 'both' if any input contains complex values.
    %
    %   Example Usage:
    %     % Plot absolute values of multiple PhasorArray objects with linear scale
    %     A1 = PhasorArray(rand(3,3,11));
    %     A2 = PhasorArray(rand(3,3,11));
    %     Tout = stem(A1, A2, 'scale', 'linear', 'display', 'abs');
    %
    %     % Inline-marker syntax: interleave marker strings between PhasorArray objects
    %     stem(A1, 'o', A2, '*', A3, 'square');  % A1->"o", A2->"*", A3->"square"
    %     stem(A1, 'o', A2);                     % A1->"o", A2->default next marker
    %
    %     % Plot real part of a PhasorArray using custom markers
    %     A = PhasorArray(rand(4,4,11));
    %     stem(A, 'display', 'real', 'marker', {"o", "^"});
    %
    %     % Combine all components of multiple PhasorArray objects in one plot
    %     stem(A1, A2, 'explosed', false, 'side', 'both');
    %
    %   See also: stemPhasor, tiledlayout.
    arguments (Repeating)
        o1
    end
    arguments
        varopt.scale {mustBeMember(varopt.scale,{'log','linear'})}='log'
        varopt.explosed = true
        varopt.display {mustBeMember(varopt.display,{'real','imag','both','abs','absangle'})} = 'abs'
        varopt.marker = {};   % empty = not user-provided (sentinel)
        varopt.side {mustBeMember(varopt.side,{'both','oneSided'})} = 'oneSided'
        varopt.parent = gcf
        varopt.uniformYLim logical = false;
    end



    % --- Inline-marker syntax: stem(A1,'o', A2,'*', ...) ---
    % Even count: strict alternating PhasorArray / string pairs.
    % Odd count:  must be all PhasorArray (normal mode).
    if mod(numel(o1), 2) == 0
        hasInlineMarkers = all(cellfun(@(x) isa(x,'PhasorArray'), o1(1:2:end))) && ...
            all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), o1(2:2:end)));
    else
        if ~all(cellfun(@(x) isa(x,'PhasorArray'), o1))
            error('PhasorArray:stem:badInput', ...
                ['Invalid input: when an odd number of arguments is provided, ' ...
                'all must be PhasorArray objects.\n' ...
                'For inline markers use an even count: stem(A1,''o'', A2,''*'', ...).']);
        end
        hasInlineMarkers = false;
    end

    if hasInlineMarkers
        % Warn if the 'marker' name-value is also set
        if ~isempty(varopt.marker)
            warning('PhasorArray:stem:markerConflict', ...
                ['Markers were specified both inline and via the ''marker'' name-value argument. ' ...
                'Inline markers take precedence; ''marker'' is ignored.']);
        end
        varopt.marker = o1(2:2:end);
        o1 = o1(1:2:end);
    end

    % Apply default marker list if none was provided (neither inline nor name-value)
    if isempty(varopt.marker)
        varopt.marker = {"o","*","x","square","diamond","^","v",">","<"};
    end

    % Check if all PhasorArray objects in o1 are real
    if ~all(cellfun(@(x) isreal(x), o1))
        varopt.side = 'both';
    end

    if isscalar(o1{1})
        varopt.explosed = false;
    end

    if ~isa(varopt.marker, "cell")
        varopt.marker = {varopt.marker};
    end
    varhold = ishold;
    T = stemPhasor(o1{1}, scale=varopt.scale, hold=varhold, explosed=varopt.explosed, display=varopt.display, marker=varopt.marker{1}, side=varopt.side, parent=varopt.parent, uniformYLim=varopt.uniformYLim);
    n = numel(o1);
    nmarker = numel(varopt.marker);
    for n_iter = 2:n
        oi = o1{n_iter};
        ni = mod(n_iter - 1, nmarker) + 1;   % 1-based cyclic index
        hold on
        stemPhasor(oi, scale=varopt.scale, hold=true, explosed=varopt.explosed, ...
            marker=varopt.marker{ni}, display=varopt.display, parent=gca, uniformYLim=varopt.uniformYLim);
    end
    hold off
    if varhold
        hold on
    else
        hold off
    end

    if nargout>0
        Tout = T;
    end
end
