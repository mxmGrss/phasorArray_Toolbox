function Tout = stem3(pA1,nvp)
    % STEM3 Generate a 3D stem plot for one or more `PhasorArray` objects.
    %
    %   STEM3(pA1, nvp) generates a 3D stem plot to visualize the phasors of one or more
    %   `PhasorArray` objects. Each input object is plotted with a distinct marker style,
    %   and the function integrates seamlessly with existing figures, axes, or tiled layouts.
    %
    %   Inputs:
    %     pA1          - (Repeating, `PhasorArray`) One or more `PhasorArray` objects to be plotted.
    %
    %   Name-Value Pair Arguments:
    %     'scale'     - (char) Y-axis scale. Options:
    %                     - 'log' (default): Logarithmic scale.
    %                     - 'linear': Linear scale.
    %     'explosed'  - (logical) Determines subplot arrangement:
    %                     - true (default): Each matrix component is plotted in its own subplot.
    %                     - false: All components are combined in a single plot.
    %     'marker'    - (cell array of chars) Marker styles for each `PhasorArray`.
    %                     - Default: {"o","*","x","square","diamond","^","v",">","<"}.
    %     'parent'    - (graphics handle) Parent figure or axes for the plot.
    %                     - Default: gcf (current figure).
    %
    %   Outputs:
    %     Tout        - (tiledlayout object) Handle to the tiled layout used for the plot.
    %
    %   Example Usage:
    %     % Plot absolute values of multiple PhasorArray objects with linear scale
    %     A1 = PhasorArray(rand(3,3,11));
    %     A2 = PhasorArray(rand(3,3,11));
    %     Tout = stem3(A1, A2, 'scale', 'linear');
    %
    %     % Combine all components of multiple PhasorArray objects in one plot
    %     stem3(A1, A2, 'explosed', false);
    %
    %   See also: quiver3, tiledlayout.
    arguments (Repeating)
        pA1
    end
    arguments
        nvp.scale {mustBeMember(nvp.scale,{'log','linear'})}='log'
        nvp.explosed = true
        nvp.marker ={"o","*","x","square","diamond","^","v",">","<"};
        nvp.parent = gcf
    end

    f = gcf;
    f.Visible = 'off';

    if ~isa(nvp.marker,"cell")
        nvp.marker= { nvp.marker};
    end

    n = numel(pA1);
    if nvp.explosed
        Tout = tiledlayout(pA1{1}.size(1),pA1{1}.size(2),'Parent',nvp.parent);
        for ii = 1:numelt(pA1{1})
            ax=nexttile;
            for jj = 1:n
                stem3Scalar(pA1{jj}{ii}, 'marker', nvp.marker{jj}, 'parent', ax, 'scale', nvp.scale);
                hold on
            end
            hold off
        end

    else
        Tout = tiledlayout(1,1,'Parent',nvp.parent);
    end

    f.Visible = 'on';
end

