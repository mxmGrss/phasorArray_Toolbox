function  T = stemPhasor(Aph,varopt)
%STEMPHASOR Plot phasors as stem plots to visualize their amplitudes against their order.
%
%   T = stemPhasor(Aph) generates a tiled layout of stem plots to visualize the amplitudes of
%   phasors against their order. The input `Aph` is a phasor array of size `N×M×(2nh+1)`, where
%   `N` and `M` are the dimensions of the matrix and `nh` is the half-length of the harmonic
%   expansion. The function returns the handle `T` to the tiled layout used for the plot.
%
%   This function generates a stem plot to visualize the amplitudes (or other parts) of 
%   phasors against their order, providing an intuitive representation of harmonic components.
%   The phasors can be plotted in a single combined plot or as individual plots for each 
%   component of the matrix.
%
%   Input Arguments:
%   - Aph (array or PhasorArray): The input phasor array, either as a `PhasorArray` object or 
%     a numeric array of size `N×M×(2nh+1)`.
%
%   Name-Value Pair Arguments:
%   - `'scale'`    (char, optional): Y-axis scale, either `'log'` (default) or `'linear'`.
%   - `'hold'`     (logical, optional): Hold state of the plot. Default: current hold state.
%   - `'explosed'` (logical, optional): If true (default), each component of `A` gets its own plot.
%   - `'display'`  (char, optional): Part of the phasor to display. Options:
%       - `'real'`: Real part of the phasors.
%       - `'imag'`: Imaginary part of the phasors.
%       - `'both'`: Both real and imaginary parts.
%       - `'abs'` (default): Magnitudes of the phasors.
%   - `'marker'`   (char, optional): Marker style for the stem plot. Default: `'o'`.
%   - `'side'`     (char, optional): Side of the plot:
%       - `'both'` (default): Displays all harmonics (negative and positive).
%       - `'oneSided'`: Displays only non-negative harmonics.
%   - `'parent'`   (graphics handle, optional): Parent figure or axes for the plot. Default: current figure.
%
%   Output Arguments:
%   - T (tiledlayout object): Handle to the tiled layout used for the plot.
%
%
%   Key Features:
%   - Allows visualization of the real, imaginary, absolute, or both parts of the phasors.
%   - Supports logarithmic or linear scaling for the y-axis.
%   - Handles multi-component phasors (`N×M×(2nh+1)` array) by either combining all components 
%     in one plot or plotting each component separately.
%   - Flexibly integrates with existing figures, axes, or tiled layouts.
%
%   Visualization Details:
%   - If `explosed` is true (default), each matrix component gets its own subplot in a tiled layout.
%   - If `explosed` is false, all components are combined into a single plot.
%   - Phasors can be displayed for all harmonics (`both`) or only non-negative harmonics (`oneSided`).
%
%   Tiled Layout Management:
%   - Automatically manages tiled layouts within existing figures, axes, or parent objects.
%   - Creates new layouts if required, or uses existing layouts if their size matches the data.
%   - Tags the created layouts for reuse and better organization.
%
%   Example Usage:
%   % Plot the absolute values of phasors with a logarithmic scale
%   T = stemPhasor(Aph, 'scale', 'log', 'display', 'abs');
%
%   % Plot real parts of phasors in a single combined plot
%   T = stemPhasor(Aph, 'explosed', false, 'display', 'real');
%
%   % Plot using an existing axes or figure
%   ax = gca;
%   T = stemPhasor(Aph, 'parent', ax, 'side', 'both');
%
%   See also: TILEDLAYOUT, STEM, PLOT.
arguments
    Aph
    varopt.scale {mustBeMember(varopt.scale,{'log','linear'})}='log'
    varopt.hold=ishold
    varopt.explosed= true
    varopt.display {mustBeMember(varopt.display,{'real','imag','both','abs','absangle'})} = 'abs'
    varopt.marker = "o"
    varopt.side {mustBeMember(varopt.side,{'both','oneSided'})} = 'oneSided'
    varopt.parent = gcf
    varopt.uniformYLim = false;
end

if ishold
    varopt.hold=true;
end

if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end

nx=size(Aph,1);
ny=size(Aph,2);
nh=(size(Aph,3)-1)/2;
if nh < 2
    Aph = cat(3,zeros(nx,ny,2-nh),Aph,zeros(nx,ny,2-nh));
    nh = 2;
end

switch varopt.side
    case 'oneSided'
        axe_x = 0:nh;
        phas_index = (nh+1):(2*nh+1);
    otherwise
        axe_x = -nh:nh;
        phas_index = 1:(2*nh+1);
end

%is parent a figure : 
%   - YES : is it empty ? 
%      - YES : create a new tiled layout
%      - NO : is it a tiled layout ?
%          - YES : is it the right size ?
%              - YES : use it
%              - NO : delete it and create a new one
%          - NO : create a new tiled layout
%   - NO : is it a tiled layout ?
%       - YES : is it tagged 'stemPhasor' ?
%           - YES : is it the right size ?
%               - YES : use it
%               - NO : create a new tiled layout inside the parent of parent and tag it 'stemPhasor'
%           - NO : create a new tiled layout inside the parent in a available space and tag it 'stemPhasor'
%       - NO : is it a graphical object able to have a tiledLayout as child ?
%           - YES : create a new tiled layout inside the parent and tag it 'stemPhasor'
%           - NO : create a new figure and a new tiled layout inside it and tag it 'stemPhasor', display a warning message



parent = varopt.parent;
ff = ancestor(parent, 'figure'); % Get the figure handle of the parent
%set(ff, 'Visible', 'off'); % Make the current figure invisible

%if ~ishold % Check if the hold state is off
    %%clf; % Clear the current figure
    %T = createTiledLayout(nx, ny); % Create a new tiled layout
%else

if varopt.explosed
    T = manageTiledLayout2(parent, nx, ny);
else
    T = manageTiledLayout2(parent, 1, 1);
end

if varopt.explosed
    for nyi=1:ny
        for nxi=1:nx
            toto=squeeze(Aph(nxi,nyi,phas_index));
            nexttile(T,sub2ind([ny, nx], nyi, nxi))
            if varopt.hold
                hold on
            end
            switch varopt.display
                case 'real'
                    stem(phas_index-nh-1,squeeze(real(toto)),varopt.marker)
                    set(gca,"yscale",varopt.scale)
                case 'imag'
                    stem(phas_index-nh-1,squeeze(imag(toto)),varopt.marker)
                    set(gca,"yscale",varopt.scale)
                case 'both'
                    stem(phas_index-nh-1,squeeze(real(toto)),varopt.marker)
                    yyaxis right
                    stem(squeeze(imag(toto)),varopt.marker)
                case 'abs'
                    stem(phas_index-nh-1,squeeze(abs(toto)),varopt.marker)
                    set(gca,"yscale",varopt.scale)
                    grid on
                    % yyaxis right
                    % stem(squeeze(angle(toto)),varopt.marker)
                    % ylim([-pi-0.1 pi+0.1])
                    % yyaxis left
                case 'absangle'
                    stem(phas_index-nh-1,squeeze(abs(toto)),varopt.marker)
                    set(gca,"yscale",varopt.scale)
                    grid on
                    ylabel('Abs')
                    yyaxis right
                    stem(phas_index-nh-1,squeeze(angle(toto)),varopt.marker)
                    ylim([-pi-0.1 pi+0.1])
                    ylabel('angle')
                    yyaxis left
            end
            xlim(gca,"padded")
            toto = xlim;
            xticks(gca,(ceil(toto(1)):1:floor(toto(2))))
            grid minor
            if varopt.hold
                hold off
            end
            if nxi == nx
                xlabel('Harmonic order')
            end
        end
    end
        try
            if varopt.uniformYLim
                allAxes = findall(T, 'Type', 'axes');
                % Get the current y-limits of all axes
                yLimits = arrayfun(@(ax) ax.YLim, allAxes, 'UniformOutput', false);
                % Find the global min and max y-limits
                minY = min(cellfun(@(yl) yl(1), yLimits));
                maxY = max(cellfun(@(yl) yl(2), yLimits));
                % Set the same y-limits for all axes
                arrayfun(@(ax) set(ax, 'YLim', [minY, maxY]), allAxes);
            end

            linkaxes(findall(gcf, 'type', 'axes'), 'x')
            ylim auto
        catch
        end
else
    toto=reshape(Aph,[nx*ny 2*nh+1]);
    toto=toto(:,phas_index);
    nexttile(T,1)
    stem(phas_index-nh-1,abs(toto)',varopt.marker)
    xlim(gca,"padded")
    toto = xlim;
    xticks(gca,(ceil(toto(1)):1:floor(toto(2))))
    grid on
    grid minor
    set(gca,"yscale",varopt.scale)
end
if varopt.explosed
    switch varopt.display
        case 'real'
            sgtitle('stem real part of phasor of Matrix')
        case 'imag'
            sgtitle('stem imag part of phasor of Matrix')
        case 'both'
            sgtitle('stem of phasor of Matrix')
        case 'abs'
            sgtitle('stem abs of phasor of Matrix')
    end
else
    switch varopt.display
        case 'real'
            title('stem real part of phasor of Matrix')
        case 'imag'
            title('stem imag part of phasor of Matrix')
        case 'both'
            title('stem of phasor of Matrix')
        case 'abs'
            title('stem abs of phasor of Matrix')
    end
end

set(gcf, 'Visible', 'on'); % Make the current figure invisible
end


function T = manageTiledLayout2(parent, nx, ny, Tag)
%MANAGETILEDLAYOUT2 Create or reuse a tiled layout for `stemPhasor` plots.
%
%   This helper function manages tiled layouts by either creating a new layout or reusing
%   an existing one if it matches the required dimensions and tag.
%
%   Key Features:
%   - Reuses existing layouts when possible for efficiency and consistency.
%   - Automatically adapts to parent objects (figures, axes, or tiled layouts).
%   - Supports tagging for better organization and reuse of layouts.
%
%   Input Arguments:
%   - parent (graphics handle): Parent object where the layout will be created or reused.
%   - nx (integer): Number of rows in the layout.
%   - ny (integer): Number of columns in the layout.
%   - Tag (char, optional): Tag for the layout. Default: `'stemPhasor'`.
%
%   Output Arguments:
%   - T (tiledlayout object): Handle to the created or reused tiled layout.
%
%   Example Usage:
%   % Create a new tiled layout in a figure
%   T = manageTiledLayout2(gcf, 2, 3, 'stemPhasor');
%
%   % Reuse an existing tiled layout
%   T = manageTiledLayout2(ax, 2, 2);
%
%   See also: TILEDLAYOUT.
    arguments
        parent 
        nx {mustBeInteger}
        ny {mustBeInteger}
        Tag string = 'stemPhasor'
    end
    
    % Check if parent is a figure
    if isa(parent, 'matlab.ui.Figure')
        boolVec(1) = true; % 1
        % Check if the figure is empty
        if isempty(parent.Children)
            % Create a new tiled layout
            T = tiledlayout(parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
            T.Tag = Tag;            
            boolVec(2) = true; % 1 1
        else
            boolVec(2) = false; % 1 0
            % Check if the first child is a tiled layout
            if isa(parent.Children(1), 'matlab.graphics.layout.TiledChartLayout')
                boolVec(3) = true; % 1 0 1
                T = parent.Children(1);
                % Check if the tiled layout is the right size
                if all(T.GridSize == [nx, ny])
                    boolVec(4) = true; % 1 0 1 1
                    % Use the existing tiled layout
                    T.Tag = Tag;
                    boolVec;
                    return;
                else
                    boolVec(4) = false; % 1 0 1 0
                    % Delete the existing tiled layout and create a new one
                    delete(T);
                    T = tiledlayout(parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                    T.Tag = Tag;
                end
            else
                boolVec(3) = false; % 1 0 0
                % Create a new tiled layout
                T = tiledlayout(parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                T.Tag = Tag;
            end
        end
    else
        boolVec(1) = false; % 0
        % Check if parent is a tiled layout
        if isa(parent, 'matlab.graphics.layout.TiledChartLayout')
            boolVec(2) = true; % 0 1
            % Check if the tiled layout is tagged Tag
            if strcmp(parent.Tag, Tag) %it is already a "Tag" tiledLayout
                boolVec(3) = true;   % 0 1 1
                % Check if the tiled layout is the right size
                if all(parent.GridSize == [nx, ny])
                    boolVec(4) = true; % 0 1 1 1
                    % Use the existing tiled layout
                    T = parent;
                    T.Tag = Tag;
                else
                    boolVec(4) = false; % 0 1 1 0
                    %get current layout of parent
                    Layout = parent.Layout;
                    % Create a new tiled layout inside the parent of parent and tag it Tag
                    T = tiledlayout(parent.Parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                    T.Tag = Tag;
                    T.Layout = Layout;
                    % Delete the existing tiled layout
                    delete(parent);
                end
            else %parent is a genereic tiledLayout, "Tag" will be a child
                % Check the number of children and if creating a new child would exceed the number of tiles
                boolVec(3) = false; % 0 1 0
                if numel(parent.Children) < prod(parent.GridSize)
                    boolVec(4) = true; % 0 1 0 1
                    % Create a new tiled layout inside the parent in an available space and tag it Tag
                    T = tiledlayout(parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                    T.Tag = Tag;
                    % Move the newly created tiled layout to an available tile
                    T.Layout.Tile = numel(parent.Children);
                else
                    boolVec(4) = false; % 0 1 0 0
                    error('No available tiles in the current tiled layout.');
                end
            end
        else
            boolVec(2) = false;;
            % Check if parent is a graphical object able to have a tiled layout as child
            if isa(parent, 'matlab.graphics.axis.Axes') || isa(parent, 'matlab.graphics.chart.Chart')
                boolVec(3) = true;
                % Create a new tiled layout inside the parent and tag it Tag
                T = tiledlayout(parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                T.Tag = Tag;
            else
                boolVec(3) = false;
                % Create a new figure and a new tiled layout inside it and tag it Tag
                fig = figure;
                T = tiledlayout(fig, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                T.Tag = Tag;
                warning('Parent is not a graphical object able to have a tiled layout as child. Created a new figure.');
            end
        end
    end
    boolVec;
end