function  T = stemPhasor(Aph,nvp)
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
    nvp.scale {mustBeMember(nvp.scale,{'log','linear'})}='log'
    nvp.hold=ishold
    nvp.explosed= true
    nvp.display {mustBeMember(nvp.display,{'real','imag','both','abs','absangle'})} = 'abs'
    nvp.marker = "o"
    nvp.side {mustBeMember(nvp.side,{'both','oneSided'})} = 'oneSided'
    nvp.parent = []
    nvp.uniformYLim = false;
end

if ishold
    nvp.hold=true;
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

switch nvp.side
    case 'oneSided'
        axe_x = 0:nh;
        phas_index = (nh+1):(2*nh+1);
    otherwise
        axe_x = -nh:nh;
        phas_index = 1:(2*nh+1);
end





% parent = nvp.parent;
% ff = ancestor(parent, 'figure'); % Get the figure handle of the parent
% set(ff, 'Visible', 'off'); % Make the current figure invisible

%if ~ishold % Check if the hold state is off
%%clf; % Clear the current figure
%T = createTiledLayout(nx, ny); % Create a new tiled layout
%else

if nvp.explosed
    T = manageTiledLayout(nvp.parent, nx, ny);
else
    T = manageTiledLayout(nvp.parent, 1, 1);
end

if nvp.explosed
    for nyi=1:ny
        for nxi=1:nx
            toto=squeeze(Aph(nxi,nyi,phas_index));
            nexttile(T,sub2ind([ny, nx], nyi, nxi))
            if nvp.hold
                hold on
            end
            switch nvp.display
                case 'real'
                    stem(phas_index-nh-1,squeeze(real(toto)),nvp.marker)
                    set(gca,"yscale",nvp.scale)
                case 'imag'
                    stem(phas_index-nh-1,squeeze(imag(toto)),nvp.marker)
                    set(gca,"yscale",nvp.scale)
                case 'both'
                    stem(phas_index-nh-1,squeeze(real(toto)),nvp.marker)
                    yyaxis right
                    stem(squeeze(imag(toto)),nvp.marker)
                case 'abs'
                    stem(phas_index-nh-1,squeeze(abs(toto)),nvp.marker)
                    set(gca,"yscale",nvp.scale)
                    grid on
                    % yyaxis right
                    % stem(squeeze(angle(toto)),nvp.marker)
                    % ylim([-pi-0.1 pi+0.1])
                    % yyaxis left
                case 'absangle'
                    stem(phas_index-nh-1,squeeze(abs(toto)),nvp.marker)
                    set(gca,"yscale",nvp.scale)
                    grid on
                    ylabel('Abs')
                    yyaxis right
                    stem(phas_index-nh-1,squeeze(angle(toto)),nvp.marker)
                    ylim([-pi-0.1 pi+0.1])
                    ylabel('angle')
                    yyaxis left
            end
            xlim(gca,"padded")
            toto = xlim;
            xticks(gca,(ceil(toto(1)):1:floor(toto(2))))
            grid minor
            if nvp.hold
                hold off
            end
            if nxi == nx
                xlabel('Harmonic order')
            end
        end
    end
    try
        if nvp.uniformYLim
            allAxes = findall(T, 'Type', 'axes');
            % Get the current y-limits of all axes
            yLimits = arrayfun(@(ax) ax.YLim, allAxes, 'UniformOutput', false);
            % Find the global min and max y-limits
            minY = min(cellfun(@(yl) yl(1), yLimits));
            maxY = max(cellfun(@(yl) yl(2), yLimits));
            % Set the same y-limits for all axes
            arrayfun(@(ax) set(ax, 'YLim', [minY, maxY]), allAxes);
        end

        linkaxes(T, 'x')
        ylim auto
    catch
    end
else
    toto=reshape(Aph,[nx*ny 2*nh+1]);
    toto=toto(:,phas_index);
    nexttile(T,1)
    if nvp.hold
        hold on
    end
    stem(phas_index-nh-1,abs(toto)',nvp.marker)
    if nvp.hold
        hold off
    end
    xlim(gca,"padded")
    toto = xlim;
    xticks(gca,(ceil(toto(1)):1:floor(toto(2))))
    grid on
    grid minor
    set(gca,"yscale",nvp.scale)
end
if nvp.explosed
    switch nvp.display
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
    switch nvp.display
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

