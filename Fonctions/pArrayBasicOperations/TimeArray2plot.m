function t = TimeArray2plot(Mt, t, T, nvp)
%TIMEARRAY2PLOT Plot time-domain data produced by PhasorArray2time.
%
%   t = TimeArray2plot(Mt, t, T, nvp)
%
%   Extracted from PhasorArray2time to isolate the plotting logic; the visual
%   result is the one PhasorArray2time produced before the split.
%
%   Mt  - Time domain data (N x M x P)
%   t   - Time vector (returned, since the layout may adjust it)
%   T   - Period
%   nvp - Options struct forwarded from PhasorArray2time.

% Ensure t is computed if empty
s = size(Mt);
if isempty(t)
    t = linspace(0, T, s(3));
end

nx = s(1);
ny = s(2);

% Extract options with defaults
explosed = true; if isfield(nvp, 'explosed'), explosed = nvp.explosed; end
is_hold = false; if isfield(nvp, 'hold'), is_hold = nvp.hold; end
plot3D = false;  if isfield(nvp, 'plot3D'), plot3D = nvp.plot3D; end
title_str = '';  if isfield(nvp, 'title'), title_str = nvp.title; end
parent = [];     if isfield(nvp, 'parent'), parent = nvp.parent; end
zero_centered = false; if isfield(nvp, 'ZeroCentered'), zero_centered = nvp.ZeroCentered; end

disp_real = true; if isfield(nvp, 'DispReal'), disp_real = nvp.DispReal; end
disp_imag = false; if isfield(nvp, 'DispImag'), disp_imag = nvp.DispImag; end
grid_opt = 'on'; if isfield(nvp, 'grid'), grid_opt = nvp.grid; end
line_style = '-'; if isfield(nvp, 'LineStyle'), line_style = nvp.LineStyle; end
global_ylim = false; if isfield(nvp, 'GlobalYLim'), global_ylim = nvp.GlobalYLim; end
link_axes = 'off'; if isfield(nvp, 'linkaxes'), link_axes = nvp.linkaxes; end
xlabelStr = 'time (sec)'; if isfield(nvp, 'xlabelStr'), xlabelStr = nvp.xlabelStr; end

singlePart = xor(disp_imag, disp_real) || plot3D;

if explosed
    if singlePart
        TL = manageTiledLayout(parent, nx, ny, "plotTimePhasor", "ishold", is_hold);
        for nxi = 1:nx
            for nyi = 1:ny
                % Explicit tile index: on a held layout the tiles are already
                % filled, and a bare nexttile would run past the grid.
                ax = nexttile(TL, sub2ind([ny, nx], nyi, nxi));
                if is_hold, hold(ax, 'on'); end

                y = squeeze(Mt(nxi,nyi,:));

                if plot3D
                    plot3(ax, real(y), imag(y), t, line_style);
                    xlabel(ax, "Re(a_{"+num2str(nxi)+num2str(nyi)+"})");
                    ylabel(ax, "Im(a_{"+num2str(nxi)+num2str(nyi)+"})");
                    zlabel(ax, xlabelStr);
                    ylim(ax, 'auto'); xlim(ax, 'auto');
                    if zero_centered
                        ylim(ax, max(abs(ylim(ax))).*[-1 1]);
                        xlim(ax, max(abs(xlim(ax))).*[-1 1]);
                    end
                else
                    if disp_imag
                        plot(ax, t, imag(y), line_style);
                    else
                        plot(ax, t, real(y), line_style);
                    end
                    if nxi == nx
                        xlabel(ax, xlabelStr);
                    end
                    ylim(ax, 'auto');
                    if zero_centered
                        ylim(ax, max(abs(ylim(ax))).*[-1 1]);
                    end
                end

                grid(ax, 'off');
                grid(ax, grid_opt);
            end
        end
    else
        TL = manageTiledLayout(parent, nx, 2*ny, "plotTimePhasor", "ishold", is_hold);
        for nxi = 1:nx
            for nyi = 1:ny
                % Real part
                ax1 = nexttile(TL, sub2ind([2*ny, nx], (nyi-1)*2+1, nxi));
                if is_hold, hold(ax1, 'on'); end
                plot(ax1, t, real(squeeze(Mt(nxi,nyi,:))), line_style);
                ylabel(ax1, "Re(a_{"+num2str(nxi)+num2str(nyi)+"})");
                ylim(ax1, 'auto');
                if zero_centered
                    ylim(ax1, max(abs(ylim(ax1))).*[-1 1]);
                end
                % Both columns of the bottom row carry the abscissa label; the
                % pre-split code labelled only the imaginary one.
                if nxi == nx
                    xlabel(ax1, xlabelStr);
                end
                grid(ax1, 'off'); grid(ax1, grid_opt);

                % Imaginary part
                ax2 = nexttile(TL, sub2ind([2*ny, nx], nyi*2, nxi));
                if is_hold, hold(ax2, 'on'); end
                plot(ax2, t, imag(squeeze(Mt(nxi,nyi,:))), line_style);
                ylabel(ax2, "Im(a_{"+num2str(nxi)+num2str(nyi)+"})");
                ylim(ax2, 'auto');
                if zero_centered
                    ylim(ax2, max(abs(ylim(ax2))).*[-1 1]);
                end
                if nxi == nx
                    xlabel(ax2, xlabelStr);
                end
                grid(ax2, 'off'); grid(ax2, grid_opt);
            end
        end
    end

    if ~isempty(title_str)
        sgtitle(TL, title_str);
    end

    ax = findall(TL, 'Type', 'axes');

    yset = {'y','xy','yx','yz','zy','xyz','yxz','xzy','zxy','zyx','yzx'};
    xset = {'x','xy','yx','xz','zx','xyz','yxz','xzy','zxy','zyx','yzx'};

    if global_ylim || any(strcmp(link_axes, yset))
        uu = max(abs(cell2mat(ylim(ax))), [], 'all');
        set(ax, 'ylim', uu*[-1, 1]);
    end

    if plot3D && (global_ylim || any(strcmp(link_axes, xset)))
        uu = max(abs(cell2mat(xlim(ax))), [], 'all');
        set(ax, 'xlim', uu*[-1, 1]);
        Link = linkprop(ax(:), {'CameraUpVector','CameraPosition','CameraTarget','ZLim'});
        setappdata(gcf, 'StoreTheLink', Link);
    end

    linkaxes(ax, link_axes);
else
    % Not exploded: every coefficient superposed on one pair of axes.
    if singlePart
        TL = manageTiledLayout(parent, 1, 1, "plotTimePhasor", "ishold", is_hold);
        ax = nexttile(TL, 1);
        if is_hold, hold(ax, 'on'); end

        y = reshape(Mt, [], numel(t));

        if plot3D
            plot3(ax, real(y).', imag(y).', repmat(t(:), 1, size(y,1)), line_style);
            xlabel(ax, 'Re(M(t))');
            ylabel(ax, 'Im(M(t))');
            zlabel(ax, xlabelStr);
        elseif disp_imag
            plot(ax, t, imag(y), line_style);
            if isempty(title_str), title(ax, 'M(t), imaginary part'); else, title(ax, title_str); end
        else
            plot(ax, t, real(y), line_style);
            if isempty(title_str), title(ax, 'M(t), real part'); else, title(ax, title_str); end
        end

        ylim(ax, 'auto');
        if zero_centered
            ylim(ax, max(abs(ylim(ax))).*[-1 1]);
        end
        grid(ax, 'off'); grid(ax, grid_opt);
    else
        TL = manageTiledLayout(parent, 1, 2, "plotTimePhasor", "ishold", is_hold);
        y = reshape(Mt, [], numel(t));

        ax1 = nexttile(TL, 1);
        if is_hold, hold(ax1, 'on'); end
        plot(ax1, t, real(y), line_style);
        ylim(ax1, 'auto');
        if zero_centered
            ylim(ax1, max(abs(ylim(ax1))).*[-1 1]);
        end
        title(ax1, 'M(t), real part');
        grid(ax1, 'off'); grid(ax1, grid_opt);

        ax2 = nexttile(TL, 2);
        if is_hold, hold(ax2, 'on'); end
        plot(ax2, t, imag(y), line_style);
        ylim(ax2, 'auto');
        if zero_centered
            ylim(ax2, max(abs(ylim(ax2))).*[-1 1]);
        end
        title(ax2, 'M(t), imag part');
        grid(ax2, 'off'); grid(ax2, grid_opt);
    end
end

end
