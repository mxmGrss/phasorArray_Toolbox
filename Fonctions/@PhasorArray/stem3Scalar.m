function ax = stem3Scalar(o1, nvp)
    % STEM3SCALAR Generate a 3D stem plot for a scalar `PhasorArray` object.
    %
    %   STEM3SCALAR(o1, marker) generates a 3D stem plot to visualize the phasors of a scalar
    %   `PhasorArray` object. The function uses `quiver3` to plot the phasors.
    %
    %   Inputs:
    %     o1     - (PhasorArray) The scalar `PhasorArray` object to be plotted.
    %     marker - (char) Marker style for the plot.
    %
    %   Example Usage:
    %     % Plot a scalar PhasorArray object
    %     A = PhasorArray(rand(1,1,11));
    %     stem3Scalar(A, 'o');
    %
    %   See also: quiver3.
    arguments
        o1
        nvp.marker (1,:) char = 'o'
        nvp.parent = gca
        nvp.scale {mustBeMember(nvp.scale,{'log','linear'})}='log'
        nvp.side {mustBeMember(nvp.side,{'both','oneSided'})}='oneSided'
    end


    if strcmp(nvp.scale,'log')
        oo1 = o1.value;
        oo1a = abs(oo1);
        oo1p = angle(oo1);
        oo1a = log10(oo1a);
        oo1a(oo1a == -Inf) = NaN;
        base = min(floor(oo1a(:)));
        oo1a = oo1a - base;
        newo1 = oo1a.*exp(1i*oo1p);
        o1 = PhasorArray(newo1);
    else
        base = 0;
    end


    h = o1.h;
    x = -h:h;
    y = squeeze(real(o1.value))';
    z = squeeze(imag(o1.value))';

    u = zeros(size(x));
    v = -y;
    w = -z;

    ax = quiver3(nvp.parent,x, y, z, u, v, w,'off', '-<','marker',nvp.marker);
    xlabel('Harmonics');
    ylabel('Real Part');
    zlabel('Imaginary Part');
    title('3D Stem Plot of Scalar PhasorArray');

    % Find all quiver objects in the current axes
    a = gca;
    quiverObjects = findall(a.Children, 'Type', 'quiver');
    clrInd = a.ColorOrderIndex;

    limvar = ceil(max(abs([quiverObjects.YData quiverObjects.ZData]),[],"all"));
    ylim([-limvar limvar]);
    zlim([-limvar limvar]);

    hmax = max([quiverObjects.XData],[],"all");


    varhold = ishold;
    if strcmp(nvp.scale,'log')
        delete(findall(a.Children, 'Type', 'line'));
        title(sprintf('3D Stem Plot of Logarithmic Scalar PhasorArray, with base 10^{%d} and max', base));
        for r=1:ceil(limvar)
            hold on
            for d = 0:pi/(20*r):2*pi
                [tks_x,tks_y] = pol2cart([d d],[r r+0.03]);
                cero = [0 0];
                plot3(cero,tks_x,tks_y,'-','Color','#606060','LineWidth',1)
                %plot3(cero-hmax,tks_x,tks_y,'-','Color','#606060','LineWidth',1)
                %plot3(cero+hmax,tks_x,tks_y,'-','Color','#606060','LineWidth',1)
            end
            hold off
        end
        a.ColorOrderIndex = clrInd;
        grid on
        a.YGrid = 'off';
        a.ZGrid = 'off';
        a.XGrid = 'on';
    else
        title('3D Stem Plot of Scalar PhasorArray');
    end
    if varhold
        hold on
    end

    switch nvp.side
        case 'both'
            xticks(-hmax:hmax);
            xlim([-hmax hmax]);
        case 'oneSided'
            xticks(0:hmax);
            xlim([0 hmax]);
    end

    %bring figure to frond
    %figure(gcf)

end

