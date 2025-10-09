function [Mt,t] = PhasorArray2time(Mph,T,t,arg)
% PHASORARRAY2TIME Evaluate a periodic matrix function A(t) or A(θ) from its phasor representation.
%
%   PHASORARRAY2TIME reconstructs the time-domain representation of a periodic 
%   matrix function A(t) given its phasor coefficients. It can evaluate:
%   - A(t) for a fixed period `T`
%   - A(θ) for a vector of angles `θ`
%   - A(t) when `T` and `t` are both vectors, in which case `T` is interpreted as θ(t).
%
%   Inputs:
%     Mph  - (m × n × (2h+1) array) Phasor coefficients of the periodic matrix A(t).
%     T    - (double or vector, optional) Period of the function or a vector of angles `θ`:
%              - If `T` is a **scalar**, it represents the period of A(t).
%              - If `T` is a **vector**:
%                 - If `t` is empty, it is interpreted as `θ`, and A(θ) is computed.
%                 - If `T` and `t` have matching sizes, `T` is interpreted as a time-dependent frequency phase `θ(t)`.
%              - Default: `1`.
%     t    - (vector, optional) Time instants at which to evaluate A(t).
%              - If `T` is a vector and `t` is empty, then `T` is interpreted as `θ`, and `t = T`.
%              - If `T` and `t` have matching sizes, `T` is interpreted as `θ(t)`, defining a time-varying phase.
%              - If empty, a default time grid is computed.
%     arg  - (struct, optional) Name-value pair arguments:
%              - 'plot' (logical): Plot matrix coefficients versus time (default: `false`).
%              - 'explosed' (logical): Plot each coefficient in a separate subplot (default: `true`).
%              - 'hold' (logical): Hold the current plot (default: `false`).
%              - 'DispImag' (logical): Display the imaginary part of the matrix (default: `false`).
%              - 'DispReal' (logical): Display the real part of the matrix (default: `true`).
%              - 'ZeroCentered' (logical): Center the plot around zero (default: `false`).
%              - 'forceReal' (logical): Force the matrix to be real-valued (default: `[]`).
%              - 'checkReal' (logical): Check if the matrix is real-valued (default: `false`).
%              - 'checkRealTol' (double): Tolerance for real-valued check (default: `1e-8`).
%              - 'title' (string): Title of the plot (default: `[]`).
%              - 'plot3D' (logical): Plot the matrix coefficients in 3D (default: `false`).
%              - 'linetype' (string): Line type for the plot (default: `'-'`).
%              - 'GlobalYLim' (logical): Set uniform y-limits for all subplots (default: `false`).
%              - 'linkaxes' (string): Link axes ('x', 'y', 'xy', etc.) (default: `'x'`).
%              - 'computationMethod' (integer): Method used for evaluation (default: `1`).
%              - 'providedPhasorForm' (string): Phasor representation ('exp' or 'SinCos') (default: `"exp"`).
%              - 'parent' (handle): Parent figure or axes (default: `[]`).
%
%   Outputs:
%     Mt   - (m × n × length(t) array) Evaluated matrix A(t) or A(θ).
%     t    - (vector) Time instants (or angles θ) at which A(t) is evaluated.
%
%   Behavior:
%     - **Evaluation Modes**:
%         - If `T` is a scalar and `t` is provided, computes **A(t) for a fixed period T**.
%         - If `T` is a vector and `t` is empty, computes **A(θ)** using `θ = T`.
%         - If `T` and `t` are both vectors of the same length, computes **A(t) with θ(t) as a time-varying frequency phase**.
%     - **Automatic Time Grid**:
%         - If `t` is empty, it is computed automatically based on `T` and the highest harmonic.
%     - **Real-Valued Constraints**:
%         - Allows optional enforcement of real-valued output.
%     - **Plotting Options**:
%         - Can generate time-domain plots of individual matrix elements.
%
%   Example Usage:
%     % Evaluate a periodic matrix at given time points
%     Mph = rand(3, 3, 5);
%     T = 1;
%     t = 0:0.01:1;
%     [Mt, t] = PhasorArray2time(Mph, T, t);
%
%     % Evaluate at a vector of angles instead of time
%     theta = linspace(0, 2*pi, 100);
%     [Mt, theta] = PhasorArray2time(Mph, theta);
%
%     % Evaluate A(t) with a time-varying phase θ(t)
%     t = linspace(0, 1, 100);
%     theta_t = 2*pi*t.^2;  % Example of a quadratic phase modulation
%     [Mt, t] = PhasorArray2time(Mph, theta_t, t);
%
%     % Generate a time-domain plot
%     [Mt, t] = PhasorArray2time(Mph, T, t, 'plot', true);
%
%   See also: timeArray2Phasors, stemPhasor.

arguments
    Mph
    T=1
    t=[]
    arg.plot logical         = false
    arg.explosed logical     = true
    arg.hold logical         = false
    arg.DispImag logical     = false
    arg.DispReal logical     = true
    arg.ZeroCentered logical = false
    arg.forceReal logical    = []
    arg.checkReal logical    = false
    arg.checkRealTol         = 1e-8
    arg.title                = []
    arg.plot3D logical       = false
    arg.linetype             = '-'
    arg.GlobalYLim logical   = false
    arg.linkaxes             = 'x';
    arg.computationMethod    = 1
    arg.providedPhasorForm {mustBeMember(arg.providedPhasorForm,["exp","SinCos"])} = "exp"
    arg.parent               = [];
end

% Use a less efficient but compatible computation method for older MATLAB releases
if isMATLABReleaseOlderThan("R2022a")
    arg.computationMethod = 3;
end

if isempty(arg.forceReal)
    if arg.providedPhasorForm == "exp"
        if isrealp(Mph)
            arg.forceReal = true;
        end
    else
        arg.forceReal = false;
    end
end

% User specified that the input PhasorArray is in SinCos form
if arg.providedPhasorForm == "SinCos"
    % Ensure that the matrix is real valued
    if ~arg.forceReal
        arg.forceReal = true;
        warning("SinCos phasor form is only compatible with real valued matrices. Forcing real valued computation.")
    end

    % Convert the matrix to SinCos form if it is not already
    if isa(Mph, 'PhasorArray')
        Mph = Mph.SinCosForm();
        warning("Converting PhasorArray to SinCos form.")
    end


end

% Convert PhasorArray or symbolic variables to numeric values
Mph = convertToNumeric(Mph);

% Get the size of Mph
nx = size(Mph,1);
ny = size(Mph,2);
h = size(Mph,3);
h = (h - 1) / 2;

% Compute the time vector if not provided
if and(numel(t) < 3 , numel(T) == 1)
    dt = computeTimeStep(T, h);
    t = computeTimeVector(t, dt, T);
end

% Ensure T and t are row vectors
T = reshape(T, 1, []);
t = reshape(t, 1, []);

% Compute theta
theta = computeTheta(T, t);

% Compute the basis a each time
% case of complex valued matrix
if ~arg.forceReal
    eit = exp(1i*(-h:h)'*theta);
    Meval=Mph;
    %case of real valued matrix
else
    %convert to real valued matrix in sin/cos form
    if arg.providedPhasorForm == "exp"
        Mphr  = real(Mph + flip(Mph,3))/2 + 1i * imag(Mph - flip(Mph,3))/2;
        Meval = real(cat(3,1i*(flip(Mphr(:,:,h+2:end),3)-Mphr(:,:,1:h)),Mphr(:,:,h+1),(Mphr(:,:,h+2:end)+flip(Mphr(:,:,1:h),3))));
    else
        Meval = Mph;
    end
    eit = [sin((h:-1:1)'*theta);cos((0:h)'*theta)];

end

if isa(t,'sym')
    arg.computationMethod=3;
end

%choose the computation method to compute the matrix at each time
switch arg.computationMethod
    case 1
        %best time but tensorprod reliant so matlab >22a
        Mt=tensorprod(Meval,double(eit),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))

    case 2
        %bad time very slow
        reM=reshape(Meval,nx,[],1);
        reEit=kron(eit,eye(ny));
        rMt=tensorprod(reM,reEit,2,1);
        Mt=reshape(rMt,nx,ny,[]);

    otherwise
        %a bit better thanks to sparse
        reM=reshape(Meval,nx,[],1);
        reEit=kron(eit,speye(ny));
        rbMt=reM*reEit;
        Mt=reshape(rbMt,nx,ny,[]);
end

% Check if the imaginary part of Mt is negligible
if arg.checkReal
    if sum(abs(imag(Mt))>arg.checkRealTol)==0
        Mt=real(Mt);
    else
        % if imag part is not negligeable, output is complex valued and a warning is displayed
        warning('Argument enforce reality of M but imag part is not negligeable, output is complex valued')
        return
    end
end

% Mt=real(Mt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Ploting function  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if arg.plot
    t = plot2(arg,t,T,Mt);

end


    function Mph = convertToNumeric(Mph)
        if isa(Mph, 'PhasorArray')
            Mph = Mph.Value;
        elseif isa(Mph, "ndsdpvar") || isa(Mph, 'sdpvar')
            Mph = value(Mph);
        elseif isa(Mph, 'sym')
            Mph = vpa(value(Mph));
            arg.computationMethod = 3;
        end
    end
    function dt = computeTimeStep(T, h)
        % This function computes the time step based on the period T and the harmonic order h
        % Inputs: T - period, h - harmonic order
        % Output: dt - time step

        dt = T / 2 ^ (nextpow2(max((h * 8), 64)));
    end

    function t = computeTimeVector(t, dt, T)
        % This function computes the time vector based on the input t, the time step dt, and the period T
        % Inputs: t - input time, dt - time step, T - period
        % Output: t - time vector

        switch numel(t)
            case 0 %default one period evaluation
                t = 0:dt:T-dt;
            case 1
                % DO NOTHING scalar case of evaluation at only one given time
            case 2
                t = t(1):dt:t(2);
        end
    end

    function theta = computeTheta(T, t)
        % This function computes the theta value based on the period T and the time vector t
        % Inputs: T - period, t - time vector
        % Output: theta - computed theta value

        if numel(T) > 1
            theta = T;
        else
            theta = 2 * pi / T * t;
        end
    end
end


function T = manageTiledLayout2(parent, nx, ny, Tag)
    arguments
        parent 
        nx 
        ny 
        Tag string = 'plotTimePhasor'
    end
    % Check if parent is a figure
    if isa(parent, 'matlab.ui.Figure')
        % Check if the figure is empty
        if isempty(parent.Children)
            % Create a new tiled layout
            T = tiledlayout(parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
            T.Tag = Tag;            
        else
            % Check if the first child is a tiled layout
            if isa(parent.Children(1), 'matlab.graphics.layout.TiledChartLayout')
                T = parent.Children(1);
                % Check if the tiled layout is the right size
                if all(T.GridSize == [nx, ny])
                    % Use the existing tiled layout
                    T.Tag = Tag;
                    return;
                else
                    % Delete the existing tiled layout and create a new one
                    delete(T);
                    T = tiledlayout(parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                    T.Tag = Tag;
                end
            else
                % Create a new tiled layout
                T = tiledlayout(parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                T.Tag = Tag;
            end
        end
    else
        % Check if parent is a tiled layout
        if isa(parent, 'matlab.graphics.layout.TiledChartLayout')
            % Check if the tiled layout is tagged Tag
            if strcmp(parent.Tag, Tag) %it is already a "Tag" tiledLayout
                % Check if the tiled layout is the right size
                if all(parent.GridSize == [nx, ny])
                    % Use the existing tiled layout
                    T = parent;
                    T.Tag = Tag;
                else
                    % Create a new tiled layout inside the parent of parent and tag it Tag
                    T = tiledlayout(parent.Parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                    T.Tag = Tag;
                end
            else %parent is a genereic tiledLayout, "Tag" will be a child
                % Check the number of children and if creating a new child would exceed the number of tiles
                try
                    % Create a new tiled layout inside the parent in an available space and tag it Tag
                    T = tiledlayout(parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                    T.Tag = Tag;
                    % Move the newly created tiled layout to an available tile
                    T.Layout.Tile = numel(parent.Children);
                catch e
                    error('No available tiles in the current tiled layout.');
                end
            end
        else
            % Check if parent is a graphical object able to have a tiled layout as child
            if isa(parent, 'matlab.graphics.axis.Axes') || isa(parent, 'matlab.graphics.chart.Chart')
                % Create a new tiled layout inside the parent and tag it Tag
                T = tiledlayout(parent, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                T.Tag = Tag;
            else
                % Create a new figure and a new tiled layout inside it and tag it Tag
                fig = figure;
                T = tiledlayout(fig, nx, ny, 'TileSpacing', 'compact', 'Padding', 'compact');
                T.Tag = Tag;
                warning('Parent is not a graphical object able to have a tiled layout as child. Created a new figure.');
            end
        end
    end
end

function t = plot1(arg,t,T,Mt)


nx = size(Mt,1);
ny = size(Mt,2);

if isempty(arg.parent) 
    % If the parent figure is not specified, use the current figure
    parent = gcf;
else
    % If the parent figure is specified, use it
    parent = arg.parent;
end



if ~arg.hold && arg.explosed
    clf
end

% If T is a vector, t must be either empty or have the same number of elements as T for plotting
if isempty(t) % case of a phase vector over T
    t=T;
    xlabelStr= 'angle (rad)';
elseif numel(T)>1 && numel(t) ~= numel(T)
    error('If T is a vector, t must be either empty or have the same number of elements as T for plotting')
else
    xlabelStr= 'time (sec)';
end


% % Calculate the number of expected subplots based on the display options
% if xor(arg.DispImag, arg.DispReal) || arg.plot3D
%     % If only one part is displayed or 3D plot is enabled, one subplot per matrix element
%     numSubplots = nx * ny;
%     doublesubplot = 0;
% else
%     % If both real and imaginary parts are displayed, two subplots per matrix element
%     numSubplots = 2 * nx * ny;
%     doublesubplot = 1;
% end
% if arg.explosed
%     T = manageTiledLayout2(parent,numSubplots,nx,ny*(1+doublesubplot));
% else
%     T = manageTiledLayout2(parent,numSubplots,1,1+doublesubplot);
% end

if arg.explosed
    % Find all axes in the current figure to prevent bugs when plotting on top of an old plot
    old_ax=findall(gcf,'Type','axes');
    if ~isempty(old_ax)
        % Unlink all axes to prevent bugs when plotting on top of an old plot
        linkaxes(old_ax,'');
    end
    % numOldA stores the number of old axes in the figure
    numOldA=numel(old_ax);

    % arg.explosed checks if the plot should be exploded or not
    nx=size(Mt,1);
    ny=size(Mt,2);
    if xor(arg.DispImag,arg.DispReal) || arg.plot3D
        % if only one part is displayed, only one subplot is needed for each coefficient
        ax = gobjects(nx,ny);
        expectedPlotNumber=nx*ny;
        if numOldA==expectedPlotNumber
            old_ax=reshape(flip(old_ax),nx,ny);
        end
    else
        % if both parts are displayed, two subplots are needed for each coefficient
        ax = gobjects(2*nx,ny);
        expectedPlotNumber=nx*ny*2;
        if numOldA==expectedPlotNumber
            old_ax=reshape(flip(old_ax),nx*2,ny);
        end
    end



    for nxi=1:nx
        for nyi=1:ny
            %   if only one part is displayed, only one subplot is needed for each coefficient
            if xor(arg.DispImag,arg.DispReal) || arg.plot3D
                ax(nxi,nyi)=subplot(nx,ny,(nxi-1)*ny+nyi);
                if arg.hold
                    hold on
                end
                if arg.plot3D
                    plot3(real(squeeze(Mt(nxi,nyi,:))),imag(squeeze(Mt(nxi,nyi,:))),t,arg.linetype)
                    xlabel("Re(a_{"+num2str(nxi)+num2str(nyi)+"})")
                    ylabel("Im(a_{"+num2str(nxi)+num2str(nyi)+"})")
                    zlabel(xlabelStr)

                    ylim('auto')
                    xlim('auto')
                    if arg.ZeroCentered
                        ylim(max(abs(ylim)).*[-1 1])
                        xlim(max(abs(xlim)).*[-1 1])
                    end
                else
                    if arg.DispImag
                        plot(t,imag(squeeze(Mt(nxi,nyi,:))),arg.linetype)
                        % ylabel("Im(a_{"+num2str(nxi)+num2str(nyi)+"})")
                    else
                        plot(t,real(squeeze(Mt(nxi,nyi,:))),arg.linetype)
                        % if isreal(Mt)
                        %     ylabel("a_{"+num2str(nxi)+num2str(nyi)+"}")
                        % else
                        %     ylabel("Re(a_{"+num2str(nxi)+num2str(nyi)+"})")
                        % end
                    end
                    if nxi == nx
                        xlabel(xlabelStr)
                    end
                    ylim('auto')
                    if arg.ZeroCentered
                        ylim(max(abs(ylim)).*[-1 1])
                    end
                end
                %  if both parts are displayed, two subplots are needed for each coefficient
            else
                % real part
                ax(2*(nxi-1)*ny+nyi)=subplot(2*nx,ny,2*(nxi-1)*ny+nyi);
                if arg.hold
                    hold on
                end
                plot(t,real(squeeze(Mt(nxi,nyi,:))),arg.linetype)
                ylabel("Re(a_{"+num2str(nxi)+num2str(nyi)+"})")
                ylim('auto')
                if arg.ZeroCentered
                    ylim(max(abs(ylim)).*[-1 1])
                end

                grid off
                grid minor

                % imaginary part
                ax((2*(nxi-1)+1)*ny+nyi)=subplot(2*nx,ny,(2*(nxi-1)+1)*ny+nyi);
                if arg.hold
                    hold on
                end
                plot(t,imag(squeeze(Mt(nxi,nyi,:))),arg.linetype)
                ylabel("Im(a_{"+num2str(nxi)+num2str(nyi)+"})")
                ylim('auto')
                if arg.ZeroCentered
                    ylim(max(abs(ylim)).*[-1 1])
                end

                if nxi == nx
                    xlabel(xlabelStr)
                end

            end
            grid off
            grid minor


        end
    end
    if isempty(arg.title)
        % sgtitle('M(t), vue explosée de la matrice')
    else
        sgtitle(arg.title)
    end
    % Link the axes if the user specified it
    % Possible values for arg.linkaxes are : 'x', 'y', 'z', 'xy', 'yx', 'yz', 'zy', 'xyz', 'yxz', 'xzy', 'zxy', 'zyx', 'yzx'
    uuu_y={'y','xy','yx','yz','zy','xyz','yxz','xzy','zxy','zyx','yzx'};
    uuu_x={'x','xy','yx','xz','zx','xyz','yxz','xzy','zxy','zyx','yzx'};
    uuu_z={'z','zy','yz','xz','zx','xyz','yxz','xzy','zxy','zyx','yzx'};

    % Set the same y limits for all the subplots if the user specified it or if the user wants to link the y axes
    if arg.GlobalYLim || any(strcmp(arg.linkaxes,uuu_y))
        uu=max(abs( cell2mat(ylim(ax))),[],'all');
        set(ax,'ylim',uu*[-1,1])
    end

    % Set the same x limits for all the subplots if the user specified it or if the user wants to link the x axes
    if arg.plot3D
        if arg.GlobalYLim || any(strcmp(arg.linkaxes,uuu_x))
            uu=max(abs( cell2mat(xlim(ax))),[],'all');
            set(ax,'xlim',uu*[-1,1])
            Link = linkprop(ax(:),{'CameraUpVector', 'CameraPosition', ...
                'CameraTarget', 'ZLim'});
            setappdata(gcf, 'StoreTheLink', Link);
        end
    end

    linkaxes(ax,arg.linkaxes);
else

    % If the plot is not exploded, plot the matrix coefficients in a single figure
    if xor(arg.DispImag,arg.DispReal)
        if arg.hold
            hold on
        end
        if arg.DispImag
            plot(t,imag(reshape(Mt,[],numel(t))),arg.linetype)
            if isempty(arg.title)
                title('M(t), imag part')
            else
                title(arg.title)
            end
        else
            plot(t,real(reshape(Mt,[],numel(t))),arg.linetype)
            if isempty(arg.title)
                title('M(t), real part')
            else
                title(arg.title)
            end
        end

        ylim('auto')
        if arg.ZeroCentered
            %arg.zeroCentered checks if the plot should be centered on 0
            ylim(max(abs(ylim)).*[-1 1])
        end

    else
        ff = gcf; 
        if isa(ff.Children(1),"TiledChartLayout")
            TOuter = ff.Children(1);
        else
            TOuter = gcf;
        end

        % If both parts are displayed, plot the real and imaginary parts of the matrix coefficients superposed in two subplots
        nexttile(TOuter,1)
        if arg.hold
            % arg.hold checks if the plot should be superposed on the current figure
            hold on
        end
        plot(t,real(reshape(Mt,[],numel(t))),arg.linetype)

        ylim('auto')
        if arg.ZeroCentered
            ylim(max(abs(ylim)).*[-1 1])
        end
        title('M(t), real part')

        nexttile(TOuter,2)
        if arg.hold
            % arg.hold checks if the plot should be superposed on the current figure
            hold on
        end
        plot(t,imag(reshape(Mt,[],numel(t))),arg.linetype)

        ylim('auto')
        if arg.ZeroCentered
            % make y limits symetric
            ylim(max(abs(ylim)).*[-1 1])
        end
        title('M(t), imag part')
    end
end
end


function t = plot2(arg,t,T,Mt)


nx = size(Mt,1);
ny = size(Mt,2);

if isempty(arg.parent) 
    % If the parent figure is not specified, use the current figure
    parent = gcf;
else
    % If the parent figure is specified, use it
    parent = arg.parent;
end



if ~arg.hold && arg.explosed
    clf
end

% If T is a vector, t must be either empty or have the same number of elements as T for plotting
if isempty(t) % case of a phase vector over T
    t=T;
    xlabelStr= 'angle (rad)';
elseif numel(T)>1 && numel(t) ~= numel(T)
    error('If T is a vector, t must be either empty or have the same number of elements as T for plotting')
else
    xlabelStr= 'time (sec)';
end


% Calculate the number of expected subplots based on the display options
if xor(arg.DispImag, arg.DispReal) || arg.plot3D
    % If only one part is displayed or 3D plot is enabled, one subplot per matrix element
    numSubplots = nx * ny;
    doublesubplot = 0;
else
    % If both real and imaginary parts are displayed, two subplots per matrix element
    numSubplots = 2 * nx * ny;
    doublesubplot = 1;
end
if arg.explosed
    T = manageTiledLayout2(parent,nx,ny*(1+doublesubplot));
else
    T = manageTiledLayout2(parent,1,1+doublesubplot);
end

if arg.explosed
    for nxi=1:nx
        for nyi=1:ny
            %   if only one part is displayed, only one subplot is needed for each coefficient
            if xor(arg.DispImag,arg.DispReal) || arg.plot3D
                nexttile(T,sub2ind([ny, nx], nyi, nxi))
                if arg.hold
                    hold on
                end
                if arg.plot3D
                    plot3(real(squeeze(Mt(nxi,nyi,:))),imag(squeeze(Mt(nxi,nyi,:))),t,arg.linetype)
                    xlabel("Re(a_{"+num2str(nxi)+num2str(nyi)+"})")
                    ylabel("Im(a_{"+num2str(nxi)+num2str(nyi)+"})")
                    zlabel(xlabelStr)

                    ylim('auto')
                    xlim('auto')
                    if arg.ZeroCentered
                        ylim(max(abs(ylim)).*[-1 1])
                        xlim(max(abs(xlim)).*[-1 1])
                    end
                else
                    if arg.DispImag
                        plot(t,imag(squeeze(Mt(nxi,nyi,:))),arg.linetype)
                        % ylabel("Im(a_{"+num2str(nxi)+num2str(nyi)+"})")
                    else
                        plot(t,real(squeeze(Mt(nxi,nyi,:))),arg.linetype)
                        % if isreal(Mt)
                        %     ylabel("a_{"+num2str(nxi)+num2str(nyi)+"}")
                        % else
                        %     ylabel("Re(a_{"+num2str(nxi)+num2str(nyi)+"})")
                        % end
                    end
                    if nxi == nx
                        xlabel(xlabelStr)
                    end
                    ylim('auto')
                    if arg.ZeroCentered
                        ylim(max(abs(ylim)).*[-1 1])
                    end
                end
                %  if both parts are displayed, two subplots are needed for each coefficient
            else
                % real part

                nexttile(T,sub2ind([2*ny, nx], (nyi-1)*2+1, nxi))
                % ax(2*(nxi-1)*ny+nyi)=subplot(2*nx,ny,2*(nxi-1)*ny+nyi);
                if arg.hold
                    hold on
                end
                plot(t,real(squeeze(Mt(nxi,nyi,:))),arg.linetype)
                ylabel("Re(a_{"+num2str(nxi)+num2str(nyi)+"})")
                ylim('auto')
                if arg.ZeroCentered
                    ylim(max(abs(ylim)).*[-1 1])
                end

                grid off
                grid minor

                % imaginary part
                % ax((2*(nxi-1)+1)*ny+nyi)=subplot(2*nx,ny,(2*(nxi-1)+1)*ny+nyi);

                nexttile(T,sub2ind([2*ny, nx], nyi*2, nxi))
                if arg.hold
                    hold on
                end
                plot(t,imag(squeeze(Mt(nxi,nyi,:))),arg.linetype)
                ylabel("Im(a_{"+num2str(nxi)+num2str(nyi)+"})")
                ylim('auto')
                if arg.ZeroCentered
                    ylim(max(abs(ylim)).*[-1 1])
                end

                if nxi == nx
                    xlabel(xlabelStr)
                end

            end
            grid off
            grid minor


        end
    end
    if isempty(arg.title)
        % sgtitle('M(t), vue explosée de la matrice')
    else
        sgtitle(arg.title)
    end


    ax = findall(gcf,'Type','axes');

    % Link the axes if the user specified it
    % Possible values for arg.linkaxes are : 'x', 'y', 'z', 'xy', 'yx', 'yz', 'zy', 'xyz', 'yxz', 'xzy', 'zxy', 'zyx', 'yzx'
    uuu_y={'y','xy','yx','yz','zy','xyz','yxz','xzy','zxy','zyx','yzx'};
    uuu_x={'x','xy','yx','xz','zx','xyz','yxz','xzy','zxy','zyx','yzx'};
    uuu_z={'z','zy','yz','xz','zx','xyz','yxz','xzy','zxy','zyx','yzx'};

    % Set the same y limits for all the subplots if the user specified it or if the user wants to link the y axes
    if arg.GlobalYLim || any(strcmp(arg.linkaxes,uuu_y))
        uu=max(abs( cell2mat(ylim(ax))),[],'all');
        set(ax,'ylim',uu*[-1,1])
    end

    % Set the same x limits for all the subplots if the user specified it or if the user wants to link the x axes
    if arg.plot3D
        if arg.GlobalYLim || any(strcmp(arg.linkaxes,uuu_x))
            uu=max(abs( cell2mat(xlim(ax))),[],'all');
            set(ax,'xlim',uu*[-1,1])
            Link = linkprop(ax(:),{'CameraUpVector', 'CameraPosition', ...
                'CameraTarget', 'ZLim'});
            setappdata(gcf, 'StoreTheLink', Link);
        end
    end

    linkaxes(ax,arg.linkaxes);
else

    % If the plot is not exploded, plot the matrix coefficients in a single figure
    if xor(arg.DispImag,arg.DispReal)
        nexttile(T,1)
        if arg.hold
            hold on
        end
        if arg.DispImag
            plot(t,imag(reshape(Mt,[],numel(t))),arg.linetype)
            if isempty(arg.title)
                title('M(t), imaginary part')
            else
                title(arg.title)
            end
        else
            plot(t,real(reshape(Mt,[],numel(t))),arg.linetype)
            if isempty(arg.title)
                title('M(t), real part')
            else
                title(arg.title)
            end
        end

        ylim('auto')
        if arg.ZeroCentered
            %arg.zeroCentered checks if the plot should be centered on 0
            ylim(max(abs(ylim)).*[-1 1])
        end

    else
        % If both parts are displayed, plot the real and imaginary parts of the matrix coefficients superposed in two subplots
        nexttile(T,1)
        if arg.hold
            % arg.hold checks if the plot should be superposed on the current figure
            hold on
        end
        plot(t,real(reshape(Mt,[],numel(t))),arg.linetype)

        ylim('auto')
        if arg.ZeroCentered
            ylim(max(abs(ylim)).*[-1 1])
        end
        title('M(t), real part')

        nexttile(T,2)
        if arg.hold
            % arg.hold checks if the plot should be superposed on the current figure
            hold on
        end
        plot(t,imag(reshape(Mt,[],numel(t))),arg.linetype)

        ylim('auto')
        if arg.ZeroCentered
            % make y limits symetric
            ylim(max(abs(ylim)).*[-1 1])
        end
        title('M(t), imag part')
    end
end
end