function r=plot3D(o1,T,t,arg)
    %PLOT3D Produce a 3D plot where x-axis is the real part, y-axis is the imaginary part, and z-axis is time
    %
    %   r = PLOT3D(o1, T, t, arg) produces a 3D plot where the x-axis is the real part of A(t),
    %   the y-axis is the imaginary part of A(t), and the z-axis is time.
    %   For each element of A, a subplot is created.
    %
    %   Input Arguments:
    %   o1 - The PhasorArray object to be plotted.
    %   T  - The period of the PhasorArray. Default is 1.
    %   t  - Time instants at which to evaluate the PhasorArray. Can be:
    %       - [] (empty), then t takes the value [0 T].
    %       - [tmin tmax], then t=tmin:dt:tmax with automatic discretization.
    %       - A vector on which A(t) is evaluated.
    %       - A single scalar, then t takes the value [0 t].
    %   arg - (Optional) Name-value pair arguments:
    %       'ZeroCentered' - Logical flag to normalize x and y axes around zero. Default is false.
    %       'title' - String to display a custom title to the figure. Default is [].
    %       'GlobalYLim' - Logical flag to enforce same Y and X limits on the axes. Default is false.
    %       'linkaxes' - String to link zoom on the x, y, z axes. Default is 'x'.
    %
    %   Output Arguments:
    %   r - The evaluated PhasorArray at the specified time instants.
    %
    %   Example:
    %   r = plot3D(o1, 2*pi, []);
    %   r = plot3D(o1, 2*pi, [0 2*pi], 'ZeroCentered', true, 'title', '3D Plot');
    %
    %   See also: PhasorArray2time
    arguments
        o1
        T=1
        t=[0 T]
        arg.ZeroCentered=false
        arg.title=[]
        arg.GlobalYLim=false
        arg.linkaxes='x'
    end
    if numel(t)==1
        t=sort([0 t]);
    end
    rr=PhasorArray2time(o1,T,t,plot=true,plot3D=true, DispImag=false, DispReal=true,explosed=true,ZeroCentered=arg.ZeroCentered,title=arg.title,GlobalYLim=arg.GlobalYLim,linkaxes=arg.linkaxes);

    if nargout
        r=rr;
    end

end

