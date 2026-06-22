function [r,t] = plot(o1,T,t,arg)
    % PLOT Evaluate and plot a T-periodic `PhasorArray`.
    %
    %   PLOT(o1, T, t, arg) computes and plots the time-domain representation of a
    %   `PhasorArray`, assuming it is `T`-periodic.
    %
    %   This method evaluates the time-domain representation of a `PhasorArray` and
    %   plots it over a specified period `T`. It provides flexible options for
    %   visualization, including real/imaginary decomposition, subplot arrangements,
    %   and axis linking.
    %
    %   Inputs:
    %     o1  - (PhasorArray) The `PhasorArray` object to be evaluated and plotted.
    %     T   - (double, optional) Period of the `PhasorArray`.
    %              - Default: `1`.
    %     t   - (vector or scalar, optional) Time instants for evaluation.
    %              - If `t = []`: A default time grid is generated as `0:dt:T-dt`, where `dt = T/(20*h)`, with `h` the highest harmonic.
    %              - If `t = [tmin tmax]`: Uses `t = tmin:dt:tmax` with `dt` computed as above.
    %              - If `t` is a vector: Evaluates `A(t)` at the specified values.
    %              - If `t` is a scalar: Evaluates `A(t)` at that single instant.
    %     Name-Value Pair Arguments:
    %              - 'plot' (logical): Display the plot (default: `true`).
    %              - 'explosed' (logical): Plot each matrix component in a separate subplot (default: `true`).
    %              - 'hold' (logical): Hold the current plot (default: `false`).
    %              - 'DispImag' (logical): Display the imaginary part of the matrix (default: `false`).
    %              - 'DispReal' (logical): Display the real part of the matrix (default: `true`).
    %              - 'ZeroCentered' (logical): Center the Y-axis around zero (default: `false`).
    %              - 'title' (string): Custom title for the figure (default: `[]`).
    %              - 'linetype' (string): Line style for the plot (default: `'-'`).
    %              - 'GlobalYLim' (logical): Apply the same Y-limits across subplots (default: `false`).
    %              - 'linkaxes' (string): Link axes of subplots (`'x'`, `'y'`, `'xy'`, etc.) (default: `'x'`).
    %              - 'forceReal' (logical): Assume `o1` is real-valued and simplify computation (default: `false`).
    %
    %   Outputs:
    %     r - (m × n × length(t) array) Evaluated time-domain representation of `o1`.
    %     t - (vector) Time instants at which `o1` is evaluated.
    %
    %   Behavior:
    %     - **Evaluation**: Computes `A(t)` using `PhasorArray2time`, which performs an inverse Fourier summation.
    %     - **Automatic Time Grid**: If `t = []`, it generates a default time grid based on `T` and harmonics.
    %     - **Real-Valued Constraints**: If `o1` is real-valued, `forceReal` is automatically set to `true`.
    %     - **Plot Customization**:
    %         - Supports separate subplots for each matrix element (`explosed = true`).
    %         - Can plot only the real or imaginary part (`DispReal`, `DispImag`).
    %         - Allows axis linking and uniform Y-axis scaling across subplots.
    %
    %   Example Usage:
    %     % Evaluate and plot a PhasorArray over one period
    %     A = PhasorArray(rand(3,3,11));
    %     plot(A, 2*pi, []);
    %
    %     % Evaluate A(t) on a custom time range
    %     t = linspace(0, 2*pi, 100);
    %     plot(A, 2*pi, t, 'plot', true, 'explosed', false);
    %
    %     % Plot only the imaginary part with a different linestyle
    %     plot(A, 2*pi, [], 'DispImag', true, 'DispReal', false, 'linetype', '--');
    %
    %   See also: PhasorArray2time.
    arguments
        o1
        T=1
        t=[]
        arg.plot logical =true
        arg.explosed logical =true
        arg.hold logical =false
        arg.DispImag logical = []
        arg.DispReal logical = []
        arg.ZeroCentered logical =false
        arg.title = [] %options : "none" or string or []
        arg.linetype='-'
        arg.GlobalYLim logical =false
        arg.linkaxes='x'
        arg.forceReal = false
        arg.grid = 'on'
    end

    if isscalar(o1)
        arg.explosed = false;
    end

    if isspecial(o1)
        o1 = PhasorArray(value(o1.value));
    end
    if ishold
        arg.hold = true;
    end
    if isreal(o1)
        arg.forceReal = true;
        if isempty(arg.DispImag)
            arg.DispImag = false;
        end
        if isempty(arg.DispReal)
            arg.DispReal = true;
        end
    else
        if isempty(arg.DispImag)
            arg.DispImag = true;
        end
        if isempty(arg.DispReal)
            arg.DispReal = true;
        end
    end
    [rr,tt]=PhasorArray2time(o1,T,t,plot=arg.plot, DispImag=arg.DispImag, ...
        DispReal=arg.DispReal,explosed=arg.explosed,hold=arg.hold,ZeroCentered=arg.ZeroCentered, ...
        title=arg.title,linetype=arg.linetype,GlobalYLim=arg.GlobalYLim,linkaxes=arg.linkaxes,forceReal=arg.forceReal,grid = arg.grid);

    if nargout>0
        r=rr;
    end
    if nargout>1
        t=tt;
    end
end

