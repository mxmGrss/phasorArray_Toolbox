function r=evalTime(o1,T,t)
    % EVALTIME Evaluate a `PhasorArray` in the time domain for a given period `T`.
    %
    %   EVALTIME(o1, T, t) computes the time-domain representation of a `PhasorArray`,
    %   assuming it is `T`-periodic. The function returns `A(t)`, evaluated at specified
    %   time instants `t`, without generating any plots.
    %
    %   Inputs:
    %     o1  - (PhasorArray) The `PhasorArray` object to be evaluated.
    %     T   - (double, optional) The period of the `PhasorArray`.
    %              - Default: `2*pi`.
    %     t   - (vector or scalar, optional) Time instants for evaluation.
    %              - If `t = []`: Uses a default grid `0:dt:T-dt`, with `dt = T/(20*h)`, where `h` is the highest harmonic.
    %              - If `t = [tmin tmax]`: Uses `t = tmin:dt:tmax` with automatic step size.
    %              - If `t` is a vector: Evaluates `A(t)` at the specified values.
    %              - If `t` is a scalar: Evaluates `A(t)` at a single instant.
    %
    %   Outputs:
    %     r - (m × n × length(t) array) Evaluated time-domain representation of `o1`.
    %
    %   Behavior:
    %     - Uses `plot(o1, T, t, ...)` internally but **without plotting** (`'plot', false`).
    %     - Automatically generates a time grid if `t = []`.
    %     - Designed for numerical evaluation without visualization.
    %
    %   Example Usage:
    %     % Evaluate A(t) over one period
    %     A = PhasorArray(rand(3,3,11));
    %     r = evalTime(A, 2*pi, []);
    %
    %     % Evaluate A(t) on a custom time range
    %     t = linspace(0, 2*pi, 100);
    %     r = evalTime(A, 2*pi, t);
    %
    %   See also: plot, PhasorArray2time.
    arguments
        o1
        T=2*pi
        t=[]
    end
    argo=struct;
    argo.plot=false;
    argo.explosed=true;
    argo.hold=false;
    argo.DispImag=false;
    argo.DispReal=true;
    argo.ZeroCentered=false;
    argo.title=[];
    C=namedargs2cell(argo);
    r=plot(o1,T,t,C{:});
end
