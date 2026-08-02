function [r,t]=sim(o1,T,t,nvp)
    % SIM Evaluate the time-domain representation of a `PhasorArray`.
    %
    %   SIM(o1, T, t, nvp) computes the time-domain representation of a `PhasorArray`,
    %   assuming it is `T`-periodic.
    %
    %   SIM is a **convenience function** for evaluating `A(t)`, leveraging `PhasorArray2time`.
    %   It includes an option to enforce real-valued computation when `A(t)` is known to be real.
    %
    %   Inputs:
    %     o1  - (PhasorArray) The `PhasorArray` object representing a periodic matrix.
    %     T   - (double, optional) The period of `A(t)`.
    %              - Default: `2*pi`.
    %     t   - (vector or scalar, optional) Time instants for evaluation.
    %              - If `t = []`: Uses a default time grid `0:dt:T-dt`, where `dt = T/(20*h)`, and `h` is the highest harmonic.
    %              - If `t = [tmin tmax]`: Uses `t = tmin:dt:tmax` with automatic step size.
    %              - If `t` is a vector: Evaluates `A(t)` at the specified values.
    %              - If `t` is a scalar: Evaluates `A(t)` at a single instant.
    %     nvp - (Optional) Name-value pair arguments:
    %              - 'isRealValued' (logical): Enforce real-valued computation (default: `false`).
    %
    %   Outputs:
    %     r - (m × n × length(t) array) Evaluated time-domain representation of `o1`.
    %     t - (vector) Time instants at which `o1` was evaluated.
    %
    %   Behavior:
    %     - If `o1` is **real-valued**, `isRealValued` is automatically set to `true`.
    %     - Calls `PhasorArray2time` with a structured argument list for evaluation.
    %
    %   Example Usage:
    %     % Evaluate A(t) over one period
    %     A = PhasorArray(rand(3,3,11));
    %     [r, t] = sim(A, 2*pi, []);
    %
    %     % Force real-valued computation
    %     [r, t] = sim(A, 2*pi, [], 'isRealValued', true);
    %
    %   See also: PhasorArray2time, plot.
    arguments
        o1
        T=2*pi
        t=[]
        nvp.isRealValued = false;
    end


    if isreal(o1)
        nvp.isRealValued = true;
    end

    argo=struct;
    argo.plot=false;
    argo.explosed=true;
    argo.hold=false;
    argo.DispImag=false;
    argo.DispReal=true;
    argo.ZeroCentered=false;
    argo.title=[];
    argo.forceReal = nvp.isRealValued;
    C=namedargs2cell(argo);
    [r,t]=PhasorArray2time(o1,T,t,C{:});
end
