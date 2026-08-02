function [y,t,dy]=initial(o1,x0,T,tfinal,nvp)
    % INITIAL Compute the system response to an initial state in a periodic state-space model.
    %
    %   INITIAL(o1, x0, T, tfinal, nvp) simulates the system response for:
    %       dx/dt = A(t)x
    %   where `A(t)` is a `T`-periodic state-space matrix. This function calls `lsim` to compute the response.
    %
    %   Syntax:
    %     [y, t] = INITIAL(A, x0, T, tfinal)
    %     [y, t, dy] = INITIAL(A, x0, T, tfinal, Name, Value)
    %
    %   Inputs:
    %     o1     - (PhasorArray) The time-varying system matrix `A(t)`, stored as a **3D phasor array**.
    %     x0     - (vector, optional) Initial state `x(0)`.
    %                - Default: `ones(size(o1,1),1)`.
    %     T      - (double, optional) Period of the system. Default: `2*pi`.
    %     tfinal - (scalar or vector, optional) Final simulation time.
    %                - Default: `10*T`.
    %                - If scalar: Simulates from `t=0` to `t=tfinal`.
    %                - If `[tmin tmax]`: Simulates from `tmin` to `tmax`.
    %                - If vector: Uses provided time grid.
    %
    %   Name-Value Pair Arguments:
    %     'opts'             - (struct) Options for the ODE solver (`odeset`). Default: `[]`.
    %     'plot'             - (logical) Plot the state trajectory. Default: `true`.
    %     'solver'           - (char) ODE solver method. Default: `'adaptative'`. Options:
    %                              - `'adaptative'`, `'forward-euler'`, `'RK4'`
    %     'FSprecpow'        - (integer) Power of 2 for frequency sampling. Default: `8`.
    %     'checkReal'        - (logical) Force real-valued output (only if system is **guaranteed** real). Default: `false`.
    %     'isRealValued'     - (logical) Force real-valued computation. Default: `false`.
    %
    %   Outputs:
    %     y  - (matrix) State trajectory of the system (`size(y,1) = size(o1,1)`).
    %     t  - (vector) Time instants at which `y(t)` is evaluated.
    %     dy - (matrix, optional) Derivative of the state trajectory (only if `nargout > 2`).
    %
    %   Behavior:
    %     - This function is a **wrapper for `lsim`**, automatically setting `U(t) = 0`.
    %     - If `isRealValued` is **not provided**, it is automatically detected from `o1`.
    %     - If `tfinal = 0`, the function **automatically sets** `tfinal = 10*T`.
    %
    %   Example Usage:
    %     % Simulate free response over one period
    %     A = PhasorArray(rand(3,3,11));
    %     [y, t] = initial(A, ones(3,1), 2*pi, 2*pi);
    %
    %     % Simulate with a fixed-step RK4 method
    %     [y, t] = initial(A, ones(3,1), 2*pi, 2*pi, 'solver', 'RK4');
    %
    %   See also: LSIM, HMQ_SIM.
    arguments
        o1
        x0=ones(size(o1,1),1)
        T=2*pi
        tfinal=0
        nvp.opts=[] %odeset option
        nvp.plot=true
        nvp.solver='adaptative'
        nvp.FSprecpow=8
        nvp.checkReal=0
        nvp.isRealValued logical = false
    end
    if isreal(o1)
        nvp.isRealValued = true;
    end
    vvarg = namedargs2cell(nvp);
    if nargout == 3
        [y,t,dy]=lsim(o1,tfinal,x0,T,vvarg{:});
    else
        [y,t]=lsim(o1,tfinal,x0,T,vvarg{:});
    end

end
