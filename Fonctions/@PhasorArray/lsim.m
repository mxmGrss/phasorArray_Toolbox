function [y,t,dy]=lsim(o1,tfinal,x0,T,Uph,varg)
    % LSIM Simulate the response of a time-periodic linear system.
    %
    %   LSIM(o1, tfinal, x0, T, Uph, varg) simulates the system:
    %       dx/dt = A(t)x + U(t)
    %   where `A(t)` and `U(t)` are `T`-periodic matrices, represented as phasor arrays.
    %
    %   The function supports **adaptive and fixed-step solvers**, allows **real-valued computation**, and
    %   provides **ODE solver customization** via `odeset` options.
    %
    %   Syntax:
    %     [y, t] = LSIM(A, tfinal, x0, T)
    %     [y, t, dy] = LSIM(A, tfinal, x0, T, U, Name, Value)
    %
    %   Inputs:
    %     o1     - (PhasorArray) The periodic system matrix `A(t)`, stored as a **3D phasor array**.
    %     tfinal - (scalar or vector, optional) Final simulation time.
    %                - Default: `10*T`
    %                - If scalar: Simulates from `t=0` to `t=tfinal`.
    %                - If `[tmin tmax]`: Simulates from `tmin` to `tmax`.
    %                - If vector: Uses provided time grid.
    %     x0     - (vector, optional) Initial condition `x(0)`.
    %                - Default: `ones(size(o1,1),1)`.
    %     T      - (double, optional) The period of `A(t)`. Default: `1`.
    %     Uph    - (PhasorArray or matrix, optional) The time-varying input matrix `U(t)`.
    %                - Default: `[]` (zero input).
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
    %     - If `isRealValued` is **not provided**, it is automatically detected from `o1` and `Uph`.
    %     - **Default solver:** `ode15s` (adaptive).
    %     - If `tfinal = 0`, the function **automatically sets** `tfinal = 10*T`.
    %     - If `Uph` is empty, the system is simulated as **homogeneous** (`dx/dt = A(t)x`).
    %     - The highest resolved harmonic `h` determines the **integration step size**.
    %
    %   Example Usage:
    %     % Simulate system response over one period
    %     A = PhasorArray(rand(3,3,11));
    %     [y, t] = lsim(A, 2*pi, ones(3,1), 2*pi);
    %
    %     % Simulate with a time-varying input U(t)
    %     U = PhasorArray(rand(3,1,11));
    %     [y, t] = lsim(A, 2*pi, ones(3,1), 2*pi, U);
    %
    %     % Use a fixed-step RK4 method
    %     [y, t] = lsim(A, 2*pi, ones(3,1), 2*pi, [], 'solver', 'RK4');
    %
    %   See also: HMQ_SIM, ode15s.
    arguments
        o1
        tfinal=0
        x0=ones(size(o1,1),1)
        T=1
        Uph=[]
        varg.opts=[] %odeset option
        varg.plot=true
        varg.solver='adaptative'
        varg.FSprecpow=8
        varg.checkReal=0
        varg.isRealValued logical = false
    end
    if isempty(x0)
        x0=ones(size(o1,1),1);
    end
    if isreal(o1)
        varg.isRealValued = true;
    end

    C=namedargs2cell(varg);
    C{1}="odeOpts";

    %asking derivative trigger more computation from hmq_sim,
    %procede with care
    if nargout>2
        [y,t,dy] = hmq_sim(o1,tfinal,x0,T,Uph,C{:});
    else %sinon
        [y,t]    = hmq_sim(o1,tfinal,x0,T,Uph,C{:});
    end

end

