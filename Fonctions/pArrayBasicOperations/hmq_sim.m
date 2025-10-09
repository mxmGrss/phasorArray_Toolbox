function [x,t,dx] = hmq_sim(Aph,tfinal,x0,T,Uph,varg)
    % HMQ_SIM Simulate the response of a time-periodic system.
    %
    %   HMQ_SIM(Aph, tfinal, x0, T, Uph, varg) simulates the system:
    %       dx/dt = A(t)x + U(t)
    %   where `A(t)` and `U(t)` are `T`-periodic matrices, represented as phasor arrays.
    %
    %   The function provides multiple **ODE solvers**, handles **real-valued or complex-valued matrices**, 
    %   and allows **automatic step size selection** based on harmonics resolution.
    %
    %   Syntax:
    %     [x, t] = HMQ_SIM(Aph, tfinal, x0, T)
    %     [x, t, dx] = HMQ_SIM(Aph, tfinal, x0, T, Uph, Name, Value)
    %
    %   Inputs:
    %     Aph    - (PhasorArray) The periodic system matrix `A(t)`, stored as a **3D phasor array**.
    %     tfinal - (scalar or vector, optional) Final simulation time.
    %                - Default: `10*T`
    %                - If scalar: Simulates from `t=0` to `t=tfinal`.
    %                - If `[tmin tmax]`: Simulates from `tmin` to `tmax`.
    %                - If vector: Uses provided time grid.
    %     x0     - (vector, optional) Initial condition `x(0)`. 
    %                - Default: `ones(size(Aph,1),1)`.
    %     T      - (double, optional) The period of `A(t)`. Default: `1`.
    %     Uph    - (PhasorArray or matrix, optional) The time-varying input matrix `U(t)`.
    %                - Default: `[]` (zero input).
    %
    %   Name-Value Pair Arguments:
    %     'odeOpts'           - (struct) Options for the ODE solver (e.g., tolerances). Default: `[]`.
    %     'plot'              - (logical) Plot the state trajectory. Default: `true`.
    %     'solver'            - (char) ODE solver method. Default: `'adaptative'`. Options:
    %                              - `'adaptative'`, `'forward-euler'`, `'RK4'`, `'heuns'`, `'midpoint'`, 
    %                                `'leapfrog'`, `'backward-euler'`, `'trapezoidal'`, `'adams-bashforth'`.
    %     'FSprecpow'         - (integer) Power of 2 for frequency sampling. Default: `8`.
    %     'checkReal'         - (logical) Check if `A(t)` and `U(t)` are real-valued. Default: `false`.
    %     'isRealValued'      - (logical) Force real-valued computation. Default: `auto-detected`.
    %     'providedPhasorForm' - (char) Phasor representation (`"exp"` or `"SinCos"`). Default: `"exp"`.
    %
    %   Outputs:
    %     x  - (matrix) State trajectory of the system (`size(x,1) = size(Aph,1)`).
    %     t  - (vector) Time instants at which `x(t)` is evaluated.
    %     dx - (matrix, optional) Derivative of the state trajectory (only if `nargout > 2`).
    %
    %   Behavior:
    %     - If `isRealValued` is **not provided**, it is automatically detected from `Aph` and `Uph`.
    %     - **Default solver:** `ode15s` (adaptive).
    %     - If `tfinal = 0`, the function **automatically sets** `tfinal = 10*T`.
    %     - If `Uph` is empty, the system is simulated as **homogeneous** (`dx/dt = A(t)x`).
    %     - The highest resolved harmonic `h` determines the **integration step size**.
    %
    %   Example Usage:
    %     % Simulate system response over one period
    %     A = PhasorArray(rand(3,3,11));
    %     [x, t] = hmq_sim(A, 2*pi, ones(3,1), 2*pi);
    %
    %     % Simulate with a time-varying input U(t)
    %     U = PhasorArray(rand(3,1,11));
    %     [x, t] = hmq_sim(A, 2*pi, ones(3,1), 2*pi, U);
    %
    %     % Use a fixed-step RK4 method
    %     [x, t] = hmq_sim(A, 2*pi, ones(3,1), 2*pi, [], 'solver', 'RK4');
    %
    %   See also: PhasorArray2time, evalTime, ode15s.
arguments
    Aph
    tfinal=0
    x0=ones(size(Aph,1),1)
    T=1
    Uph=[]
    varg.odeOpts = {}
    varg.plot logical = true
    varg.solver {mustBeMember(varg.solver,{'adaptative','forward-euler','RK4', 'heuns', 'midpoint', 'leapfrog', 'backward-euler', 'trapezoidal', 'adams-bashforth'})} ='adaptative'
    varg.FSprecpow = 8
    varg.checkReal  logical = false
    varg.isRealValued logical = []
    varg.providedPhasorForm {mustBeMember(varg.providedPhasorForm,["exp","SinCos"])} = "exp"
end

if isempty(varg.isRealValued)
    try
        AA = PhasorArray(Aph);
        UU = PhasorArray(Uph);
        if isreal(AA) && isreal(UU)
            varg.isRealValued=true;
            warning('HMQ_SIM:RealValuedMatrix',"In call to hmq_sim, Matrices appears to be real-valued, to enforce complex-valued computation, set 'isRealValued' to false")
        else
            varg.isRealValued=false;
        end
    catch
        varg.isRealValued=false;
    end
end



if tfinal==0
    tfinal=10*T;
end

% ensure that the input matrices are valid and convert them to numerical arrays if they are symbolic or sdp variables
[Aph, Uph] = validateAphUph(Aph, Uph);

if  varg.isRealValued
    switch varg.providedPhasorForm
        case "exp"
            % Compute the cosine-sine decomposition of the matrices Aph and Uph and use it to define the function to be integrated by the ODE solver
            % it will result in a faster computation of the time-varying matrix-vector product
            h = (size(Aph,3)-1)/2;
            Aphr = real(Aph + flip(Aph,3))/2 + 1i * imag(Aph - flip(Aph,3))/2;
            Acs=real(cat(3,1i*(flip(Aphr(:,:,h+2:end),3)-Aphr(:,:,1:h)),Aphr(:,:,h+1),(Aphr(:,:,h+2:end)+flip(Aphr(:,:,1:h),3))));

            h = (size(Uph,3)-1)/2;
            Uphr = real(Uph + flip(Uph,3))/2 + 1i * imag(Uph - flip(Uph,3))/2;
            Ucs=real(cat(3,1i*(flip(Uphr(:,:,h+2:end),3)-Uphr(:,:,1:h)),Uphr(:,:,h+1),(Uphr(:,:,h+2:end)+flip(Uphr(:,:,1:h),3))));
        otherwise
            Acs=Aph;
            Ucs=Uph;
    end
    % Define the function to be integrated by the ODE solver and the functions to evaluate the time-varying matrices Aph and Uph at a given angle or time
    fsim = @(t,y) f_cs(t,y,Acs,T,Ucs);
    ATt =@(T,t) evalTimeSinCos(Acs,2*pi/T*t);
    UTt =@(T,t) evalTimeSinCos(Ucs,2*pi/T*t);
else
    % Define the function to be integrated by the ODE solver and the functions to evaluate the time-varying matrices Aph and Uph at a given angle or time
    ATt =@(T,t) evalTimeCmplxTt(Aph,T,t);
    UTt =@(T,t) evalTimeCmplxTt(Uph,T,t);
    fsim = @(t,y) f_ph2(t,y,Aph,T,Uph );
end

%
nh=(size(Aph,3)-1)/2;

% Define the time vector for the simulation
% define time step for simulation, higher than the highest frequency of the system and 2^FSprecpow
% selection of power of 2 for the frequency sampling
p = nextpow2(max(nh*8,2^varg.FSprecpow));
dt_sim=T/2^p;
if isscalar(tfinal)
    t=0:dt_sim:tfinal;
elseif numel(tfinal)==2
    t=tfinal(1):dt_sim:tfinal(2);
else
    t= tfinal;
end

switch varg.solver
    case 'adaptative'
        odeOpts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',dt_sim,'Stats','on','Jacobian',@(t,y) ATt(T,t));
        if nargout>2
            rDotDot = [];
            odeOpts=odeset(odeOpts,OutputFcn=@outputFcn);
        end
        if ~isempty(varg.odeOpts)
            if isa(varg.odeOpts,"cell")
                varg.odeOpts=odeset(varg.odeOpts{:});
            end
            odeOpts=odeset(odeOpts,varg.odeOpts);
        end

        %future version of matlab will have a better way to handle this problem with the new ode class
        %F = handleOdeProblem();

        [t,x] = ode15s(fsim,t,double(x0),odeOpts);

        x=x.';
    case 'forward-euler'
        x=zeros(numel(x0),numel(t));
        x(:,1)=x0;
        for ti =2:numel(t)
            dyA = ATt(T,t(ti-1))*x(:,ti-1);
            dyB = UTt(T,t(ti-1));
            x(:,ti)=(dyA+dyB)*dt_sim+x(:,ti-1);
        end

    case 'RK4'
        x=zeros(numel(x0),numel(t));
        x(:,1)=x0;
        for ti =2:numel(t)
            tn = t(ti-1);
            xn = x(:,ti-1);

            k1 = fsim(tn,xn);
            k2 = fsim(tn+dt_sim/2,xn+dt_sim/2*k1);
            k3 = fsim(tn+dt_sim/2,xn+dt_sim/2*k2);
            k4 = fsim(tn+dt_sim,xn+dt_sim/2*k3);

            x(:,ti)=xn+dt_sim/6*(k1+2*k2+2*k3+k4);
        end

    case 'heuns'
        x=zeros(numel(x0),numel(t));
        x(:,1)=x0;
        for ti =2:numel(t)
            tn = t(ti-1);
            xn = x(:,ti-1);

            k1 = ATt(T,tn)*xn+UTt(T,tn);
            k2 = ATt(T,tn+dt_sim)*(xn+dt_sim*k1)+UTt(T,tn+dt_sim);

            x(:,ti)=xn+dt_sim/2*(k1+k2);
        end

    case 'midpoint'
        x=zeros(numel(x0),numel(t));
        x(:,1)=x0;
        for ti =2:numel(t)
            tn = t(ti-1);
            xn = x(:,ti-1);

            k1 = ATt(T,tn)*xn+UTt(T,tn);
            k2 = ATt(T,tn+dt_sim/2)*(xn+dt_sim/2*k1)+UTt(T,tn+dt_sim/2);

            x(:,ti)=xn+dt_sim*k2;
        end

    case 'leapfrog'
        x=zeros(numel(x0),numel(t));
        x(:,1)=x0;
        x_half = x(:,1) + 0.5*dt_sim*(ATt(T,0)*x(:,1)+UTt(T,0));
        for ti =2:numel(t)
            tn = t(ti-1);
            xn = x_half;

            k = ATt(T,tn)*xn+UTt(T,tn);
            x(:,ti)=x(:,ti-1)+dt_sim*k;
            x_half = x(:,ti) + 0.5*dt_sim*(ATt(T,tn+dt_sim)*x(:,ti)+UTt(T,tn+dt_sim));
        end
    case 'backward-euler'
        x=zeros(numel(x0),numel(t));
        x(:,1)=x0;
        for ti =2:numel(t)
            dyA = ATt(T,t(ti))*x(:,ti-1);
            dyB = UTt(T,t(ti));
            x(:,ti)=(dyA+dyB)*dt_sim+x(:,ti-1);
        end

    case 'trapezoidal'
        x=zeros(numel(x0),numel(t));
        x(:,1)=x0;
        for ti =2:numel(t)
            tn = t(ti-1);
            xn = x(:,ti-1);

            k1 = ATt(T,tn)*xn+UTt(T,tn);
            k2 = ATt(T,tn+dt_sim)*(xn+dt_sim*k1)+UTt(T,tn+dt_sim);

            x(:,ti)=xn+dt_sim/2*(k1+k2);
        end

    case 'adams-bashforth'
        x=zeros(numel(x0),numel(t));
        x(:,1)=x0;
        for ti =3:numel(t)
            tn = t(ti-1);
            xn = x(:,ti-1);

            k1 = ATt(T,tn-dt_sim)*x(:,ti-2)+UTt(T,tn-dt_sim);
            k2 = ATt(T,tn)*xn+UTt(T,tn);

            x(:,ti)=xn+dt_sim/2*(3*k2-k1);
        end
end



if varg.plot
    if ishold
        holdvar='on';
    else
        holdvar='off';
    end
    for iii = 1:(1):size(x,1)
        subplot(size(x,1),1+(nargout>2),(iii-1)*(1+(nargout>2))+1)
        hold off
        hold(holdvar);
        plot(t,x(iii,:))
        hold off
        if nargout>2
            subplot(size(x,1),1+(nargout>2),(iii-1)*(1+(nargout>2))+2)
            hold(holdvar);
            plot(t,rDotDot(iii,:))
            hold off
        end
    end

end

if nargout>2
    dx=rDotDot;
end

    function status = outputFcn(t,q,flag)
        try
            status = 0;
            % at initialization, and after each succesful step
            if  strcmp(flag, 'init')
                for ii = 1
                    deriv = fsim(t(ii),q(:,ii));
                    rDotDot(:,end+1) = deriv;
                end
            end
            if isempty(flag)
                for ii = 1:numel(t)
                    deriv = fsim(t(ii),q(:,ii));
                    rDotDot(:,end+1) = deriv;
                end
            end
        catch e
            t
            q
            flag
            error(e)
        end
    end % output function


    function F = handleOdeProblem()
        %for future version of matlab with class ODE, will allow more flexibility in the solver choice and the output
        F = ode ;
        F.ODEFcn = fsim;
        F.InitialValue = x0 ;
        F.Jacobian = @(t,y) ATt(T,t);
    end

end


function dydt = f_ph(t,y,Aph,T,Uph,checkReal )
y=y(:);
dydt = PhasorArray2time(Aph,T,t,checkReal =checkReal )*y+PhasorArray2time(Uph,T,t,checkReal =checkReal );
dydt=dydt(:);
end
function dydt = f_ph2(t,y,Aph,T,Uph )
y=y(:);
dydt = evalTimeCmplxTt(Aph,T,t)*y+evalTimeCmplxTt(Uph,T,t);
dydt=dydt(:);
end


function dydt = f_cs(t,y,Aph,T,Uph)
y=y(:);
angle = 2*pi/T*t;
dydt = evalTimeSinCos(Aph,angle)*y+evalTimeSinCos(Uph,angle);
dydt=dydt(:);
end

function Mt = evalTimeCmplxTt(Phas,T,t)
angle = 2*pi/T*t;
h = (size(Phas,3)-1)/2;
eit=exp(1i*(-h:h)'*angle);
Mt=tensorprod(Phas,double(eit),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
end

function Mt = evalTimeCmplxA(Phas,angle)
h = (size(Phas,3)-1)/2;
eit=exp(1i*(-h:h)'*angle);
Mt=tensorprod(Phas,double(eit),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
end

function Mt = evalTimeSinCos(PhasSC,angle)
h = (size(PhasSC,3)-1)/2;
eit=[sin((h:-1:1)'*angle); cos((0:h)'*angle) ];
Mt=tensorprod(PhasSC,double(eit),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
end

function [Aph, Uph] = validateAphUph(Aph, Uph)
% Validate Aph
if isa(Aph, 'PhasorArray')
    Aph = Aph.Value;
elseif isa(Aph, "ndsdpvar") || isa(Aph, 'sdpvar')
    Aph = value(Aph);
end

% Validate Uph
if isempty(Uph)
    Uph = zeros(size(Aph, 1), 1);
else
    if isa(Uph, 'PhasorArray')
        Uph = Uph.Value;
    elseif isa(Uph, "ndsdpvar") || isa(Uph, 'sdpvar')
        Uph = value(Uph);
    end
end
end

