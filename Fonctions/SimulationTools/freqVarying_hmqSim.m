function [x,t,dx] = freqVarying_hmqSim(Aph,tfinal,x0,omega,Bph,u,varg)
%freqVarying_hmqSim simulate free-response of system \dot x = A(t) x where A(theta) is
%a 2pi periodique matrices, with dot theta = omega(t) ans arbitrary function
%, from initial condition x0.
%   freqVarying_hmqSim(Aph, tmax,x0,omega)

arguments
    Aph
    tfinal=0
    x0=ones(size(Aph,1),1)
    omega = @(t) 2 * pi
    Bph=[]
    u = @(y,t) 1;
    varg.opts = {}
    varg.plot = true
    varg.solver {mustBeMember(varg.solver,{'adaptative','forward-euler','RK4'})} ='adaptative'
    varg.FSprecpow = 8
    varg.checkReal  logical = false
    varg.isRealValued logical = false
end

opts=varg.opts;

if isscalar(tfinal)&& tfinal==0
    tfinal=10;
end
if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end

if isa(Aph,"ndsdpvar") || isa(Aph,'sdpvar')
    Aph=value(Aph);
end

if isempty(Bph)
    Bph=zeros(size(Aph,1),1);
else
    if isa(Bph,'PhasorArray')
        Bph=Bph.Value;
    end
    if isa(Bph,"ndsdpvar") || isa(Bph,'sdpvar')
        Bph=value(Bph);
    end
end
switch varg.solver
    case 'adaptative'
        nh=(size(Aph,3)-1)/2;
        dt_sim=1/(2^nextpow2(max(nh*8,2^varg.FSprecpow)))
        opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',dt_sim,'Stats','on');
        if nargout>2
            rDotDot = [];
            opts=odeset(opts,OutputFcn=@outputFcn);
        end
        if ~isempty(varg.opts)
            if isa(varg.opts,"cell")
            varg.opts=odeset(varg.opts{:});
            end
            opts;
            opts=odeset(opts,varg.opts);
        end

        if isscalar(tfinal)
            tfinal=[0:dt_sim:tfinal];
        end
        if  varg.isRealValued

            h = (size(Aph,3)-1)/2;
            Aphr = real(Aph + flip(Aph,3))/2 + 1i * imag(Aph - flip(Aph,3))/2;
            Acs=real(cat(3,1i*(flip(Aphr(:,:,h+2:end),3)-Aphr(:,:,1:h)),Aphr(:,:,h+1),(Aphr(:,:,h+2:end)+flip(Aphr(:,:,1:h),3))));

            h = (size(Bph,3)-1)/2;
            Uphr = real(Bph + flip(Bph,3))/2 + 1i * imag(Bph - flip(Bph,3))/2;
            Ucs=real(cat(3,1i*(flip(Uphr(:,:,h+2:end),3)-Uphr(:,:,1:h)),Uphr(:,:,h+1),(Uphr(:,:,h+2:end)+flip(Uphr(:,:,1:h),3))));
            
            fsim = @(t,y) f_cs(t,y,Acs,omega,Ucs,u);
        else
            fsim = @(t,y) f_ph(t,y,Aph,omega,Bph,u,varg.checkReal );
        end
        [t,x] = ode15s(fsim,tfinal,double(x0),opts);

        x=x.';
    case 'forward-euler'
        nh=(size(Aph,3)-1)/2;
        dt_sim=(1)/(2^(max(ceil(log((nh+1))/log(2)+1),varg.FSprecpow)));
        Tmax=tfinal;%-dt_sim;
        t=0:dt_sim:Tmax;
        x=zeros(numel(x0),numel(t));
        phas=zeros(1,numel(t));
        x(:,1)=x0(1:end-1);
        phas(1)=x0(end);
        w(:,1)=u(phas(1),t(1));
        for ti =2:numel(t)
            dyA=evalTimeCmplx(Aph,phas(ti-1))*x(:,ti-1);
            dyB=evalTimeCmplx(Bph,phas(ti-1))*w(:,ti-1);
            x(:,ti)=(dyA+dyB)*dt_sim+x(:,ti-1);
            w(:,ti)=u([x(:,ti);phas(ti)],t(ti));
            phas(ti)=phas(ti-1)+dt_sim*omega(t(ti-1));
        end

    case 'RK4'
        nh=(size(Aph,3)-1)/2;
        %     max(ceil(log((nh+1))/log(2)+1),8)
        dt_sim= ((1)/((2^(max(ceil(log((nh+1))/log(2)+1),varg.FSprecpow)))+0));
        Tmax=tfinal;%-dt_sim;
        t=0:dt_sim:Tmax;
        x=zeros(numel(x0)-1,numel(t));
        y=zeros(numel(x0),numel(t));
        phas=zeros(1,numel(t));
        y(:,1)=x0;
        x(:,1)=x0(1:end-1);
        phas(1)=x0(end);
        w(:,1)=u(phas(1),t(1));
        for ti =2:numel(t)
            tn = t(ti-1);
            xn = x(:,ti-1);
            phin = phas(ti-1);
            wn = w(:,ti-1);
            yn = [xn;phin];

            k1_x = evalTimeCmplx(Aph,phin)*xn+evalTimeCmplx(Uph,phin)*u(yn,tn);
            k1_ph = omega(tn);
            k1 = [k1_x;k1_ph];

            tn_k23 = tn + dt_sim/2;
            yn_k2 = yn + dt_sim/2*k1;
            xn_k2 = yn(1:end-1);
            phin_k2 = yn(end);

            k2_x  = evalTimeCmplx(Aph,phin_k2)*(xn_k2) + evalTimeCmplx(Uph,phin_k2)*u(yn_k2,tn_k23);
            k2_ph = omega(tn_k23);
            k2 = [k2_x;k2_ph];

            yn_k3 = yn + dt_sim/2*k2;
            xn_k3 = yn(1:end-1);
            phin_k3 = yn(end);

            k3_x  = evalTimeCmplx(Aph,phin_k3)*(xn_k3) + evalTimeCmplx(Uph,phin_k3)*u(yn_k3,tn_k23);
            k3_ph = omega(tn_k23);
            k3 = [k3_x;k3_ph];

            
            tn_k4 = tn + dt_sim;
            yn_k4 = yn + dt_sim*k3;
            xn_k4 = yn(1:end-1);
            phin_k4 = yn(end);

            k4_x  = evalTimeCmplx(Aph,phin_k4)*(xn_k4) + evalTimeCmplx(Uph,phin_k4)*u(yn_k4,tn_k4);
            k4_ph = omega(tn_k4);
            k4 = [k4_x;k4_ph];


            dy=1/6*(k1+2*k2+2*k3+k4);
            
            y(:,ti)  = dy*dt_sim+yn;
            x(:,ti)  = y(1:end-1,ti);
            phas(ti) = y(end,ti);

        end
end


if varg.plot
    for ii = 1:(1):size(x,1)
        subplot(size(x,1),1+(nargout>2),(ii-1)*(1+(nargout>2))+1)
        plot(t,x(ii,:))
        if nargout>2
            subplot(size(x,1),1+(nargout>2),(ii-1)*(1+(nargout>2))+2)
            plot(t,rDotDot(ii,:))
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

end


function dydt = f_ph(t,y,Aph,omega,Bph,u,checkReal )
y=y(:);
angle=y(end);
dangle = omega(t);
dydt = evalTimeCmplx(Aph,angle)*y(1:end-1)+evalTimeCmplx(Bph,angle)*u(y,t);
dydt(end+1)=dangle;
dydt=dydt(:);
end


function dydt = f_cs(t,y,Aph,omega,Bph,u)
y=y(:);
angle=y(end);
dangle = omega(t);
dydt = evalTimeSinCos(Aph,angle)*y(1:end-1)+evalTimeSinCos(Bph,angle)*u(y,t);
dydt(end+1)=dangle;
dydt=dydt(:);
end

function Mt = evalTimeCmplx(Phas,angle)
h = (size(Phas,3)-1)/2;
eit=exp(1i*(-h:h)'*angle);
Mt=tensorprod(Phas,double(eit),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
end

function Mt = evalTimeSinCos(PhasSC,angle)
h = (size(PhasSC,3)-1)/2;
eit=[sin((h:-1:1)'*angle); cos((0:h)'*angle) ];
Mt=tensorprod(PhasSC,double(eit),3,1); %est un 3D array dont Mt(:,:,k) est M(t(k))
end
