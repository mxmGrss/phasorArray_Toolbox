function [Phi0T,E,V,D,W] = TransitionMatrixOverT(Aph,T,varg)
%TransitionMatrixOverT simulate dot x = A(t) x from 0 to T
% with x0 = ei and concatenate the result to forme transition matrix
%the resulting vectors

arguments
    Aph
    T = 1
    varg.simutype {mustBeMember(varg.simutype,{'adaptative','forward-euler','RK4'})} ='adaptative'
    varg.FSprecpow = 0
    varg.enforce_real logical = []
end

if isempty(varg.enforce_real)
    try
        AA = PhasorArray(Aph);
        if isreal(AA)
            varg.enforce_real=true;
            warning('TransitionMatrixOverT:RealValuedMatrix',"In call to TransitionMatrixOverT, Matrix A appears to be real-valued, to enforce complex-valued computation, set 'isRealValued' to false")
        else
            varg.enforce_real=false;
        end
    catch
        varg.enforce_real=false;
    end
end


[nx,ny]=size(Aph,[1 2]);
nh = (size(Aph,3)-1)/2;

if nx~=ny
    error('A(t) must be square matrix')
end

Phi0T = zeros(ny,nx);
p = nextpow2(max(nh*4,2^varg.FSprecpow));
dt_sim=T/2^p;
for nxi = 1:nx
    ei=zeros(nx,1);
    ei(nxi)=1;
    switch varg.simutype
        case 'adaptative'
            opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',dt_sim,'Stats','off');
            [y]=hmq_sim(Aph,0:dt_sim:T,ei,T,odeOpts=opts,plot=false,isRealValued=varg.enforce_real);
            Phi0T(:,nxi)=y(:,end);
        case 'forward-euler'
            [y]=hmq_sim(Aph,0:dt_sim:T,ei,T,plot=false,solver='forward-euler',FSprecpow=varg.FSprecpow,isRealValued=varg.enforce_real);
            Phi0T(:,nxi)=y(:,end);
        case 'RK4'
            [y]=hmq_sim(Aph,0:dt_sim:T,ei,T,plot=false,solver='RK4',FSprecpow=varg.FSprecpow,isRealValued=varg.enforce_real);
            Phi0T(:,nxi)=y(:,end);
    end
end

if nargout>1
    E=eig(Phi0T);
end
if nargout>2
    [V,D,W]=eig(Phi0T);

end