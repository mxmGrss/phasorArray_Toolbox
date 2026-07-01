function [Phi0T,E,V,D,W] = TransitionMatrixOverT(Aph,T,varg)
%TransitionMatrixOverT Compute the state transition matrix over one period.
%   This function simulates the linear periodic system \dot{x} = A(t)x 
%   from t=0 to t=T for each fundamental initial condition (x0 = e_i)
%   and concatenates the results to form the transition matrix Phi(T,0).
%
%   Syntax:
%       Phi0T = TransitionMatrixOverT(Aph, T)
%       [Phi0T, E, V, D, W] = TransitionMatrixOverT(Aph, T, Name, Value)
%
%   Input Arguments:
%       Aph : PhasorArray or 3D double - The periodic state matrix A(t).
%       T   : double - The period of the system (default: 1).
%
%   Name-Value Arguments:
%       simutype     : char - Time integration solver {'adaptative', 'forward-euler', 'RK4'} (default: 'adaptative').
%       FSprecpow    : double - Power defining the fixed step size dt = T/2^p for fixed solvers (default: 0).
%       enforce_real : logical - Force complex-valued ODE integration even if A(t) appears real.
%
%   Output Arguments:
%       Phi0T : matrix - The state transition matrix evaluated at t=T (monodromy matrix).
%       E     : vector - Eigenvalues of Phi0T (Floquet multipliers).
%       V     : matrix - Right eigenvectors of Phi0T.
%       D     : matrix - Diagonal matrix of eigenvalues of Phi0T.
%       W     : matrix - Left eigenvectors of Phi0T.
%
%   See also: FloquetDec, freqVarying_hmqSim

arguments
    Aph
    T (1,1) double = 1
    varg.simutype (1,:) char {mustBeMember(varg.simutype,{'adaptative','forward-euler','RK4'})} = 'adaptative'
    varg.FSprecpow (1,1) double = 0
    varg.enforce_real logical = []
end

if isempty(varg.enforce_real)
    try
        AA = PhasorArray(Aph);
        if isreal(AA)
            varg.enforce_real=true;
            warning('TransitionMatrixOverT:appearsRealValued', 'Matrix A appears to be real-valued. To enforce complex-valued computation, set ''isRealValued'' to false.')
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