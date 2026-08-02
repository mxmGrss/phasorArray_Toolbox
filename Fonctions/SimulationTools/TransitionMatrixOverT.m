function [Phi0T,E,V,D,W] = TransitionMatrixOverT(Aph,T,nvp)
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
%       T   : double - The period of the system (default: 2*pi).
%
%   Name-Value Arguments:
%       simutype     : char - Time integration solver {'adaptative', 'forward-euler', 'RK4'} (default: 'adaptative').
%       FSprecpow    : double - Power defining the fixed step size dt = T/2^p for fixed solvers (default: 0).
%       enforce_real : logical - Integrate over the reals. Left empty it is detected from
%                      A: true when A(t) is real, false otherwise. Forwarded to hmq_sim as
%                      its 'isRealValued', which is the same switch under another name.
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
    T (1,1) double = 2*pi
    nvp.simutype (1,:) char {mustBeMember(nvp.simutype,{'adaptative','forward-euler','RK4'})} = 'adaptative'
    nvp.FSprecpow (1,1) double = 0
    nvp.enforce_real logical = []
end

if isempty(nvp.enforce_real)
    try
        AA = PhasorArray(Aph);
        if isreal(AA)
            nvp.enforce_real=true;
            warning('TransitionMatrixOverT:appearsRealValued', ['Matrix A appears to be ' ...
                'real-valued, so the integration is run over the reals. Pass ' ...
                'enforce_real=false to force the complex path.'])
        else
            nvp.enforce_real=false;
        end
    catch
        nvp.enforce_real=false;
    end
end


[nx,ny]=size(Aph,[1 2]);
nh = (size(Aph,3)-1)/2;

if nx~=ny
    error('A(t) must be square matrix')
end

Phi0T = zeros(ny,nx);
p = nextpow2(max(nh*4,2^nvp.FSprecpow));
dt_sim=T/2^p;
odeRelTol = 1e-10;   % the tolerance the adaptive branch asks for, below
for nxi = 1:nx
    ei=zeros(nx,1);
    ei(nxi)=1;
    switch nvp.simutype
        case 'adaptative'
            opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',dt_sim,'Stats','off');
            [y]=hmq_sim(Aph,0:dt_sim:T,ei,T,"odeOpts", opts,"plot", false,"isRealValued", nvp.enforce_real);
            Phi0T(:,nxi)=y(:,end);
        case 'forward-euler'
            [y]=hmq_sim(Aph,0:dt_sim:T,ei,T,"plot", false,"solver", 'forward-euler',"FSprecpow", nvp.FSprecpow,"isRealValued", nvp.enforce_real);
            Phi0T(:,nxi)=y(:,end);
        case 'RK4'
            [y]=hmq_sim(Aph,0:dt_sim:T,ei,T,"plot", false,"solver", 'RK4',"FSprecpow", nvp.FSprecpow,"isRealValued", nvp.enforce_real);
            Phi0T(:,nxi)=y(:,end);
    end
end

% The fast Floquet exponents are lost long before the integration is at fault.
% Phi(T) carries the multipliers exp(mu_i*T); once the smallest sits below the
% roundoff of the largest, no solver and no tolerance can recover it. Measured
% on a 2x2 with prescribed exponents: at ratio 1.9e-03 the exponent error is
% 1.5e-08, at 3.5e-06 it is 8.3e-06, at 1.1e-11 it is 1.5e-02 -- the same for
% ode15s and ode78. Use HmqNEig instead, which diagonalises T_h(A)-N and never
% exponentiates: exact to 1.4e-13 on the same cases, including mu = -200.
% The estimate below is calibrated on the adaptive branch; the fixed-step
% solvers are no better, so it stays a lower bound on the error for those.
Ewarn = eig(Phi0T);
ratio = min(abs(Ewarn)) / max(abs(Ewarn));
muFloor = odeRelTol / max(ratio * T, realmin);
if muFloor > 1e-6
    % muFloor triggers the warning but is not printed: it tracks the right
    % magnitude only loosely (optimistic by 50x at ratio 3.5e-06, pessimistic
    % by 1e4 at 1.2e-14), so quoting it as an accuracy would mislead. The
    % span is the hard fact.
    warning('TransitionMatrixOverT:multiplierUnderflow', ...
        ['Floquet multipliers span %.1e, so the fastest exponents are below the ' ...
         'roundoff of the slowest and cannot be recovered from Phi(T) at any ' ...
         'tolerance. Prefer HmqNEig(A,h,T,''fundamental''), which diagonalises ' ...
         'T_h(A)-N and is unaffected.'], ratio);
end

if nargout>1
    E=Ewarn;
end
if nargout>2
    [V,D,W]=eig(Phi0T);

end