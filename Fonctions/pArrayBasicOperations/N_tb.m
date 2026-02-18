function [N,Nw] = N_tb(n,nh,T,varg)
%N_tb(n,nh,T) output the derivation matrix N = eye(n) kron diag(jkw) associated to
%the TB harmonic matrix, n is the size of the state, nh the order of
%truncature, T the period
% 
arguments
    n
    nh
    T=[]
    varg.omega = []
end

if isempty(varg.omega)
    if isempty(T)
        T = 2 * pi; % Default period if not provided
        w=2*pi/T;
    elseif isinf(T)
        w=0;
    else
        w=2*pi/T;
    end
else
    w = varg.omega; % Use provided omega if available
    if ~isempty(T)
        error("In PhasoArray/N_tb Function : Don't provide T and omega at the same time")
    end
end

if isnan(w)
    error("In PhasoArray/N_tb Function : omega cannot be NaN ")
end

if isphasor(n)
    n=size(n,1);
end
if ~isscalar(n) %a matrix or 3D array of phasor is provided
    n=size(n,1); %we match the first dim of the 3D array phasor
end



k=(-nh:nh)';
Nw=diag(1i*k*w);
N = kron(eye(n),Nw) ;
end

