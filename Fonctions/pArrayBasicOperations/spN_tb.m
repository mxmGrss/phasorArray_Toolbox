function [N,Nw] = spN_tb(n,nh,T)
%spN_tb(n,nh,T) output the derivation matrix N = eye(n) kron diag(jkw) associated to
%the TB harmonic matrix, n is the size of the state, nh the order of
%truncature, T the period
% 
arguments
    n
    nh
    T=1
end

if ~isscalar(n) %a matrix or 3D array of phasor is provided
    n=size(n,1); %we match the first dim of the 3D array phasor
end


k=(-nh:nh)';
w=2*pi/T;
Nw=sparse(1:2*nh+1,1:2*nh+1,1i*k*w);
N = kron(speye(n),Nw) ;
end

