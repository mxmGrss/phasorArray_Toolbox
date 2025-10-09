function [N,Nw] = spNBT(n,nh,T)
%NBT output the derivation matrix N = diag(jkw) kron eye(n)  associated to
%the BT harmonic matrix
% 
arguments
    n
    nh
    T=1
end


if ndims(n)>1 %a matrix or 3D array of phasor is provided
    n=size(n,1); %we match the first dim of the 3D array phasor
end


k=(-nh:nh)';
w=2*pi/T;
Nw=sparse(1:2*nh+1,1:2*nh+1,1i*k*w);
N = kron(Nw,speye(n)) ;
end

