function [out] = S_tb(h,phi,n)
%S_tb produce a Toeplitz Block dephasing matrix until the order h
arguments
    h (1,1) double
    phi
    n (1,1) double
end
temp=cell(numel(phi),1);
for phi_i=1:numel(phi)
temp{phi_i} = kron(eye(n),diag(exp((-h:h)*1i*phi(phi_i))));
end
out= cat(1,temp{:});
end
