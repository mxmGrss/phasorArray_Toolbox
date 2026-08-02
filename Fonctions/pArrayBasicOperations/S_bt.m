function [out] = S_bt(h,phi,n)
%S_bt produce a Block toeplitz dephasing matrix until the order h
arguments
    h (1,1) double
    phi
    n (1,1) double
end
if ~iscolumn(phi)
    phi=phi.';
end

temp=cell(2*h+1,1);
for h_i=1:numel(-h:h)
temp{h_i} = (exp((-h+h_i-1)*1i*phi))*eye(n);
end
out= blkdiag(temp{:});
end
