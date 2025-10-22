function [out] = S_tb(h,phi)
%S_tb produce a Toeplitz Block dephasing matrix until the order h
temp=cell(numel(phi),1);
for phi_i=1:numel(phi)
temp{phi_i} = diag(exp((-h:h)*1i*phi(phi_i)));
end
out= cat(1,temp{:});
end