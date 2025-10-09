function [out] = SBT(h,phi)
%STB produce a Block toeplitz dephasing matrix until the order h
if ~iscolumn(phi)
    phi=phi.';
end

temp=cell(2*h+1,1);
for h_i=1:numel(-h:h)
temp{h_i} = (exp((-h+h_i-1)*1i*phi));
end
out= blkdiag(temp{:});
end
