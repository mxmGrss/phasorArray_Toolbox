function [AoBph] = PhasorArrayOplus(Aph,Bph)
%OPLUS realise the kronecker sum of A(t) and B(t)
%   A(t) oplus B(t) = A(t) kron I_nB + I_nA kron B(t)
%   Matrix must be square. 
%   3D array of phasors are accepted
%   if one arg is empty, it is replaced with a zero square matrix of same
%   dimension as the other provided argument.

arguments
    Aph=[]
    Bph=[]
end

if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end

if isa(Bph,'PhasorArray')
    Bph=Bph.Value;
end


funcA =  (size(Aph, 1) ~= size(Aph, 2));
funcB =  (size(Bph, 1) ~= size(Bph, 2));

if funcB || funcA
    error('Matrix arguments must be either empty matrix ([]) or 3D array of square matrix')
end
if isempty(Aph)
    Aph=Bph*0;
elseif isempty(Bph)
    Bph=Aph*0;
elseif isempty(Aph) && isempty(Bph)
    error('Function called without argument, at least one matrix must be non empty)')
else 
end

na=size(Aph,1);
nb=size(Bph,1);



hA=(size(Aph,3)-1)/2;
hB=(size(Bph,3)-1)/2;
h=max(hA,hB);
Ia=speye(na);
Ib=speye(nb);



if hA<hB
    Aph=padarray(Aph,[0 0 hB-hA]);
elseif hB<hA
    Bph=padarray(Bph,[0 0 hA-hB]);
end

%More obvious but less efficient implementation
% AoBph=zeros(na*nb,na*nb,size(Aph,3));
% t2=tic;
for ii=1:size(Aph,3)
    AoBph(:,:,ii)=full(kron(Aph(:,:,ii),Ib)+kron(Ia,Bph(:,:,ii)));
end
% toc(t2)

% t3=tic;

%On average more efficient implementation
%Code proposed by Luis Mendo on https://stackoverflow.com/questions/28547585/kronecker-product-between-two-tensors
% x=Aph;
% y=repmat(Ib,1,1,2*h+1);
% z1 = bsxfun(@times, permute(x, [4 1 5 2 3]), permute(y, [1 4 2 5 3]));  %// step 1
% z1 = reshape(z1, size(x,1)*size(y,1), size(x,2)*size(y,2), size(x,3));
% 
% x=repmat(Ia,1,1,2*h+1);
% y=Bph;
% z2 = bsxfun(@times, permute(x, [4 1 5 2 3]), permute(y, [1 4 2 5 3]));  %// step 1
% z2 = reshape(z2, size(x,1)*size(y,1), size(x,2)*size(y,2), size(x,3));
% 
% AoBph=z1+z2;

% toc(t3)

end

