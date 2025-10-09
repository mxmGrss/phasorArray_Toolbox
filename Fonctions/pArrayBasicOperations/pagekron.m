function [Cph] = pagekron(Aph,Bph)
%TensorKron realise slicewise kron product between 3D array A and B
%   ie C(:,:,i) = kron(Aph(:,:,i),Bph(:,:,i))
%   if A or B is a matrix (size(A,3)==1), matrix is duplicated to match
%   other arg size. 
%   
if isa(Bph,'PhasorArray')
    Bph=Bph.Value;
end
if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end
sa=size(Aph);
sb=size(Bph);

if numel(sa)==2
    Aph=repmat(Aph,1,1,sb(3));
elseif numel(sb)==2
    Bph=repmat(Bph,1,1,sa(3));
end


z = bsxfun(@times, permute(Aph, [4 1 5 2 3]), permute(Bph, [1 4 2 5 3]));  %// step 1
Cph = reshape(z, size(Aph,1)*size(Bph,1), size(Aph,2)*size(Bph,2), size(Aph,3));


end