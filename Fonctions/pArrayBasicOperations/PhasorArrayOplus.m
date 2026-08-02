function [AoBph] = PhasorArrayOplus(Aph,Bph)
%PHASORARRAYOPLUS realise the kronecker sum of A(t) and B(t)
%   A(t) oplus B(t) = A(t) kron I_nB + I_nA kron B(t)
%   Matrix must be square. 
%   3D array of phasors are accepted
%   if one arg is empty, it is replaced with a zero square matrix of same
%   dimension as the other provided argument.

arguments
    Aph=[]
    Bph=[]
end

Aph = pvalue(Aph);
Bph = pvalue(Bph);


funcA =  (size(Aph, 1) ~= size(Aph, 2));
funcB =  (size(Bph, 1) ~= size(Bph, 2));

if funcB || funcA
    error('PhasorArrayOplus:invalidInput', 'Arguments must be either empty ([]) or 3-D arrays of square matrices.')
end
if isempty(Aph)
    Aph=Bph*0;
elseif isempty(Bph)
    Bph=Aph*0;
elseif isempty(Aph) && isempty(Bph)
    error('PhasorArrayOplus:noInput', 'At least one non-empty matrix argument is required.')
else 
end

na=size(Aph,1);
nb=size(Bph,1);



hA=nHarm(Aph);
hB=nHarm(Bph);
h=max(hA,hB);
Ia=speye(na);
Ib=speye(nb);



if hA<hB
    Aph=phasorPad(Aph,[0 0 hB-hA]);
elseif hB<hA
    Bph=phasorPad(Bph,[0 0 hA-hB]);
end

%More obvious but less efficient implementation
% AoBph=zeros(na*nb,na*nb,size(Aph,3));
% t2=tic;
for ii=1:size(Aph,3)
    AoBph(:,:,ii)=full(kron(Aph(:,:,ii),Ib)+kron(Ia,Bph(:,:,ii)));
end

end

