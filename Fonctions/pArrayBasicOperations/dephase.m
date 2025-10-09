function [B] = dephase(A,angle)
%DEPHASE change phase of periodic matrix defined by its phasor A, by an
%angle of angle. basically every phasor is multiplied by exp jk angle
h=(size(A,3)-1)/2;
base(1,1,:)=exp((-h:h)*1i*angle);
try
    B=pagemtimes(base,value(A));
catch 
    for ii=1:h
        B(:,:,h+1+ii)=A(:,:,h+1+ii)*base(1,1,h+1+ii);
        B(:,:,h+1-ii)=A(:,:,h+1-ii)*base(1,1,h+1-ii);
    end
        ii=0;
        B(:,:,h+1-ii)=A(:,:,h+1+ii)*base(1,1,h+1-ii);
end
if isa(A,'PhasorArray')
    B=PhasorArray(B,'reduce',0);
end
end

