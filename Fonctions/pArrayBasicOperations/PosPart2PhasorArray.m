function [Aph] = PosPart2PhasorArray(A0,Apos)
%POSPART2PHASORARRAY expect a matrix and a 3D Array, the first arg is the 0
%phasor, the second are the positive phasors, they are conjugated to
%determine the negatives ones
arguments
    A0 
    Apos 
end

if isa(A0,'PhasorArray')
    A0=A0.value;
end
if isa(Apos,'PhasorArray')
    Apos=Apos.value;
end


n1 =size(Apos,1);
n2 =size(Apos,2);
n3 =size(Apos,3);
[m1,m2]=size(A0);

%Check dimensionality
if n3==1 %case Apos is not a 3D array, can be the case for convenience if A(t) is scalar or column/row vector, usign the other dimension to store phasor
    if n1 == m1 && n2 == m2
        %nothing to see there
    elseif isvector(A0)
        if isscalar(A0)
            if isrow(Apos)
                Apos=Apos.';
            elseif iscolumn(Apos)
            else
                warning('In the scalar case, Apos needs to be a vector')
                return
            end
            Apos=permute(Apos,[2,3,1]);
        elseif iscolumn(A0) || n1==m1
            Apos=permute(Apos,[1,3,2]);
        elseif isrow(A0) || n2==m2
            Apos=permute(Apos,[3,2,1]);
        else
            warning('Erreur de dimensions')
            return
        end
    end
end



B=real(flip(Apos,3))-1i*imag(flip(Apos,3));
Aph=cat(3,B,A0,Apos);
Aph=PhasorArray(Aph,reduce=0);
end

