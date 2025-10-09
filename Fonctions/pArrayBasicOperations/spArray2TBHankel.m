function [HpJ,JHm,Hp,Hm] = spArray2TBHankel(Aph,m,varg)
%ARRAY2HANKELP Summary 3D array to Hankel * J positive Toeplitz-blocs matrix
%   recoit une liste de matrice (3D array) [A-nh1 ... A0 ... Anh1] et
%   construit la matrice faites de Blocs Toeplitz, chacun représentant la
%   toeplitz de chacun des coeefficients. 
%   if specified, nhcible pad / truncate the 3D array with zeros so that output is of len (nh+1) blocks. 
%cell2mat is sligthly more efficient than concatenation of blocks, but not
%compatible with sdp var, hence the argument method.

arguments
    Aph
    m=[]
    varg.method='cell2mat'
end

if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end

if isa(Aph,"ndsdpvar") || isa(Aph,"sdpvar")
    varg.method='cat';
end
if nargin == 1
    [n1,n2,nhlenbis]=size(Aph);
    nh=(nhlenbis-1)/2;
    m=nh;




varg.method

elseif nargin>1

    [n1,n2,nhlenbis]=size(Aph);
    nh=(nhlenbis-1)/2;
    if nh>2*m+1
        Aph=Aph(:,:,(nh+1+(-2*m-1:2*m+1)));
        nh=2*m+1;
        nhlenbis=2*nh+1;
    elseif nh<2*m+1
%         Aph=padarray(Aph,[0 0 (m-nh)]);
%             nh=m
    end

end

switch varg.method
    case 'cell2mat'
% OutM=zeros(n1*(m+1),n2*(m+1));
cp=cell(n1,n2);
cm=cell(n1,n2);
        for xi = 1:n1
            for yi=1:n2
                ui=sparse((2*m+2)+(-nh:nh),ones(nhlenbis,1),squeeze(Aph(xi,yi,:)),4*m+3,1);
%                 ui=sparse([zeros(2*m+1-nh,1) ; squeeze(Aph(xi,yi,:)) ; zeros(2*m+1-nh,1)]);
                TAijp = toeplitz(ui((3*m+3):end), ui((3*m+3):-1:(2*m+3)) );
%                 TAij
                cp{xi,yi}=TAijp;
                TAijm = toeplitz(ui((1*m+1):2*m+1), ui((1*m+1):-1:(1)) );
%                 TAij
                cm{xi,yi}=TAijm;
%                 OutM((xi-1)*(m+1)+(1:(m+1)),(yi-1)*(m+1)+(1:(m+1)))=TAij;
            end
        end
        % c
        HpJ=cell2mat(cp);
        JHm=cell2mat(cm);

    case 'cat'
        HpJ=[];
        JHm=[];
        for xi = 1:n1
            HpJ1=[];
            JHm1=[];
            for yi=1:n2
                ui=sparse((2*m+2)+(-nh:nh),ones(nhlenbis,1),squeeze(Aph(xi,yi,:)),4*m+3,1);
                TAijp = toeplitz(ui((3*m+3):end), ui((3*m+3):-1:(2*m+3)) );
                TAijm = toeplitz(ui((1*m+1):2*m+1), ui((1*m+1):-1:(1)) );
                HpJ1=[HpJ1, TAijp];
                JHm1=[JHm1, TAijm];
            end
            HpJ=[HpJ; HpJ1];
            JHm=[JHm; JHm1];
        end
end

[II,JJ]=meshgrid(1:size(HpJ,1),1:size(HpJ,2));
HpJ=sparse(II(:),JJ(:),HpJ(:));

[II,JJ]=meshgrid(1:size(JHm,1),1:size(JHm,2));
JHm=sparse(II(:),JJ(:),JHm(:));

J=flip(speye(m+1));
Jp=kron(speye(n2),J);
Jm=kron(speye(n1),J);
Hp=HpJ*Jp;
Hm=Jm*JHm;


end