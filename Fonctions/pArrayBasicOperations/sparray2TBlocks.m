function [OutM,N] = sparray2TBlocks(Aph,m,varg)
%array2TBlocks Summary : 3D array to Toeplitz-blocs matrix SPARSE VERSION
%   recoit une liste de matrice (3D array) [A-nh1 ... A0 ... Anh1] et
%   construit la matrice faites de Blocs Toeplitz, chacun représentant la
%   toeplitz de chacun des coeefficients. 
%   if specified, nhcible pad / truncate the 3D array with zeros so that output is of len (nh+1) blocks. 
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

elseif nargin>1
    [n1,n2,nhlenbis]=size(Aph);
    nh=(nhlenbis-1)/2;
    if nh>m
        Aph=Aph(:,:,(nh+1+(-m:m)));
        nh=m;
    elseif nh<m
%         Aph=padarray(Aph,[0 0 (m-nh)]);
    end

end

switch varg.method
    case 'cell2mat'
        % OutM=zeros(n1*(m+1),n2*(m+1));
        % OutM=spalloc(n1*(m+1),n2*(m+1),((m+1)*(m+1)-(m-nh)*(m-nh+1))*n1*n2);
        c=cell(n1,n2);
        % OutM=[];
        for xi = 1:n1
        %     MT1=[];
            for yi=1:n2
        %         ui=zeros(2*m+1,1);
        %         ui(m+1+(-nh:nh))=squeeze(Aph(xi,yi,:));
                ui = sparse(m+1+(-nh:nh),ones(1,2*nh+1),squeeze(Aph(xi,yi,:)),2*m+1,1);
        %         uip=sparse(1+(0:nh),1,squeeze(Aph(xi,yi,nh+1:end)),m+1,1);
        %         uim=sparse(1+(0:nh),1,squeeze(Aph(xi,yi,nh+1:-1:1)),m+1,1);
        %         TAij = toeplitz(uip, uim);
                TAij = toeplitz(ui((m+1):end), ui((m+1):-1:1));
                c{xi,yi}=TAij;
        %         OutM((xi-1)*(m+1)+(1:(m+1)),(yi-1)*(m+1)+(1:(m+1)))=TAij;
        %         MT1=[MT1 TAij];
            end
        %     OutM=[OutM; MT1];
        end
        OutM=cell2mat(c);


    case 'cat'
        OutM=[];
        for xi = 1:n1
            MT1=[];
            for yi=1:n2
                ui = sparse(m+1+(-nh:nh),ones(1,2*nh+1),squeeze(Aph(xi,yi,:)),2*m+1,1);
                TAij = toeplitz(ui((m+1):end), ui((m+1):-1:1));
                MT1=[MT1 TAij];
            end
            OutM=[OutM; MT1];
        end
end



if nargout==2
    N = NTB(n1,m/2,1);
end
end