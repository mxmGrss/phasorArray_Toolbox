function [OutM,N] = array2TBlocks(Aph,m,varg)
    %ARRAY2TBLOCKS Convert a 3D array into a Toeplitz-blocks matrix
    %   This function takes a 3D array representing periodic matrices and converts
    %   it into a matrix of Toeplitz blocks. Each block corresponds to the Toeplitz
    %   representation of each coefficient of the array.
    %
    %   The function allows for padding or truncating the input array to match the 
    %   specified number of blocks (m), ensuring the output has a length of (nh+1) blocks.
    %   The method used for matrix construction can be chosen using the 'method' option.
    %   
    %   Usage:
    %   [OutM, N] = array2TBlocks(Aph, m, varg)
    %
    %   Inputs:
    %       Aph       - 3D array [A_{-nh}, ..., A_0, ..., A_{nh}] to convert into a Toeplitz-blocks matrix.
    %       m         - Optional padding or truncation to define the number of blocks (default: nh*2).
    %       varg      - Optional name-value arguments for method selection.
    %
    %   Outputs:
    %       OutM      - The resulting matrix of Toeplitz blocks.
    %       N         - Differenciation matrix associated with TB structure and compatible size.
    %
    %   Options:
    %       varg.method  - Method used to create the Toeplitz blocks matrix ('cell2mat' or 'cat'). Default is 'cell2mat'.
    %
    %   Notes:
    %       - The 'cell2mat' method is more efficient and uses `cell2mat` to construct the blocks.
    %       - The 'cat' method uses concatenation with `toeplitz` to build the blocks, making it compatible with `sdpvar`.
    %
    %   See also: toeplitz, mat2cell, sdpvar
    
arguments
    Aph
    m=[]
    varg.method='cell2mat';
end
if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end

if isa(Aph,"ndsdpvar") || isa(Aph,"sdpvar") || isa(Aph,"sym")
    varg.method='cat';
end

if ismatrix(Aph)
    [n1,n2]=size(Aph);
    nhlenbis=1;
else
    [n1,n2,nhlenbis]=size(Aph);
end
    nh=(nhlenbis-1)/2;

    
if isempty(m)
    m=nh*2;

else 
    if nh>m
        Aph=Aph(:,:,(nh+1+(-m:m)));
        nh=m;
    elseif nh<m
%         Aph=padarray(Aph,[0 0 (m-nh)]);
%             nh=m
    end

end

switch varg.method
    case 'cell2mat'
% OutM=zeros(n1*(m+1),n2*(m+1));
    c=cell(n1,n2);
        for xi = 1:n1
            for yi=1:n2
                ui=[zeros(m-nh,1) ; squeeze(Aph(xi,yi,:)) ; zeros(m-nh,1)];
                TAij = toeplitz(ui((m+1):end), ui((m+1):-1:1) );
%                 TAij
                c{xi,yi}=TAij;
%                 OutM((xi-1)*(m+1)+(1:(m+1)),(yi-1)*(m+1)+(1:(m+1)))=TAij;
            end
        end
        % c
        OutM=cell2mat(c);

    case 'cat'
%         OutM=zeros(n1*(m+1),n2*(m+1));
        OutM=[];
        % c=cell(n1,n2);
        for xi = 1:n1
            MT1=[];
            for yi=1:n2
                try
                ui=[zeros(m-nh,1) ; squeeze(Aph(xi,yi,:)) ; zeros(m-nh,1)];
                catch e
                    xi
                    yi
                    Aph
                    Aph(xi,yi,:)
                    error(e)
                end
                ui1=ui((m+1):end);
                ui2=ui((m+1):-1:1);
                TAij = toeplitz(ui1, ui2 );
%                 TAij
        %         c{xi,yi}=TAij;
                MT1=[MT1, TAij];
%                 OutM((xi-1)*(m+1)+(1:(m+1)),(yi-1)*(m+1)+(1:(m+1)))=TAij;
            end
            OutM=[OutM; MT1];
        end
end


if nargout==2
    N = NTB(n1,m/2,1);
end

end