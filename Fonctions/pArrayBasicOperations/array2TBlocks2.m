function [OutM] = array2TBlocks2(Aph,m)
    %ARRAY2TBLOCKS2 Convert a 3D array to a Toeplitz-blocks matrix
    %   This function takes a 3D array [A_{-nh}, ..., A_0, ..., A_{nh}] and constructs
    %   a matrix made of Toeplitz blocks, each representing the Toeplitz matrix
    %   of each coefficient. If specified, the 3D array is padded or truncated with zeros
    %   so that the output consists of (nh+1) blocks.
    %
    %   Usage:
    %   [OutM] = array2TBlocks2(Aph, m)
    %
    %   Inputs:
    %       Aph - 3D array [A_{-nh}, ..., A_0, ..., A_{nh}] to convert into a Toeplitz-blocks matrix.
    %       m   - Optional padding or truncation to define the number of blocks (default: 2 * nh).
    %
    %   Outputs:
    %       OutM - The resulting matrix of Toeplitz blocks.
    %
    %   Notes:
    %       - The method builds a matrix of Toeplitz blocks, with each block representing
    %         a Toeplitz matrix for the corresponding coefficient.
    %       - Padding or truncating the 3D array is handled automatically to fit the desired
    %         number of blocks (nh+1).
    %
    %   See also: toeplitz, padarray

if isa(Aph,'PhasorArray')
    Aph=Aph.Value;
end


if nargin == 1
    [n1,n2,nhlenbis]=size(Aph);
    nh=(nhlenbis-1)/4;
    m=2*nh;

OutM=zeros(n1*(m+1),n2*(m+1));
for xi = 1:n1
    for yi=1:n2
        ui = squeeze(Aph(xi,yi,:));
        TAij = toeplitz(ui((m+1):end), ui((m+1):-1:1) );
        OutM((xi-1)*(m+1)+(1:(m+1)),(yi-1)*(m+1)+(1:(m+1)))=TAij;
    end
end





elseif nargin==2



    [n1,n2,nhlenbis]=size(Aph);
    nh=(nhlenbis-1)/2;

        dA=permute(Aph,[3,1,2]);
        colA2=reshape(dA,[],n1*n2);

    OutM=zeros(n1*(m+1),n2*(m+1));
    for xi = 1:n1
        for yi=1:n2
        
        ui = colA2(:,(yi-1)*n1+xi);
    
        if nh>m
            ui=ui((nh+1+(-m:m)));
        elseif nh<m
            ui=padarray(ui,(m-nh));
        end
            TAij = toeplitz(ui((m+1):end), ui((m+1):-1:1) );
            OutM((xi-1)*(m+1)+(1:(m+1)),(yi-1)*(m+1)+(1:(m+1)))=TAij;
        end
    end
end
end