function [Tmat,N] = array2BToepliz(Aph,m,varg)
    %ARRAY2BTOEPLIZ Convert a 3D array into a Block-Toeplitz matrix
    %   This function takes a 3D array representing a periodic matrix, where each slice
    %   corresponds to a phasor at a specific frequency, and converts it into a Block-Toeplitz matrix.
    %   The input array should be in the form [A_{-nh}, ..., A_0, ..., A_{nh}], and the resulting
    %   Block-Toeplitz matrix has the following structure:
    %
    %   Tmat = [A_0  A_1  A_2 ... A_{nh} ;
    %           A_{-1}  A_0  A_1 ... A_{nh-1} ;
    %           A_{-2}  A_{-1} A_0 ... A_{nh-2} ;
    %           ... ;
    %           A_{-nh}  A_{-nh+1} ... A_0]
    %
    %   The function can pad or truncate the input array depending on the value of `m` to ensure that 
    %   the output matrix has the desired number of blocks.
    %
    %   Usage:
    %   [Tmat, N] = array2BToepliz(Aph, m, varg)
    %
    %   Inputs:
    %       Aph       - 3D array representing the phasors [A_{-nh}, ..., A_0, ..., A_{nh}].
    %       m         - Optional padding or truncation for the number of blocks.
    %       varg      - Optional name-value arguments for method selection.
    %
    %   Outputs:
    %       Tmat      - Block-Toeplitz matrix formed from the input phasor array.
    %       N         - Differenciation matrix associated with BT structure and compatible size.
    %
    %   Options:
    %       varg.method  - Method used to create the Toeplitz matrix ('cell2mat' or 'cat'). Default is 'cell2mat'.
    %
    %   Notes:
    %       - The method 'cell2mat' constructs the matrix using `mat2cell` and `toeplitz` functions.
    %       - The method 'cat' uses the `kron` and `diag` functions to construct the matrix with blocks.
    %
    %   See also: toeplitz, mat2cell, kron, diag
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


if nargin == 1 || isempty(m)
    [n1,n2,nhlenbis]=size(Aph);
    nh=(nhlenbis-1)/2;
    m=[2*nh 2*nh];
end
    if numel(m)==1
        m=[m m];
    end

    [n1,n2,nhlenbis]=size(Aph);
    nh=(nhlenbis-1)/2;

    if nh>max(m)
        Aph=Aph(:,:,(nh+1+(-max(m):max(m))));
    elseif nh<max(m)
%             Aphd=ndsdpvar(size(Aphi,1),size(Aphi,2),max(m));
%             Aphd(:,:,(end+1)/2+(-nh:nh))=Aph;
%             Aphd(:,:,1:(end+1)/2+(-nh-1))=0;
%             Aphd(:,:,(end+1)/2+(+nh+1):end)=0;
%             Aph=Aphd;
            Aphd=cat(3,zeros(n1,n2,(max(m)-nh)),Aph,zeros(n1,n2,(max(m)-nh)));
            Aph=Aphd;
    
end


switch varg.method
    case 'cell2mat'
    Aph2=mat2cell(Aph,n1,n2,ones(2*max(m)+1,1));
    Toe=toeplitz((max(m)+1):(max(m)+m(1)+1),(max(m)+1):-1:(max(m)+1-m(2)));
    Tmat = cell2mat(Aph2(Toe));
    otherwise
        Tmat=zeros(n1*(max(m)+1),n2*(max(m)+1));
        p=max(m)+1;
        for kk=0:max(m)
            Imk=diag(ones(p-kk,1),kk);
            Ik=diag(ones(p-kk,1),-kk);
            Tmat=Tmat+kron(Ik,Aph(:,:,max(m)+1+kk));
            if kk>0
            Tmat=Tmat+kron(Imk,Aph(:,:,max(m)+1-kk));
            end
        end
end


if nargout==2
    N = N_bt(n1,m/2,1);
end
end