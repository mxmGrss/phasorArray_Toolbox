function [pA] = TB2array(ATB,nx,size_type)
    %TB2ARRAY Convert a Toeplitz Block matrix to a 3D array by extracting anti-diagonal elements.
    %
    %   pA = TB2array(ATB, nx, size_type) converts a Toeplitz Block matrix `ATB` into a 3D array 
    %   by extracting the elements from the anti-diagonal of the Toeplitz blocks. The function uses 
    %   the provided matrix dimensions to organize the matrix into a 3D format with specified truncation.
    %
    %   Inputs:
    %     ATB         - (matrix) The Toeplitz Block matrix to be converted. Size: [nxd, nyd].
    %     nx          - (scalar) The first dimension of the time-dependent matrix A(t).
    %     size_type   - (string, optional) Specifies which dimension of A(t) is provided:
    %                   - 'n1' (default) uses the first dimension of A(t).
    %                   - 'n2' uses the second dimension of A(t).
    %                   - 'h' specifies that ATB is the Toeplitz matrix of order h of A(t).
    %                   - '2h+1' specifies that ATB represents the full Toeplitz of A(t).
    %
    %   Outputs:
    %     pA          - (3D array) The extracted anti-diagonal phasors from the Toeplitz Block matrix,
    %                   reshaped into a 3D array of size [nx, ny, 2*h+1].
    %
    %   Description:
    %     The function takes a Toeplitz Block matrix `ATB` and reshapes it into a 3D array `pA`. The 
    %     elements are extracted from the anti-diagonal of the Toeplitz blocks, where each phasor 
    %     corresponds to a particular harmonic order. The `size_type` argument allows the user to specify
    %     which dimension of the matrix `A(t)` is given. The function assumes that the matrix dimensions
    %     follow the conventions defined by the `size_type` argument.
    %
    %   Example Usage:
    %     % Convert a Toeplitz Block matrix into a 3D phasor array
    %     pA = TB2array(ATB, 5, 'n1');
    %
    %   See also: kron, diag, reshape, error.

arguments
    ATB
    nx
    size_type string {mustBeMember(size_type,["n1","n2","h","2h+1"])} = "n1"
end

[nxd,nyd]=size(ATB);
switch size_type
    case "n1"
        m=nxd/nx; %number of truncature, = 2*h+1
        ny=nyd/m; %deduced second dim
    case "n2"
        m=nyd/nx; %number of truncature = 2*h+1
        nx=nxd/m; %deduced first dim
    case "h"
        m=2*nx+1;
        ny=nyd/m;
        nx=nxd/m;
    case "2h+1"
        m=nx;
        ny=nyd/m;
        nx=nxd/m;
    otherwise
        error('error in TB2array dim spec')
end



pA=zeros(nx,ny,2*m-1)*1i;
for nxi=1:nx
    for nyi=1:ny
        Eij=zeros(nx,ny);
        Eij(nxi,nyi)=1;
        dum=kron(Eij,ones(m));
        D=reshape(ATB(boolean(dum)),m,m);
        ph_aij1=diag(fliplr(D));
        ph_aij2=diag(fliplr(D),1);
        ph_aij=zeros(2*m-1,1)*1i;
        ph_aij(1:2:end)=ph_aij1;
        ph_aij(2:2:(end-1))=ph_aij2;
        pA(nxi,nyi,:)=ph_aij;
    end
end


end