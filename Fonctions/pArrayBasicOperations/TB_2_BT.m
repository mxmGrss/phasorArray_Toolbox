function [outM] = TB_2_BT(chm,nbloc)
    %TB_2_BT Convert a harmonic matrix from Toeplitz-block form to block Toeplitz form.
    %
    %   outM = TB_2_BT(chm, nbloc) converts a harmonic matrix from its Toeplitz-block form, 
    %   where each block is the Toeplitz matrix of the phasor for each coefficient of the matrix A(t),
    %   into its block Toeplitz form, where each block is the k-th phasor of A(t).
    %
    %   In the affine case, `chm = kron(Ad, D)` and the function aims to convert it into `kron(D, Ad)`.
    %
    %   Inputs:
    %     chm      - (matrix) The harmonic matrix in Toeplitz-block form, where each block is the 
    %                Toeplitz matrix of the phasor for each coefficient of A(t).
    %     nbloc    - (scalar) The number of blocks used to reshape the matrix.
    %
    %   Outputs:
    %     outM     - (matrix) The transformed matrix in block Toeplitz form, where each block is the 
    %                k-th phasor of A(t).
    %
    %   Description:
    %     The function performs matrix manipulation by reshaping and applying Kronecker products to 
    %     transform the input harmonic matrix `chm` from a Toeplitz-block form into the block Toeplitz 
    %     form. The `nbloc` parameter defines the number of blocks, which is used to partition the input 
    %     matrix `chm` into smaller blocks, then permute the matrix appropriately to obtain the desired 
    %     block Toeplitz structure.
    %
    %   Example Usage:
    %     % Convert harmonic matrix from Toeplitz-block form to block Toeplitz form
    %     outM = TB_2_BT(chm, 3);
    %
    %   See also: kron, shuffle_matrix.
[ni,nj]=size(chm);
nia=ni/nbloc;
nja=nj/nbloc ;
nid=nbloc;
njd=nbloc;
outM=shuffle_matrix(nja,njd)*chm*shuffle_matrix(nia,nid)';
end