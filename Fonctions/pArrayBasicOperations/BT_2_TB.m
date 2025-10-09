function [outM] = BT_2_TB(Ahm,nbloc)
    %BT_2_TB Convert from block Toeplitz form to Toeplitz-block form
    %   This function converts a harmonic matrix from its block Toeplitz form, 
    %   where each block is the k-nth phasor of A(t), to its Toeplitz-block form, 
    %   where each block represents the Toeplitz matrix of the phasor of each coefficient of A(t).
    %
    %   In the affine case where Ahm = kron(D, Ad), this function transforms it to 
    %   kron(Ad, D) format.
    %
    %   Usage:
    %   [outM] = BT_2_TB(Ahm, nbloc)
    %
    %   Inputs:
    %       Ahm   - Input matrix in block Toeplitz form.
    %       nbloc - The size of the blocks used in the Toeplitz matrix.
    %
    %   Outputs:
    %       outM  - The resulting matrix in Toeplitz-block form.
    %
    %   See also: kron, shuffle_matrix
    
[ni,nj]=size(Ahm);
nia=ni/nbloc;
nja=nj/nbloc ;
nid=nbloc;
njd=nbloc;
outM=shuffle_matrix(nid,nia)*Ahm*shuffle_matrix(njd,nja)';
end