function [Sm] = shuffle_matrix(p,q)
    %SHUFFLE_MATRIX Construct a matrix for converting between block Toeplitz (BT) and Toeplitz-block (TB) forms.
    %
    %   Sm = shuffle_matrix(p, q) constructs a matrix used to convert between the block Toeplitz (BT)
    %   form and the Toeplitz-block (TB) form. Specifically, this function is used for reshaping matrices
    %   from one form to the other in the context of harmonic matrix transformations.
    %
    %   Inputs:
    %     p - (scalar) The number of rows in the resulting matrix.
    %     q - (scalar) The number of columns in the resulting matrix.
    %
    %   Outputs:
    %     Sm - (matrix) The resulting matrix that performs the conversion between the two forms.
    %          The matrix is constructed by selecting rows from the identity matrix.
    %
    %   Description:
    %     This function creates a matrix `Sm` by extracting rows from the identity matrix `Ir` in a 
    %     particular pattern, effectively creating a reshaping matrix that can convert a block Toeplitz 
    %     matrix into a Toeplitz-block matrix or vice versa. The matrix `Sm` is used to manipulate the 
    %     structure of matrices in the context of harmonic analysis.
    %
    %   Example Usage:
    %     % Construct a shuffle matrix with p=3 and q=2
    %     Sm = shuffle_matrix(3, 2);
    %
    %   See also: kron, eye, TB_2_BT, BT_2_TB.
r=p*q;
Sm=[];
Ir=eye(r);
for qi=1:q
    Sm=[Sm; Ir(qi:q:r,:)];
end