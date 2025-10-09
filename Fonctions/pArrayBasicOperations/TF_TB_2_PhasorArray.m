function [PA] = TF_TB_2_PhasorArray(TF_TB,n1,n2)
    %TF_TB_2_PHASORARRAY Convert Fourier Transform (TB format) to PhasorArray form.
    %
    %   PA = TF_TB_2_PhasorArray(TF_TB, n1, n2) converts a vector `TF_TB`, which represents
    %   the Fourier decomposition in a Toeplitz block (TB) format, into a PhasorArray. The
    %   function reshapes the input vector and permutes it to match the structure of the 
    %   corresponding periodic matrix with dimensions `n1` and `n2`.
    %
    %   Inputs:
    %     TF_TB   - (vector) A vector containing the Fourier coefficients in Toeplitz block
    %               (TB) format, representing the vectorization of a periodic matrix.
    %     n1      - (scalar, optional) The first dimension of the matrix. Default: 1.
    %     n2      - (scalar, optional) The second dimension of the matrix. Default: 1.
    %
    %   Outputs:
    %     PA      - (PhasorArray) The converted PhasorArray representing the Fourier 
    %               decomposition of the periodic matrix.
    %
    %   Description:
    %     This function takes a vector `TF_TB` containing Fourier coefficients and
    %     reshapes it into a 3D array of size `n1 x n2 x (2h + 1)` corresponding to the
    %     periodic matrix in question. The transformation relies on the assumption that 
    %     the input vector `TF_TB` represents the Fourier decomposition in Toeplitz block 
    %     format. The resulting 3D array is then converted into a PhasorArray.
    %
    %   Example Usage:
    %     % Convert Fourier decomposition in TB format to PhasorArray form
    %     PA = TF_TB_2_PhasorArray(TF_TB, 3, 3);
    %
    %   See also: PhasorArray, reshape, permute, error.
arguments
    TF_TB
    n1=1
    n2=1
end

HB=reshape(TF_TB,[],n1,n2);
PA=permute(HB,[2 3 1]);
PA = PhasorArray(PA);
end

