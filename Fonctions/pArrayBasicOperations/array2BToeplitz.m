function [Tmat, N] = array2BToeplitz(Aph, m, nvp)
    %ARRAY2BTOEPLITZ Convert a 3D array into a Block-Toeplitz (BT) matrix.
    %
    %   Usage:
    %       [Tmat, N] = array2BToeplitz(Aph)
    %       [Tmat, N] = array2BToeplitz(Aph, m)
    %       [Tmat, N] = array2BToeplitz(Aph, [h1, h2])
    %
    %   Inputs:
    %       Aph   - Input data. 3D double array or PhasorArray object.
    %       m     - (Optional) Order of the output block matrix.
    %               - If scalar h: Square block matrix of size (2*h+1) x (2*h+1).
    %               - If vector [h1, h2]: Rectangular matrix of size (2*h1+1) x (2*h2+1).
    %               - Default: Current order of Aph.
    %       nvp  - Name-Value arguments.
    %
    %   Outputs:
    %       Tmat  - Block-Toeplitz matrix.
    %       N     - Differentiation matrix (only for square case).

    arguments
        Aph
        m = []
        nvp.method {mustBeMember(nvp.method, {'cell2mat', 'cat'})} = 'cell2mat';
    end

    % BT layout: block (i,j) holds harmonic k = i - j + (h2 - h1), i.e.
    % T(A)_{k,m} = A_{k-m} with harmonics ascending (-h..+h) as in
    % array2TBlocks. BT and TB then differ by a pure perfect shuffle.
    [Aph, h1, h2, nvp.method] = prepToeplitzInput(Aph, m, nvp.method);
    n1 = size(Aph, 1);
    n2 = size(Aph, 2);
    h_req = h1 + h2;

    switch nvp.method
        case 'cell2mat'
            Aph2 = mat2cell(Aph, n1, n2, ones(2 * h_req + 1, 1));
            % Formula: idx = i - j + 2*h2 + 1  (same as array2TBlocks)
            col = (1 : 2 * h1 + 1) + 2 * h2;
            row = 2 - (1 : 2 * h2 + 1) + 2 * h2;
            Toe = toeplitz(col, row);
            Tmat = cell2mat(Aph2(Toe));

        case 'cat'
            Tmat = [];
            for xi = 1:2*h1+1
                line = [];
                for yi = 1:2*h2+1
                    k = xi - yi + h2 - h1;
                    line = [line, Aph(:, :, k + h_req + 1)];
                end
                Tmat = [Tmat; line];
            end
    end

    if nargout == 2
        if h1 == h2
            N = N_bt(n1, h1, 1);
        else
            warning('array2BToeplitz:nonSquareN', 'N matrix is only defined for square Block-Toeplitz truncations (h1=%d ~= h2=%d); returning empty N.', h1, h2)
            N = [];
        end
    end
end