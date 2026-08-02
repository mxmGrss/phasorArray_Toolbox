function [OutM, N] = array2TBlocks(Aph, m, nvp)
    %ARRAY2TBLOCKS Convert a 3D array into a Toeplitz-blocks (TB) matrix.
    %
    %   Usage:
    %       [OutM, N] = array2TBlocks(Aph)
    %       [OutM, N] = array2TBlocks(Aph, m)
    %       [OutM, N] = array2TBlocks(Aph, [h1, h2])
    %
    %   Inputs:
    %       Aph   - Input data. 3D double array or PhasorArray object.
    %       m     - (Optional) Order of the output block matrix.
    %               - If scalar h: Square block matrix of size (2*h+1) x (2*h+1).
    %               - If vector [h1, h2]: Rectangular matrix of size (2*h1+1) x (2*h2+1).
    %               - Default: Current order of Aph.
    %       nvp  - Name-Value arguments:
    %               * 'method': 'cell2mat' (fast) or 'cat' (SDP/SYM compatible).
    %
    %   Outputs:
    %       OutM  - Toeplitz-Blocks matrix.
    %       N     - Differentiation matrix (only for square case).

    arguments
        Aph
        m = []
        nvp.method {mustBeMember(nvp.method, {'cell2mat', 'cat'})} = 'cell2mat';
    end

    [Aph, h1, h2, nvp.method] = prepToeplitzInput(Aph, m, nvp.method);
    n1 = size(Aph, 1);
    n2 = size(Aph, 2);

    switch nvp.method
        case 'cell2mat'
            c = cell(n1, n2);
            for xi = 1:n1
                for yi = 1:n2
                    val_vec = squeeze(Aph(xi, yi, :));
                    % Formula idx = i - j + 2h2 + 1
                    col = val_vec((1 : 2 * h1 + 1) + 2 * h2);
                    row = val_vec(2 - (1 : 2 * h2 + 1) + 2 * h2);
                    c{xi, yi} = toeplitz(col, row);
                end
            end
            OutM = cell2mat(c);

        case 'cat'
            OutM = [];
            for xi = 1:n1
                line = [];
                for yi = 1:n2
                    val_vec = squeeze(Aph(xi, yi, :));
                    col = val_vec((1 : 2 * h1 + 1) + 2 * h2);
                    row = val_vec(2 - (1 : 2 * h2 + 1) + 2 * h2);
                    line = [line, toeplitz(col, row)];
                end
                OutM = [OutM; line];
            end
    end

    if nargout == 2
        if h1 == h2
            N = N_tb(n1, h1, 1);
        else
            warning('array2TBlocks:nonSquareN', 'N matrix is only defined for square Toeplitz-Blocks truncations (h1=%d ~= h2=%d); returning empty N.', h1, h2)
            N = [];
        end
    end
end