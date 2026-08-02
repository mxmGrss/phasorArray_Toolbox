function [OutM, N] = sparray2TBlocks(Aph, m, nvp)
    %SPARRAY2TBLOCKS Convert a 3D array into a sparse Toeplitz-blocks (TB) matrix.
    %
    %   Usage:
    %       [OutM, N] = sparray2TBlocks(Aph, [h1, h2])
    %
    %   Inputs:
    %       Aph   - Input data (3D array or PhasorArray).
    %       m     - Output order(s) [h1, h2].
    %
    %   See also: array2TBlocks.

    arguments
        Aph
        m = []
        nvp.method = 'cell2mat';
    end

    [Aph, h1, h2, nvp.method] = prepToeplitzInput(Aph, m, nvp.method);
    n1 = size(Aph, 1);
    n2 = size(Aph, 2);

    % Keyed on the DATA, not on the method: a caller may legitimately ask for
    % 'cat' on numeric input. Only sym/sdpvar are unsupported here, and saying
    % so beats letting cell2mat fail with an opaque message.
    if isFunny(Aph)
        error('PhasorArray:sparray2TBlocks:specialType', ...
            ['sparray2TBlocks does not support sym/sdpvar coefficients ', ...
             '(sparse storage cannot hold them). Use array2TBlocks instead.'])
    end

    % Both methods build the same blocks (idx = i - j + 2*h2 + 1) and differ
    % only in whether the coefficient vector is made sparse first.
    useSparse = strcmp(nvp.method, 'cell2mat');
    c = cell(n1, n2);
    for xi = 1:n1
        for yi = 1:n2
            val_vec = squeeze(Aph(xi, yi, :));
            if useSparse, val_vec = sparse(val_vec); end
            col = val_vec((1 : 2 * h1 + 1) + 2 * h2);
            row = val_vec(2 - (1 : 2 * h2 + 1) + 2 * h2);
            c{xi, yi} = toeplitz(col, row);
        end
    end
    OutM = cell2mat(c);

    if nargout == 2
        if h1 == h2
            N = N_tb(n1, h1, 1);
        else
            N = [];
        end
    end
end