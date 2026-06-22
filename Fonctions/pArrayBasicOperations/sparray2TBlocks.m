function [OutM, N] = sparray2TBlocks(Aph, m, varg)
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
        varg.method = 'cell2mat';
    end

    % Detect special types
    is_special = isFunny(Aph);

    if isa(Aph, 'PhasorArray')
        data = Aph.value;
        is_special = is_special || isFunny(data);
        Aph = data;
    end

    if is_special
        varg.method = 'cat';
    end

    % Get input dimensions
    [n1, n2, nh_len] = size(Aph);
    nh_in = (nh_len - 1) / 2;

    % Resolve output orders [h1, h2]
    if isempty(m)
        h1 = nh_in;
        h2 = nh_in;
    elseif isscalar(m)
        h1 = m;
        h2 = m;
    else
        h1 = m(1);
        h2 = m(2);
    end

    % TB Logic: k = i - j - (h1 - h2)
    % Required range: k in [-(h1+h2), h1+h2]
    h_req = h1 + h2;
    k_min = -h_req;
    k_max =  h_req;
    
    % Pad/truncate using cat to preserve special types
    overlap_k_min = max(k_min, -nh_in);
    overlap_k_max = min(k_max,  nh_in);
    
    if overlap_k_min <= overlap_k_max
        data_sliced = Aph(:, :, (nh_in + 1 + overlap_k_min) : (nh_in + 1 + overlap_k_max));
        n_lead = overlap_k_min - k_min;
        n_trail = k_max - overlap_k_max;
        
        z1 = zeros(n1, n2, n_lead);
        z2 = zeros(n1, n2, n_trail);
        Aph = cat(3, z1, data_sliced, z2);
    else
        Aph = zeros(n1, n2, 2 * h_req + 1);
    end
    
    % harmonic k is at index (k + h_req + 1)

    switch varg.method
        case 'cell2mat'
            c = cell(n1, n2);
            for xi = 1:n1
                for yi = 1:n2
                    val_vec = sparse(squeeze(Aph(xi, yi, :)));
                    % idx = i - j + 2*h2 + 1
                    col = val_vec((1 : 2 * h1 + 1) + 2 * h2);
                    row = val_vec(2 - (1 : 2 * h2 + 1) + 2 * h2);
                    c{xi, yi} = toeplitz(col, row);
                end
            end
            OutM = cell2mat(c);
        otherwise
            % Safe Assembly for special types (sdpvar, sym, etc.)
            p1 = 2 * h1 + 1;
            p2 = 2 * h2 + 1;
            c = cell(n1, n2);
            for xi = 1:n1
                for yi = 1:n2
                    val_vec = squeeze(Aph(xi, yi, :));
                    col = val_vec((1 : 2 * h1 + 1) + 2 * h2);
                    row = val_vec(2 - (1 : 2 * h2 + 1) + 2 * h2);
                    c{xi, yi} = toeplitz(col, row);
                end
            end
            OutM = cell2mat(c);
    end

    if nargout == 2
        if h1 == h2
            N = N_tb(n1, h1, 1);
        else
            N = [];
        end
    end
end