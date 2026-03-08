function [Tmat, N] = array2BToepliz(Aph, m, varg)
    %ARRAY2BTOEPLIZ Convert a 3D array into a Block-Toeplitz (BT) matrix.
    %
    %   Usage:
    %       [Tmat, N] = array2BToepliz(Aph)
    %       [Tmat, N] = array2BToepliz(Aph, m)
    %       [Tmat, N] = array2BToepliz(Aph, [h1, h2])
    %
    %   Inputs:
    %       Aph   - Input data. 3D double array or PhasorArray object.
    %       m     - (Optional) Order of the output block matrix.
    %               - If scalar h: Square block matrix of size (2*h+1) x (2*h+1).
    %               - If vector [h1, h2]: Rectangular matrix of size (2*h1+1) x (2*h2+1).
    %               - Default: Current order of Aph.
    %       varg  - Name-Value arguments.
    %
    %   Outputs:
    %       Tmat  - Block-Toeplitz matrix.
    %       N     - Differentiation matrix (only for square case).

    arguments
        Aph
        m = []
        varg.method {mustBeMember(varg.method, {'cell2mat', 'cat'})} = 'cell2mat';
    end

    % Detect special types
    is_special = isa(Aph, 'ndsdpvar') || isa(Aph, 'sdpvar') || isa(Aph, 'sym');

    if isa(Aph, 'PhasorArray')
        data = Aph.value;
        is_special = is_special || isa(data, 'ndsdpvar') || isa(data, 'sdpvar') || isa(data, 'sym');
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

    % BT Logic: Block (i, j) corresponds to harmonic k = j - i + (h1 - h2)
    % Required range: k in [-(h1+h2), h1+h2]
    h_req = h1 + h2;
    k_min = -h_req;
    k_max =  h_req;
    
    % Pad/truncate Aph using cat to preserve special types
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
            Aph2 = mat2cell(Aph, n1, n2, ones(2 * h_req + 1, 1));
            % Formula: idx = j - i + 2*h1 + 1
            col = 2 - (1 : 2 * h1 + 1) + 2 * h1;
            row = (1 : 2 * h2 + 1) + 2 * h1;
            Toe = toeplitz(col, row);
            Tmat = cell2mat(Aph2(Toe));

        case 'cat'
            Tmat = [];
            for xi = 1:2*h1+1
                line = [];
                for yi = 1:2*h2+1
                    k = yi - xi + h1 - h2;
                    line = [line, Aph(:, :, k + h_req + 1)];
                end
                Tmat = [Tmat; line];
            end
    end

    if nargout == 2
        if h1 == h2
            N = N_bt(n1, h1, 1);
        else
            warning('N matrix is only defined for square Block-Toeplitz truncations.');
            N = [];
        end
    end
end