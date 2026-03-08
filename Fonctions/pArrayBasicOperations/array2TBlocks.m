function [OutM, N] = array2TBlocks(Aph, m, varg)
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
    %       varg  - Name-Value arguments:
    %               * 'method': 'cell2mat' (fast) or 'cat' (SDP/SYM compatible).
    %
    %   Outputs:
    %       OutM  - Toeplitz-Blocks matrix.
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

    % Required range: k in [-(h1+h2), h1+h2]
    h_req = h1 + h2;
    k_min = -h_req;
    k_max =  h_req;
    
    % Map overlap
    overlap_k_min = max(k_min, -nh_in);
    overlap_k_max = min(k_max,  nh_in);
    
    if overlap_k_min <= overlap_k_max
        % Pad/truncate to exactly match [-h_req, h_req]
        % Use cat to preserve potential special types (ndsdpvar/sdpvar/sym)
        data_sliced = Aph(:, :, (nh_in + 1 + overlap_k_min) : (nh_in + 1 + overlap_k_max));
        n_lead = overlap_k_min - k_min;
        n_trail = k_max - overlap_k_max;
        
        z1 = zeros(n1, n2, n_lead);
        z2 = zeros(n1, n2, n_trail);
        Aph = cat(3, z1, data_sliced, z2);
    else
        Aph = zeros(n1, n2, 2 * h_req + 1);
    end
    
    % Index in Aph is k + h_req + 1
    
    switch varg.method
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
            warning('N matrix is only defined for square Toeplitz-Blocks truncations.');
            N = [];
        end
    end
end