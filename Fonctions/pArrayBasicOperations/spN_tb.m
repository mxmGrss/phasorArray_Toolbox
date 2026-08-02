function [N, Nw] = spN_tb(n, h, T)
    %SPN_TB Sparse differentiation matrix for Toeplitz-Block (TB) representation.
    %
    %   N = spN_tb(n, h, T) returns a sparse differentiation matrix N.
    %
    %   Inputs:
    %       n   - Size of the spatial dimension.
    %       h  - Harmonic order(s). Scalar h for square case, [h1, h2] for rectangular.
    %       T   - Period (default 2*pi, matching N_tb).
    %
    %   See also: N_tb, spN_bt.

    arguments
        n
        h
        T = 2*pi
    end

    % Resolve n
    if ~isscalar(n)
        n = size(n, 1);
    end

    % Resolve omega
    w = 2 * pi / T;

    % Resolve nh1, nh2
    if isscalar(h)
        h1 = h;
        h2 = h;
    else
        h1 = h(1);
        h2 = h(2);
    end

    % Construct sparse Nw: (2*h1+1) x (2*h2+1)
    % Harmonics match if k_out = k_in
    k1 = (-h1:h1)';
    k2 = (-h2:h2)';
    
    % Intersection of harmonic ranges
    k_overlap = intersect(k1, k2);
    if isempty(k_overlap)
        if isFunny(w)
            Nw = zeros(2 * h1 + 1, 2 * h2 + 1, 'like', w);
        else
            Nw = sparse(2 * h1 + 1, 2 * h2 + 1);
        end
    else
        idx1 = find(ismember(k1, k_overlap));
        idx2 = find(ismember(k2, k_overlap));
        if isFunny(w)
            % For special types, build dense and let kron handle it
            Nw = zeros(2 * h1 + 1, 2 * h2 + 1, 'like', w);
            for ii = 1:length(k_overlap)
                Nw(idx1(ii), idx2(ii)) = 1i * k_overlap(ii) * w;
            end
        else
            Nw = sparse(idx1, idx2, 1i * k_overlap * w, 2 * h1 + 1, 2 * h2 + 1);
        end
    end

    N = kron(speye(n), Nw);
end
