function [N, Nw] = N_bt(n, h, T)
    %N_BT Differentiation matrix for Block-Toeplitz (BT) representation.
    %
    %   N = N_bt(n, h, T) returns the differentiation matrix N such that 
    %   if X is a BT vector, N*X corresponds to dX/dt.
    %
    %   N = diag(jkω) ⊗ eye(n)
    %
    %   Inputs:
    %       n   - Size of the spatial dimension.
    %       h  - Harmonic order(s). Scalar h for square case, [h1, h2] for rectangular.
    %       T   - Period (default 2*pi, matching N_tb).
    %
    %   See also: N_tb, spN_bt, array2BToeplitz.

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

    % Construct Nw: (2*h1+1) x (2*h2+1)
    k_ranges = -h1:h1;
    k_inputs = -h2:h2;
    [K_in, K_out] = meshgrid(k_inputs, k_ranges);
    Nw = zeros(2 * h1 + 1, 2 * h2 + 1);
    mask = (K_in == K_out);
    Nw(mask) = 1i * K_in(mask) * w;

    N = kron(Nw, eye(n));
end
