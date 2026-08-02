function [N, Nw] = N_tb(n, h, T, nvp)
    %N_TB Differentiation matrix for Toeplitz-Block (TB) representation.
    %
    %   N = N_tb(n, h, T) returns the differentiation matrix N such that 
    %   if X is a TB vector, N*X corresponds to dX/dt.
    %
    %   N = eye(n) ⊗ diag(jkω)
    %
    %   Inputs:
    %       n   - Size of the spatial dimension.
    %       h  - Harmonic order(s). Scalar h for square case, [h1, h2] for rectangular.
    %       T   - Period (default 2*pi).
    %       nvp - Name-Value:
    %           * omega: Provide omega directly instead of T.
    %
    %   See also: N_bt, spN_tb, array2TBlocks.

    arguments
        n
        h
        T = []
        nvp.omega = []
    end

    % Resolve omega
    if isempty(nvp.omega)
        if isempty(T)
            T = 2 * pi;
        end
        if isinf(T)
            w = 0;
        else
            w = 2 * pi / T;
        end
    else
        w = nvp.omega;
    end

    if isnan(w)
        error('PhasorArray:N_tb:NaN', 'omega cannot be NaN');
    end

    % Resolve n
    if ~isscalar(n)
        n = size(n, 1);
    end

    % Resolve nh1, nh2
    if isscalar(h)
        h1 = h;
        h2 = h;
    else
        h1 = h(1);
        h2 = h(2);
    end

    % Construct Nw: (2*h1+1) x (2*h2+1)
    % Entries jkω where harmonics match
    k_ranges = -h1:h1;
    k_inputs = -h2:h2;
    
    [K_in, K_out] = meshgrid(k_inputs, k_ranges);
    Nw = zeros(2 * h1 + 1, 2 * h2 + 1);
    mask = (K_in == K_out);
    Nw(mask) = 1i * K_in(mask) * w;

    N = kron(eye(n), Nw);
end
