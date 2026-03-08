function [N, Nw] = N_tb(n, nh, T, varg)
    %N_TB Differentiation matrix for Toeplitz-Block (TB) representation.
    %
    %   N = N_tb(n, nh, T) returns the differentiation matrix N such that 
    %   if X is a TB vector, N*X corresponds to dX/dt.
    %
    %   N = eye(n) ⊗ diag(jkω)
    %
    %   Inputs:
    %       n   - Size of the spatial dimension.
    %       nh  - Harmonic order(s). Scalar h for square case, [h1, h2] for rectangular.
    %       T   - Period (default 2*pi).
    %       varg - Name-Value:
    %           * omega: Provide omega directly instead of T.
    %
    %   See also: N_bt, spN_tb, array2TBlocks.

    arguments
        n
        nh
        T = []
        varg.omega = []
    end

    % Resolve omega
    if isempty(varg.omega)
        if isempty(T)
            T = 2 * pi;
        end
        if isinf(T)
            w = 0;
        else
            w = 2 * pi / T;
        end
    else
        w = varg.omega;
    end

    if isnan(w)
        error('PhasorArray:N_tb:NaN', 'omega cannot be NaN');
    end

    % Resolve n
    if ~isscalar(n)
        n = size(n, 1);
    end

    % Resolve nh1, nh2
    if isscalar(nh)
        h1 = nh;
        h2 = nh;
    else
        h1 = nh(1);
        h2 = nh(2);
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
