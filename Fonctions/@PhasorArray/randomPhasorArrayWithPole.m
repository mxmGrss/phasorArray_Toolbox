function obj = randomPhasorArrayWithPole(nx,poles,T,nvp)
    % RANDOMPHASORARRAYWITHPOLE Generate a PhasorArray with prescribed poles.
    %
    %   OBJ = RANDOMPHASORARRAYWITHPOLE(NX, POLES, T, VARG) constructs a
    %   time-periodic matrix `A(t)` with the specified eigenvalues (`poles`).
    %   A first random matrix `A` is generated, and a B*G term is computed to
    %   ensure the desired eigenvalues by performing pole placement.
    %   The resulting PhasorArray `OBJ` is the difference `A-BK`, where `K` is
    %   the obtained feedback gain.
    %
    %   Inputs:
    %     NX    - (integer) The size of the square matrix.
    %     POLES - (vector) Desired eigenvalues for `A(t)`, must have length NX.
    %     T     - (scalar, optional) Period of the matrix. Default is `2*pi`.
    %     VARG  - (optional) Name-value pair arguments:
    %               - 'h'      (integer) Harmonic order (default: `5`).
    %               - 'BG'     (PhasorArray) Custom B*G term (default: random).
    %               - 'isReal' (logical) Force a real-valued matrix (default: `true`).
    %
    %   Outputs:
    %     OBJ - (PhasorArray) Generated PhasorArray with specified poles.
    %
    %   See also: SylvHarmonic, random
    arguments
        nx
        poles
        T  = 2*pi
        nvp.h = []
        nvp.BG = []
        nvp.isReal = true
        nvp.badPoleTol = 25/100; % Tolerance for the number of bad poles
        nvp.poleTol = 1e-2; % Tolerance for the pole values compared to the desired poles
    end
    if isempty(nvp.h) && isempty(nvp.BG)
        nvp.h=5;
        A  = PhasorArray.random(nx,nx,nvp.h);
        BG = PhasorArray.random(nx,nx,nvp.h);
    elseif isempty(nvp.BG)
        A  = PhasorArray.random(nx,nx,nvp.h);
        BG = PhasorArray.random(nx,nx,nvp.h);
    elseif isempty(nvp.h)
        nvp.h = nvp.BG.h;
        A  = PhasorArray.random(nx,nx,nvp.h);
    else
        A  = PhasorArray.random(nx,nx,nvp.h);
    end
    assert(numel(poles)==nx)
    La = PhasorArray(diag(poles));

    %Solve the appropriate Sylvester equation
    P = PhasorArray(SylvHarmonic(-A,La,BG,4*nvp.h,2*pi/T));
    %Compute K
    K = BG/P;
    %compute the new A with appropriate eigen values
    obj = A-K;

    % Check if the generated PhasorArray has the desired poles
    E = HmqNEig(trunc(obj,nvp.h),4*nvp.h,T);

    % Check if the eigenvalues match the desired poles
    matchedPoles = sum(any(abs(real(E) - poles) < nvp.poleTol, 2));

    prop = matchedPoles / numel(E);
    if prop < 1 - nvp.badPoleTol
        warning('PhasorArray:randomWithPole:desiredPolesMissed', 'Generated PhasorArray does not match desired poles (%.1f%% match).', prop*100);

        %recall the function with the same parameters
        %convert nvp to cell
        nvp = namedargs2cell(nvp);
        obj = PhasorArray.randomPhasorArrayWithPole(nx,poles,T,nvp{:});
    end
end

