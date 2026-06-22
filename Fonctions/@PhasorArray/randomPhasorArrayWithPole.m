function obj = randomPhasorArrayWithPole(nx,poles,T,varg)
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
        varg.h = []
        varg.BG = []
        varg.isReal = true
        varg.badPoleTol = 25/100; % Tolerance for the number of bad poles
        varg.poleTol = 1e-2; % Tolerance for the pole values compared to the desired poles
    end
    if isempty(varg.h) && isempty(varg.BG)
        varg.h=5;
        A  = PhasorArray.random(nx,nx,varg.h);
        BG = PhasorArray.random(nx,nx,varg.h);
    elseif isempty(varg.BG)
        A  = PhasorArray.random(nx,nx,varg.h);
        BG = PhasorArray.random(nx,nx,varg.h);
    elseif isempty(varg.h)
        varg.h = varg.BG.h;
        A  = PhasorArray.random(nx,nx,varg.h);
    else
        A  = PhasorArray.random(nx,nx,varg.h);
    end
    assert(numel(poles)==nx)
    La = PhasorArray(diag(poles));

    %Solve the appropriate Sylvester equation
    P = PhasorArray(SylvHarmonic(-A,La,BG,4*varg.h,2*pi/T));
    %Compute K
    K = BG/P;
    %compute the new A with appropriate eigen values
    obj = A-K;

    % Check if the generated PhasorArray has the desired poles
    E = HmqNEig(trunc(obj,varg.h),4*varg.h,T);

    % Check if the eigenvalues match the desired poles
    matchedPoles = sum(any(abs(real(E) - poles) < varg.poleTol, 2));

    prop = matchedPoles / numel(E);
    if prop < 1 - varg.badPoleTol
        warning('PhasorArray:randomWithPole:desiredPolesMissed', 'Generated PhasorArray does not match desired poles (%.1f%% match).', prop*100);

        %recall the function with the same parameters
        %convert varg to cell
        varg = namedargs2cell(varg);
        obj = PhasorArray.randomPhasorArrayWithPole(nx,poles,T,varg{:});
    end
end

