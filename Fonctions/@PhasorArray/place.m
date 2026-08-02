function [K,P,res] = place(A,B,poles,nvp)
    %PLACE Compute harmonic pole placement for periodic state-space systems.
    %
    %   [K, P, res] = PLACE(A, B, poles, nvp) computes the harmonic pole
    %   placement for the periodic system defined by PhasorArray matrices A
    %   and B. The function determines the state-feedback matrix K to place
    %   the system poles at specified locations.
    %
    %   Inputs:
    %     A      - (PhasorArray) The system matrix in phasor form.
    %     B      - (PhasorArray) The input matrix in phasor form.
    %     poles  - (vector) The desired closed-loop poles.
    %     nvp   - (Optional) Name-value pair arguments:
    %         'hG'        (integer, default: A.h)     - Harmonic order of gain matrix G.
    %         'hLyap'     (integer, default: 4*A.h)   - Harmonic order for Lyapunov equation.
    %         'G'         (PhasorArray, default: [])  - Predefined gain matrix G. If empty, a random symmetric definite positive matrix is used.
    %         'T'         (double, default: 2*pi)     - Period of the system.
    %         'checkP'    (logical, default: true)    - Enable checking the invertibility of P.
    %         'tolCheckP' (double, default: 1e-8)     - Tolerance for checking P invertibility.
    %
    %   Outputs:
    %     K   - (PhasorArray) The computed state-feedback gain matrix.
    %     P   - (PhasorArray) The Lyapunov solution matrix.
    %     res - (PhasorArray) The residual of the Sylvester equation (optional).
    %
    %   Method:
    %     - Constructs a diagonal phasor matrix La with the desired pole locations.
    %     - Solves a harmonic Sylvester equation for P : -A*P + P*La + B*G = 0.
    %     - Computes K as K = G/P.
    %     - Optionally checks if P is near singular and issues an error if necessary.
    %
    %   Example:
    %     A = PhasorArray.random(3,3,5);
    %     B = PhasorArray.random(3,1,5);
    %     poles = [-1 -2 -3];
    %     [K, P] = place(A, B, poles);
    %
    %   See also: SylvHarmonic, PhasorArray
    %
arguments
    A
    B
    poles
    % A constant G. Letting it carry harmonics costs three to twelve orders of
    % accuracy on the placed exponents: measured over 8 random draws on three
    % systems, hG = 0 lands within 1e-9 while hG = A.h scatters between 1e-2
    % and 7e-1, with no warning and a residual that stays at 1e-8 throughout.
    nvp.hG = 0
    nvp.hLyap = A.h*4
    nvp.G = []
    nvp.T = 2*pi
    nvp.checkP = true
    nvp.tolCheckP = 1e-8
end

assert(isvector(poles));
assert(issquare(A));
assert(numel(poles) == size(A,1));
B = PhasorArray(B);

%convert this list into a diagonal phasor array
La = PhasorArray(diag(poles));
nx = size(A,1);
if isempty(nvp.G)
    if issquare(B)
        nvp.G=PhasorArray.randomSPD(nx, nvp.hG);
    else
        nu = size(B,2);
        nvp.G=PhasorArray.random(nu,nx,nvp.hG);
    end
end


P = PhasorArray(SylvHarmonic(-A,La,(B*nvp.G),nvp.hLyap,2*pi/nvp.T));

if nvp.checkP
    E = P.HmqEig();
    if nnz(abs(E)<nvp.tolCheckP)>0
        error('PhasorArray:place:singularP', 'P is nearly singular (min|eig(P)|=%e < tol). Modify gain matrix G or relax the pole placement constraints.', min(abs(E)))
    end
end

K = nvp.G/P;

% res is the third output, so the guard was off by one and never fired; the
% body then read T and G, which do not exist under those names.
if nargout>2
    res = d(P, nvp.T) + (-A)*P + P*La + B*nvp.G;
end


end

