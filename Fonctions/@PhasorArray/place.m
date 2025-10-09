function [K,P,res] = place(A,B,poles,varg)
    %PLACE Compute harmonic pole placement for periodic state-space systems.
    %
    %   [K, P, res] = PLACE(A, B, poles, varg) computes the harmonic pole
    %   placement for the periodic system defined by PhasorArray matrices A
    %   and B. The function determines the state-feedback matrix K to place
    %   the system poles at specified locations.
    %
    %   Inputs:
    %     A      - (PhasorArray) The system matrix in phasor form.
    %     B      - (PhasorArray) The input matrix in phasor form.
    %     poles  - (vector) The desired closed-loop poles.
    %     varg   - (Optional) Name-value pair arguments:
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
    %   See also: Sylv_harmonique, PhasorArray
    %
arguments
    A
    B
    poles
    varg.hG = A.h
    varg.hLyap = A.h*4
    varg.G = []
    varg.T = 2*pi
    varg.checkP = true
    varg.tolCheckP = 1e-8
end

assert(isvector(poles));
assert(issquare(A));
assert(numel(poles) == size(A,1));
B = PhasorArray(B);

%convert this list into a diagonal phasor array
La = PhasorArray(diag(poles));
nx = size(A,1);
if isempty(varg.G)
    if issquare(B)
        varg.G=PhasorArray.random(nx,nx,varg.hG,"time_structure","sdp");
    else
        nu = size(B,2);
        varg.G=PhasorArray.random(nu,nx,varg.hG);
    end
end


P = PhasorArray(Sylv_harmonique(-A,La,(B*varg.G),varg.hLyap,2*pi/varg.T));

if varg.checkP
    E = P.HmqEig();
    if nnz(abs(E)<varg.tolCheckP)>0
        error("P is close to non invertible, change G")
    end
end

K = varg.G/P;

if nargout>3
    res = P.d(T) + (-A)*P + P*La + B*G;
end


end

