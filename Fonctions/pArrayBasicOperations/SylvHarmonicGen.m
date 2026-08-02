function [Xph, M, colQ, colX] = SylvHarmonicGen(Ahmad, Bhmad, Chmad, Ea, Eb, h, omega, method, derivativeForm, direction)
% SYLVHARMONICGEN Solve generalized (descriptor) harmonic Sylvester equations.
%
%   SYNTAX:
%   X = SylvHarmonicGen(Ahmad, Bhmad, Chmad, Ea, Eb, h, omega) solves
%   (derivativeForm = 'product', default):
%       d/dt( Ea*X*Eb ) + Ahmad*X*Eb + Ea*X*Bhmad + Chmad = 0
%   or, with derivativeForm = 'sandwich':
%       Ea*(dX/dt)*Eb  + Ahmad*X*Eb + Ea*X*Bhmad + Chmad = 0
%
%   With direction = 'forward' the sign of the derivative term flips
%   (covariance-type integration, cf. SylvHarmonic), e.g. for 'product':
%       d/dt( Ea*X*Eb ) = Ahmad*X*Eb + Ea*X*Bhmad + Chmad
%   The direction is orthogonal to derivativeForm: it selects the sign of
%   the derivative operator, derivativeForm selects what it acts on.
%
%   The two forms differ by the dEa/dt, dEb/dt terms and coincide when the
%   mass matrices are constant. 'product' arises from descriptor systems
%   with the mass OUTSIDE the derivative (E*x' = A*x, e.g. Lyapunov of the
%   primal); 'sandwich' arises from the ADJOINT / mass-inside form
%   (d/dt(E*x) = A*x, e.g. covariance equations of the dual filter).
%   In the harmonic domain the difference is only the operator ordering:
%   N*T(M_E) ('product') versus T(M_E)*N ('sandwich').
%   See the "CHOOSING derivativeForm" section in PhasorArray/lyapG for the
%   full discussion (associated Lyapunov functions, duality, use cases).
%
%   With Ea = Eb = I both reduce to the standard Sylvester equation solved
%   by SylvHarmonic:  dX/dt + Ahmad*X + X*Bhmad + Chmad = 0.
%
%   MATHEMATICAL STRATEGY (validated in scratch/wip_lyapunov_generalized):
%   Time-domain vectorisation via vec(U*V*W) = (W.' kron U) vec(V):
%       d/dt[ M_E(t) p ] + M_S(t) p + vec(C) = 0
%       M_E = Eb.' kron Ea
%       M_S = Eb.' kron A  +  B.' kron Ea
%   In the harmonic (Toeplitz) domain the product-rule derivative is handled
%   implicitly — neither inv(E) nor dE/dt is ever formed:
%       [ N * T(M_E) + T(M_S) ] * F(vec(X)) = -F(vec(C))
%   N * T(M_E) is NOT diagonal (unlike the standard case), hence the dedicated
%   dense-operator construction below instead of the blkdiag path of SylvHarmonic.
%
%   INPUTS:
%       Ahmad, Bhmad, Chmad : PhasorArray (or 3D harmonic arrays / numeric).
%                             A: nA x nA, B: nB x nB, C: nA x nB.
%       Ea, Eb              : Mass matrices, PhasorArray or numeric (possibly
%                             time-varying). [] defaults to identity.
%                             Ea: nA x nA, Eb: nB x nB.
%       h                   : Harmonic truncation order for unknown X (hIn).
%       omega               : Fundamental pulsation [rad/s].
%       method              : 'rectangle' (overdetermined, default) or 'square'
%                             (Galerkin truncation, hOut = hIn).
%       derivativeForm      : 'product' (default) or 'sandwich' — what the
%                             derivative acts on (see above).
%       direction           : 'backward' (default) or 'forward' — sign of the
%                             derivative term.
%
%   OUTPUTS:
%       Xph  : 3D harmonic array of the solution (same raw format as SylvHarmonic).
%       M    : Harmonic operator, M = -T(M_S) - N*T(M_E).
%       colQ : Vectorised RHS, F(vec(C)) at order hOut.
%       colX : Solution vector, colX = M \ colQ (least-squares if rectangle).
%
%   See also: SylvHarmonic, PhasorArray/lyapG, PhasorArray/lyap, N_tb.

    arguments
        Ahmad
        Bhmad
        Chmad
        Ea                                                    % [] → identity
        Eb                                                    % [] → identity
        h     (1,1) double {mustBeInteger, mustBeNonnegative}
        omega (1,1)
        method {mustBeMember(method, {'rectangle','square'})} = 'rectangle'
        derivativeForm {mustBeMember(derivativeForm, {'product','sandwich'})} = 'product'
        direction {mustBeMember(direction, {'backward','forward'})} = 'backward'
    end

    % Normalise everything to PhasorArray (kron/TB algebra requires it).
    if ~isa(Ahmad, 'PhasorArray'), Ahmad = PhasorArray(Ahmad); end
    if ~isa(Bhmad, 'PhasorArray'), Bhmad = PhasorArray(Bhmad); end
    if ~isa(Chmad, 'PhasorArray'), Chmad = PhasorArray(Chmad); end

    nA = size(Ahmad, 1);
    nB = size(Bhmad, 1);

    % Empty mass matrices default to (constant) identity → standard Sylvester.
    if isempty(Ea), Ea = PhasorArray(eye(nA)); end
    if isempty(Eb), Eb = PhasorArray(eye(nB)); end
    if ~isa(Ea, 'PhasorArray'), Ea = PhasorArray(Ea); end
    if ~isa(Eb, 'PhasorArray'), Eb = PhasorArray(Eb); end

    % Kronecker operators in the phasor domain (Cauchy products inside kron).
    % NOTE: .' is the time-domain transpose required by the vec identity.
    ME = kron(Eb.', Ea);                            % (nA*nB) x (nA*nB)
    MS = kron(Eb.', Ahmad) + kron(Bhmad.', Ea);     % (nA*nB) x (nA*nB)

    % Truncation: unknown on hIn, equations enforced on hOut.
    % 'rectangle' inflates hOut by the operator spectral width so that all
    % product harmonics (incl. those of d/dt[M_E p]) are captured in L2.
    hIn = h;
    hC  = Chmad.h;
    if strcmp(method, 'square')
        hOut = hIn;
    else
        hOut = max([hIn + MS.h, hIn + ME.h, hC]);
    end

    % Rectangular Toeplitz-Block operators acting on F(vec(X)) at order hIn.
    T_MS = sparse(TB(MS, [hOut, hIn]));
    T_ME = sparse(TB(ME, [hOut, hIn]));

    % Harmonic derivative operator: N = I_{nA*nB} kron diag(j*k*omega).
    % 'product' : d/dt[M_E p]  →  N (at hOut) LEFT of T(M_E)
    % 'sandwich': M_E dp/dt    →  N (at hIn)  RIGHT of T(M_E)
    if strcmp(derivativeForm, 'sandwich')
        Nin    = sparse(N_tb(nA * nB, hIn, [], "omega", omega));
        Mderiv = T_ME * Nin;
    else
        Nout   = sparse(N_tb(nA * nB, hOut, [], "omega", omega));
        Mderiv = Nout * T_ME;
    end

    % 'backward':  derivative term + M_S p + vec(C) = 0  →  M = -T_MS - Mderiv
    % 'forward' : -derivative term + M_S p + vec(C) = 0  →  M = -T_MS + Mderiv
    if strcmp(direction, 'forward')
        M = -T_MS + Mderiv;
    else
        M = -T_MS - Mderiv;
    end

    % Vectorised RHS: F(vec(C)) padded/truncated to hOut.
    colQ = FvTB(Chmad, hOut);

    % Least-squares solve of the (possibly rectangular) harmonic system.
    colX = M \ colQ;

    % Reconstruct the solution and return the raw 3D harmonic array
    % (same output contract as SylvHarmonic).
    Xpa = PhasorArray.fromF_tb(full(colX), nA, nB);
    Xph = pvalue(Xpa);
end
