function [F, G, cert] = gramConstraint(M, margin, basis)
%GRAMCONSTRAINT  Exact LMI certificate of M(t) >= margin*I for all t.
%
%   [F, G, CERT] = GRAMCONSTRAINT(M, MARGIN) returns YALMIP constraints F that
%   hold if and only if the Hermitian trigonometric polynomial M satisfies
%
%       M(t) - margin*I = V(t)' * G * V(t) ,   G >= 0 ,   for all t,
%
%   with V(t) = [I; e^{jt} I; ...; e^{jdt} I] and d = M.h (matrix Fejer-Riesz
%   theorem, Gram form). Unlike a finite Toeplitz section T_tb(M,r) >= 0,
%   which is only a necessary condition for any r, this is exact and involves
%   no section order: its size is set by the degree of M. Reference: Vernerey,
%   Riedinger & Daafouz, "A TBLMI framework for harmonic robust control",
%   corrected version, Theorem 9 and Proposition 3.
%
%   SYNTAX
%     [F, G, cert] = M.gramConstraint(margin)
%     [F, G, cert] = M.gramConstraint(margin, basis)
%
%   INPUTS
%     M       Square n x n PhasorArray, Hermitian in time (M(t) = M(t)').
%             Coefficients may be numeric or YALMIP expressions, e.g. an LMI
%             built from PhasorArray.ndsdpvar unknowns.
%     margin  Scalar double or sdpvar: the bound in M(t) >= margin*I.
%     basis   "hermitian" (default): complex Hermitian Gram matrix of order
%             n(d+1), handed to the solver as a real cone of order 2n(d+1).
%             "real": real symmetric Gram matrix of order n(d+1) in the real
%             trigonometric basis (half-integer frequencies when d is odd),
%             Proposition 3. Requires M real-valued and symmetric in time;
%             halves the order of the cone for the same set of polynomials.
%
%   OUTPUTS
%     F     YALMIP constraints: G >= 0 and the d+1 (hermitian) or 2d+1 (real)
%           coefficient equalities.
%     G     Gram decision matrix (sdpvar) of order n(d+1).
%     cert  struct(basis, degree, n, basisChange) for CHECKGRAMCERTIFICATE.
%
%   EXAMPLES
%     P = PhasorArray.ndsdpvar(n, n, q);
%     L = -(P.d(T) + A.'*P + P*A + Q);
%     [F1, G1, c1] = L.gramConstraint(1e-6, "real");
%     [F2, G2, c2] = P.gramConstraint(1e-6, "real");
%     optimize([F1, F2], trace(phas(P, 0)));
%     r = sdpval(L).checkGramCertificate(1e-6, value(G1), c1);
%
%   Requires YALMIP.
%
%   See also: checkGramCertificate, ndsdpvar, T_tb.

arguments
    M (:,:) PhasorArray
    margin
    basis (1,1) string {mustBeMember(basis, ["hermitian", "real"])} = "hermitian"
end

if ~exist('sdpvar', 'file')
    error('PhasorArray:gramConstraint:noYalmip', 'gramConstraint requires YALMIP.');
end
n = size(M, 1);
if size(M, 2) ~= n
    error('PhasorArray:gramConstraint:notSquare', 'M must be square.');
end
if ~isscalar(margin)
    error('PhasorArray:gramConstraint:margin', 'margin must be a scalar.');
end
d = h(M);
coefficients = value(M);
if isnumeric(coefficients)
    % Numeric data can be checked here; YALMIP data is the caller's contract.
    scale = max(1, norm(coefficients(:)));
    hermitianDefect = coefficients - permute(conj(flip(coefficients, 3)), [2 1 3]);
    if norm(hermitianDefect(:)) > 1e-11*scale
        error('PhasorArray:gramConstraint:notHermitian', 'M(t) must be Hermitian.');
    end
    if basis == "real" && norm(reshape(coefficients - conj(flip(coefficients, 3)), [], 1)) > 1e-11*scale
        error('PhasorArray:gramConstraint:notReal', 'basis "real" requires a real-valued M.');
    end
end

gramOrder = n*(d + 1);
switch basis
    case "hermitian"
        % Coefficient k >= 0 of V'GV is the sum of the blocks of G on block
        % diagonal k: M_k - margin*delta_k0*I = sum_r G(r, r+k).
        % YALMIP's sdpvar parser needs char type flags, not strings.
        G = sdpvar(gramOrder, gramOrder, 'hermitian', 'complex');
        F = G >= 0;
        for k = 0:d
            gramCoefficient = 0;
            for rowBlock = 0:d-k
                rows = rowBlock*n + (1:n);
                cols = (rowBlock + k)*n + (1:n);
                gramCoefficient = gramCoefficient + G(rows, cols);
            end
            target = phas(M, k);
            if k == 0, target = target - margin*eye(n); end
            F = [F, gramCoefficient == target]; %#ok<AGROW>
        end
        basisChange = speye(gramOrder);
    case "real"
        % Real Gram matrix G in the real trigonometric basis U(t); the
        % Hermitian Gram in the Fourier basis is H = B*G*B', B unitary, and
        % its coefficient map is applied in independent real coordinates
        % (real/imaginary Fourier pairs x symmetric matrix entries).
        G = sdpvar(gramOrder, gramOrder, 'symmetric');
        basisChange = kron(conj(centeredRealBasis(d)), speye(n));
        coefficientMap = gramMap(n, d);
        independentBasis = kron(fourierBasis(d), symmetricBasis(n));
        realGramMap = real(independentBasis'*coefficientMap*kron(conj(basisChange), basisChange));
        targets = cell(2*d + 1, 1);
        for k = -d:d
            target = phas(M, k);
            if k == 0, target = target - margin*eye(n); end
            targets{k + d + 1} = target(:);
        end
        targetCoefficients = real(independentBasis'*vertcat(targets{:}));
        F = [G >= 0, realGramMap*G(:) == targetCoefficients];
end
cert = struct('basis', basis, 'degree', d, 'n', n, 'basisChange', basisChange);
end

function coefficientMap = gramMap(n, d)
% Sparse map vec(H) -> stacked vec(coefficient k), k = -d..d, of V'HV:
% block (r, c) of H contributes to coefficient c - r.
gramOrder = n*(d + 1);
rows = zeros(gramOrder^2, 1); cols = zeros(gramOrder^2, 1);
offset = 0;
[rowIndex, colIndex] = ndgrid(1:n, 1:n);
for rowBlock = 0:d
    for colBlock = 0:d
        ids = offset + (1:n^2); offset = offset + n^2;
        rows(ids) = reshape(rowIndex + (colIndex - 1)*n + (colBlock - rowBlock + d)*n^2, [], 1);
        cols(ids) = reshape(rowBlock*n + rowIndex + (colBlock*n + colIndex - 1)*gramOrder, [], 1);
    end
end
coefficientMap = sparse(rows, cols, 1, n^2*(2*d + 1), gramOrder^2);
end

function basis = fourierBasis(d)
% Unitary map from independent real coefficients (constant, then cos/sin
% pairs) to the conjugate Fourier coefficients -d..d.
basis = sparse(2*d + 1, 2*d + 1);
basis(d + 1, 1) = 1;
for k = 1:d
    basis(d + 1 + k, 2*k) = 1/sqrt(2); %#ok<SPRIX>
    basis(d + 1 - k, 2*k) = 1/sqrt(2); %#ok<SPRIX>
    basis(d + 1 + k, 2*k + 1) = -1i/sqrt(2); %#ok<SPRIX>
    basis(d + 1 - k, 2*k + 1) = 1i/sqrt(2); %#ok<SPRIX>
end
end

function basis = symmetricBasis(n)
% Orthonormal map from the n(n+1)/2 independent entries to vec of a
% symmetric matrix.
basis = sparse(n^2, n*(n + 1)/2);
column = 0;
for col = 1:n
    for row = 1:col
        column = column + 1;
        if row == col
            basis(row + (col - 1)*n, column) = 1; %#ok<SPRIX>
        else
            basis(row + (col - 1)*n, column) = 1/sqrt(2); %#ok<SPRIX>
            basis(col + (row - 1)*n, column) = 1/sqrt(2); %#ok<SPRIX>
        end
    end
end
end

function basis = centeredRealBasis(d)
% Minimal real trigonometric basis of the d+1 Gram slots, centred on d/2:
% integer frequencies 0..d/2 if d is even, half-integers 1/2..d/2 if d is odd.
if mod(d, 2) == 0
    basis = fourierBasis(d/2);
    return
end
pairCount = (d + 1)/2;
rows = zeros(4*pairCount, 1); cols = zeros(4*pairCount, 1);
values = complex(rows);
for pairIndex = 1:pairCount
    ids = 4*pairIndex + (-3:0);
    negativeRow = pairCount - pairIndex + 1;
    positiveRow = pairCount + pairIndex;
    rows(ids) = [negativeRow, positiveRow, negativeRow, positiveRow];
    cols(ids) = [2*pairIndex - 1, 2*pairIndex - 1, 2*pairIndex, 2*pairIndex];
    values(ids) = [1, 1, 1i, -1i]/sqrt(2);
end
basis = sparse(rows, cols, values, d + 1, d + 1);
end
