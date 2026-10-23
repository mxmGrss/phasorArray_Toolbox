function report = checkGramCertificate(M, margin, G, cert)
%CHECKGRAMCERTIFICATE  Validate a numerical Gram certificate a posteriori.
%
%   REPORT = CHECKGRAMCERTIFICATE(M, MARGIN, G, CERT) bounds min_t of the
%   smallest eigenvalue of M(t) from a numerical solution of the program of
%   GRAMCONSTRAINT. A solver returns G and M only up to its tolerances: G may
%   be slightly indefinite and the coefficient equalities slightly violated.
%   With E = M - margin*I - V'GV the coefficient defect,
%
%       M(t) >= ( margin + (d+1)*lambda_min(G) - sum_k ||E_k|| ) * I ,
%
%   because V(t)'V(t) = (d+1) I and ||E(t)|| <= sum_k ||E_k||. Reference:
%   Vernerey, Riedinger & Daafouz, "A TBLMI framework for harmonic robust
%   control", corrected version, Proposition 4.
%
%   SYNTAX
%     report = Mnum.checkGramCertificate(margin, Gnum, cert)
%
%   INPUTS
%     M       Numeric Hermitian PhasorArray, e.g. sdpval(M) after optimize.
%             Any harmonic of M beyond the Gram degree counts in E.
%     margin  Numeric scalar used in GRAMCONSTRAINT.
%     G       Numeric Gram matrix, value(G).
%     cert    Third output of GRAMCONSTRAINT.
%
%   OUTPUTS
%     report  struct with fields
%       lower                lower bound of min_t lambda_min(M(t))
%       accepted             lower > 0
%       gramMinimum          lambda_min(G)
%       coefficientResidual  sum_k ||E_k||
%       margin, degree, basis
%
%   Floating-point evaluation, no interval-arithmetic enclosure.
%
%   EXAMPLES
%     [F, G, c] = L.gramConstraint(1e-6);
%     optimize(F);
%     report = sdpval(L).checkGramCertificate(1e-6, value(G), c);
%
%   See also: gramConstraint.

arguments
    M (:,:) PhasorArray
    margin (1,1) double {mustBeReal, mustBeFinite}
    G (:,:) double {mustBeFinite}
    cert (1,1) struct
end

n = cert.n; d = cert.degree; B = cert.basisChange;
if size(M, 1) ~= n || size(M, 2) ~= n || ~isequal(size(G), [n*(d+1) n*(d+1)])
    error('PhasorArray:checkGramCertificate:dimensions', 'M, G and cert are inconsistent.');
end
C = value(M);
if ~isnumeric(C) || ~all(isfinite(C(:)))
    error('PhasorArray:checkGramCertificate:payload', 'M must have finite numeric coefficients.');
end
hermitianDefect = C - permute(conj(flip(C, 3)), [2 1 3]);
if norm(hermitianDefect(:)) > 1e-11*max(1, norm(C(:)))
    error('PhasorArray:checkGramCertificate:notHermitian', 'M(t) must be Hermitian.');
end

% Hermitian Gram in the Fourier basis, then its coefficients (block
% diagonals of H), on the larger of the two degrees.
G = (G + G')/2;
H = B*G*B';
outDegree = max(d, h(M));
reconstructed = zeros(n, n, 2*outDegree + 1);
for row = 0:d
    for col = 0:d
        k = col - row + outDegree + 1;
        reconstructed(:,:,k) = reconstructed(:,:,k) + H(row*n + (1:n), col*n + (1:n));
    end
end
target = zeros(n, n, 2*outDegree + 1);
target(:,:,outDegree - h(M) + 1:outDegree + h(M) + 1) = C;
target(:,:,outDegree + 1) = target(:,:,outDegree + 1) - margin*eye(n);
defect = target - reconstructed;
coefficientResidual = 0;
for k = 1:size(defect, 3)
    coefficientResidual = coefficientResidual + norm(defect(:,:,k), 2);
end

% B is unitary, so V'B G B'V has the same bound (d+1)*lambda_min(G).
gramMinimum = min(eig(G));
lower = margin + (d + 1)*gramMinimum - coefficientResidual;
report = struct('lower', lower, 'accepted', lower > 0, 'gramMinimum', gramMinimum, ...
    'coefficientResidual', coefficientResidual, 'margin', margin, ...
    'degree', d, 'basis', cert.basis);
end
