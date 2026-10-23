function report = certifyHare(A, B, Q, R, P, nvp)
%CERTIFYHARE  A posteriori certification of a harmonic Riccati candidate.
%
%   REPORT = CERTIFYHARE(A, B, Q, R, P) tests whether the real symmetric
%   trigonometric polynomial P is close, in operator norm, to the stabilizing
%   solution of the periodic Riccati equation
%
%       dP/dt + A'P + PA - P S P + Q = 0 ,   S = B R^-1 B' ,
%
%   and returns the error bound when it is. This is step 2 of Algorithm 1
%   (Kleinman-least squares) of Riedinger & Daafouz, "Solving
%   infinite-dimensional harmonic Lyapunov and Riccati equations", corrected
%   version: Proposition 7 and Theorem 8. The candidate typically comes from
%   HARE; the certificate holds whatever the accuracy of the iterations that
%   produced it.
%
%   Procedure, with Ac = A - S P and F the full Riccati residual of P:
%     1. two least-squares Lyapunov solves of order h for the closed loop,
%        dW/dt + Ac'W + W Ac = -I  and  dZ/dt + Ac'Z + Z Ac = -F ;
%     2. rW = sum_k ||RW_k||, with RW the residual of W, and the majorants
%        beta  = sum_k ||W_k|| / (1 - rW)          (bound on the inverse operator)
%        delta = sum_k ||Z_k|| + beta * sum_k ||RZ_k||
%        s     = sum_k ||S_k|| ,  theta = 4 beta s delta ;
%     3. if theta < 1, ||P - P+||_inf <= eps = 2 delta / (1 + sqrt(1 - theta)).
%
%   Closed-loop stability (test (c) of the paper) is established without
%   integrating the monodromy matrix, which is unreliable on stiff systems:
%   W(t) > 0 together with rW < 1 gives dW/dt + Ac'W + W Ac <= -(1 - rW) I,
%   since the supremum of RW(t) is bounded by rW. W is then a Lyapunov
%   function of the closed loop.
%
%   SYNTAX
%     report = certifyHare(A, B, Q, R, P)
%     report = certifyHare(A, B, Q, R, P, T=T, h=h)
%
%   INPUTS
%     A, B   Numeric real PhasorArray, n x n and n x m.
%     Q      Numeric real symmetric PhasorArray, Q(t) > 0.
%     R      Constant symmetric positive definite m x m matrix (double or
%            constant PhasorArray).
%     P      Candidate, numeric real symmetric PhasorArray.
%   Name-value:
%     T      Period (default 2*pi).
%     h      Order of the two auxiliary Lyapunov solves (default: degree of
%            the residual F). Raise it when beta or delta are pessimistic.
%     points Phase samples of the positivity tests (default 8192).
%
%   OUTPUTS
%     report  struct with fields
%       accepted             true when Theorem 8 applies
%       status               "accepted" | "stability_not_established" |
%                            "small_gain_failed" | "nonfinite"
%       errorBound           eps, bound on sup_t ||P(t) - P+(t)|| (NaN if rejected)
%       costGapBound         beta*s*eps^2, bound on the cost excess of K = R^-1 B'P
%       residualBound        sum_k ||F_k||
%       residualOrder        degree of F
%       beta, delta, s, theta, rW, rZ, wLower, auxiliaryOrder
%       stabilityEstablished true when W > 0 and rW < 1
%       lyapunovInfoW, lyapunovInfoZ   info structs of the two LYAP solves
%
%   Bounds are evaluated in floating point; positivity of W and Q is tested on
%   a phase grid widened by an analytic interpolation error. This is not an
%   interval-arithmetic enclosure. The certificate concerns the supplied
%   polynomial model of A, B, Q.
%
%   EXAMPLES
%     [K, P] = hare(A, B, Q, R);
%     report = certifyHare(A, B, Q, R, P);
%     if report.accepted, fprintf('sup|P - P+| <= %.2e\n', report.errorBound), end
%
%   See also: hare, lyap, RicHarmonicKlein.

arguments
    A (:,:) PhasorArray
    B (:,:) PhasorArray
    Q (:,:) PhasorArray
    R
    P (:,:) PhasorArray
    nvp.T (1,1) double {mustBePositive, mustBeFinite} = 2*pi
    nvp.h {mustBeScalarOrEmpty} = []
    nvp.points (1,1) double {mustBeInteger, mustBePositive} = 8192
end

if isa(R, 'PhasorArray')
    if h(R) > 0
        error('PhasorArray:certifyHare:R', 'R must be constant.');
    end
    R = value(R);
end
n = size(A, 1); m = size(B, 2);
% size of a PhasorArray also reports the harmonic dimension: compare dims 1-2.
if size(A, 2) ~= n || size(B, 1) ~= n || size(P, 1) ~= n || size(P, 2) ~= n ...
        || size(Q, 1) ~= n || size(Q, 2) ~= n || ~isequal(size(R), [m m])
    error('PhasorArray:certifyHare:dimensions', 'Incompatible dimensions.');
end
checkRealPolynomial(A, false); checkRealPolynomial(B, false);
checkRealPolynomial(Q, true);  checkRealPolynomial(P, true);
if ~isnumeric(R) || norm(R - R.', 'fro') > 1e-12*max(1, norm(R, 'fro')) ...
        || ~isreal(R) || chol_fails(R)
    error('PhasorArray:certifyHare:R', 'R must be real symmetric positive definite.');
end
if symbolLowerBound(Q, nvp.points) <= 0
    error('PhasorArray:certifyHare:Q', 'Cannot establish Q(t) > 0.');
end

% Full residual and closed loop: all products are exact, nothing is truncated.
S  = B*(R\eye(m))*B.';
F  = P.d(nvp.T) + A.'*P + P*A - P*S*P + Q;
Ac = A - S*P;
hAux = nvp.h;
if isempty(hAux), hAux = h(F); end

lyapArgs = {'T', nvp.T, 'h', hAux, 'autoUpdateh', false, 'systemType', 'rectangle'};
[W, wInfo] = lyap(Ac, PhasorArray(eye(n)), lyapArgs{:});
[Z, zInfo] = lyap(Ac, F, lyapArgs{:});
% Remove round-off skew and imaginary defects of the auxiliary solutions;
% their residuals are recomputed after this projection, so the bounds hold.
W = realSymmetricPart(W); Z = realSymmetricPart(Z);
RW = W.d(nvp.T) + Ac.'*W + W*Ac + eye(n);
RZ = Z.d(nvp.T) + Ac.'*Z + Z*Ac + F;

report = struct('accepted', false, 'status', "stability_not_established", ...
    'errorBound', NaN, 'costGapBound', NaN, ...
    'residualBound', coefficientSum(F), 'residualOrder', h(F), ...
    'beta', NaN, 'delta', NaN, 's', coefficientSum(S), 'theta', NaN, ...
    'rW', coefficientSum(RW), 'rZ', coefficientSum(RZ), ...
    'wLower', symbolLowerBound(W, nvp.points), 'auxiliaryOrder', hAux, ...
    'stabilityEstablished', false, 'lyapunovInfoW', wInfo, 'lyapunovInfoZ', zInfo);
if ~all(isfinite([report.rW report.rZ report.wLower report.s report.residualBound]))
    report.status = "nonfinite"; return
end
if report.wLower <= 0 || report.rW >= 1, return, end
report.stabilityEstablished = true;

% Proposition 7 majorants and Theorem 8 small-gain condition.
report.beta  = coefficientSum(W)/(1 - report.rW);
report.delta = coefficientSum(Z) + report.beta*report.rZ;
report.theta = 4*report.beta*report.s*report.delta;
report.status = "small_gain_failed";
if ~isfinite(report.theta) || report.theta >= 1, return, end
report.errorBound   = 2*report.delta/(1 + sqrt(1 - report.theta));
report.costGapBound = report.beta*report.s*report.errorBound^2;
if ~all(isfinite([report.errorBound report.costGapBound]))
    report.status = "nonfinite"; return
end
report.accepted = true;
report.status = "accepted";
end

function checkRealPolynomial(X, symmetric)
% Finite numeric coefficients, real in time, and symmetric in time if asked.
C = value(X);
if ~isnumeric(C) || ~all(isfinite(C(:)))
    error('PhasorArray:certifyHare:payload', 'Finite numeric coefficients required.');
end
scale = max(1, norm(C(:)));
reality = C - conj(flip(C, 3));
if norm(reality(:)) > 1e-11*scale
    error('PhasorArray:certifyHare:real', 'Real-in-time data required.');
end
if symmetric
    skew = C - permute(C, [2 1 3]);
    if norm(skew(:)) > 1e-11*scale
        error('PhasorArray:certifyHare:symmetric', 'Symmetric-in-time data required.');
    end
end
end

function tf = chol_fails(R)
[~, flag] = chol(R); tf = flag ~= 0;
end

function X = realSymmetricPart(X)
C = value(X); C = (C + conj(flip(C, 3)))/2;
X = PhasorArray((C + permute(C, [2 1 3]))/2);
end

function total = coefficientSum(X)
% sum_k ||X_k||_2, an upper bound on sup_t ||X(t)||_2.
C = value(X); total = 0;
for k = 1:size(C, 3), total = total + norm(C(:,:,k), 2); end
end

function lower = symbolLowerBound(X, points)
% Lower bound of min_t lambda_min(X(t)) for a real symmetric polynomial:
% minimum over a uniform phase grid, minus the interpolation error between
% grid points, bounded by min(pi/N * sum|k|||X_k||, pi^2/(2N^2) * sum k^2||X_k||).
C = value(X); degree = h(X); N = max(points, 16*degree + 1);
slope = 0; curvature = 0;
for k = -degree:degree
    normK = norm(C(:,:,degree+1+k), 2);
    slope = slope + abs(k)*normK; curvature = curvature + k^2*normK;
end
samples = evalp(X, (0:N-1)*2*pi/N); lower = Inf;
for j = 1:N
    page = real(samples(:,:,j));
    lower = min(lower, min(eig((page + page.')/2)));
end
lower = lower - min(pi/N*slope, pi^2/(2*N^2)*curvature);
end
