function [mu0, lambda, quality, ev_m, ev_l] = findTrueFloquetV2(o1, T, varg)
% FINDTRUEFLOQUET_V2  Floquet exponent identification via T0→A-→A→S pipeline (spec v1.3).
%
%   SYNTAX:
%       [mu0, lambda, quality]             = findTrueFloquet_v2(o1, T)
%       [mu0, lambda, quality, ev_m, ev_l] = findTrueFloquet_v2(o1, T, Name, Value)
%
%   INPUTS:
%       o1   - (PhasorArray) Square periodic matrix object.
%       T    - (double)      Period of the system. Default: 2*pi.
%
%   NAME-VALUE:
%       hinit              (int)    Initial harmonic truncation order. Default: o1.h.
%       hmax               (int)    Maximum harmonic truncation order. Default: 10*o1.h.
%       thresholdResidual  (double) Convergence threshold on Liouville trace residual.
%                                   Default: 1e-3.
%       stagnationWindow   (int)    Look-back window for stagnation detection. Default: 5.
%       stagnationRatio    (double) Min relative improvement to avoid stagnation. Default: 0.05.
%       verbose            (int)    Verbosity: 0=silent, 1=table, 2=same, 3=verbose.
%       method             (string) Clustering method:
%                                     'auto'    T0 → DBSCAN → GMM → peigne (default)
%                                     'band'    Band filter only (Tier 0)
%                                     'hdbscan' DBSCAN seed only (A-)
%                                     'gmm_vm'  T0 → DBSCAN → GMM von Mises (A)
%                                     'peigne'  Full pipeline up to comb fit (S)
%       n_multistart       (int)    Peigne multi-start runs (GUD-004). Default: 20.
%       card_thresh        (double) Fusion detection cardinality factor. Default: 1.5.
%       tol_liouville      (double) Liouville residual tolerance. Default: 1e-3.
%       alpha_sep          (double) FM-1 split acceptance: sep > alpha*sigma. Default: 2.0.
%
%   OUTPUTS:
%       mu0     - (n×1 complex) Principal Floquet exponents.
%       lambda  - (n×1 complex) Floquet multipliers exp(T*mu0).
%       quality - (struct) Full diagnostics:
%                   Spec v1.3 fields:
%                     .sigma          : n×1 intra-cluster dispersion
%                     .liouville_res  : |sum(mu0) - trA0| after WLS projection
%                     .multiplicity   : n×1 detected algebraic multiplicity
%                     .outlier_frac   : fraction of M_h rejected as outliers
%                     .method_used    : string — method that produced the result
%                     .fm2_suspected  : logical — band filter failed (Im(mu)~±ω/2)
%                     .n_hat          : integer — coherence check (should equal n)
%                     .fallback_skipped : logical — n_hat≠n, FM-1 fallback skipped
%                     .fallback_reason  : string — explanation if fallback_skipped
%                     .low_confidence   : logical — FM-3 indicators active
%                   Convergence fields (from outer h-loop):
%                     .status         : 0=converged, 1=stagnated, 2=hmax_reached
%                     .statusMsg      : human-readable exit summary
%                     .hhistory       : truncation order at each iteration
%                     .resTrace       : Liouville residual at each iteration
%                     .hBest          : h at which best residual was achieved
%                     .resBest        : best residual achieved
%       ev_m    - (vector) Full HSS spectrum — Floquet exponents (optional).
%       ev_l    - (vector) Full HSS spectrum — multipliers exp(T*ev_m) (optional).
%
%   ALGORITHM:
%       Implements spec-algorithm-floquet-exponent-identification.md v1.3.
%       Outer h-loop (Floquet-Liouville trace identity as convergence criterion) drives
%       an inner clustering pipeline at each h:
%         Tier 0 (T0) : Band filter — select ν with |Im(ν)| < ω/2. Exact if FM-2 absent.
%         Rank A-     : DBSCAN with cylindrical metric + k̂ pre-weighting. Seed for A/S.
%         Rank A      : GMM von Mises-Gaussian EM + noise component (Rang A, spec §8).
%         Rank S      : Parametric comb fit, k̂-aware, multi-start (GUD-004).
%       All paths enforce Liouville WLS projection (REQ-004, spec §6).
%
%   REFERENCE:
%       spec-algorithm-floquet-exponent-identification.md v1.3
%       pseudocode-floquet-identification.md v1.3
%       audit-clustering-floquet.md v1.3
%
%   See also: findTrueFloquet, HmqNEig.

arguments
    o1    PhasorArray
    T    (1,1) double                               = 2*pi
    varg.hinit             (1,1) {mustBeInteger, mustBePositive} = o1.h
    varg.hmax              (1,1) {mustBeInteger, mustBePositive} = o1.h * 10
    varg.thresholdResidual (1,1) double             = 1e-3
    varg.stagnationWindow  (1,1) {mustBeInteger, mustBePositive} = 5
    varg.stagnationRatio   (1,1) double             = 0.05
    varg.verbose           (1,1) {mustBeInteger, mustBeNonnegative} = 0
    varg.method            {mustBeMember(varg.method, {'auto','band','hdbscan','gmm_vm','peigne'})} = 'auto'
    varg.n_multistart      (1,1) {mustBeInteger, mustBePositive} = 20
    varg.card_thresh       (1,1) double             = 1.5
    varg.tol_liouville     (1,1) double             = 1e-3
    varg.alpha_sep         (1,1) double             = 2.0
end
nx = o1.size(1);
assert(o1.size(1) == o1.size(2), 'PhasorArray must be square for eigenvalue analysis.');
assert(varg.hmax >= varg.hinit, ...
    'findTrueFloquet_v2: hmax (%d) must be >= hinit (%d).', varg.hmax, varg.hinit);

capacity = varg.hmax - varg.hinit + 1;

%% --- Trace of A0 ---
trA0            = trace(o1.phas(0));
useAbsCriterion = abs(trA0) < 1e-12 * nx;
if useAbsCriterion && varg.verbose >= 1
    warning('findTrueFloquet_v2:zeroTrace', ...
        'tr(A0) ≈ 0 (|tr|=%.2e): switching to absolute trace criterion.', abs(trA0))
end

%% --- Options struct for inner pipeline ---
opts.method        = varg.method;
opts.n_multistart  = varg.n_multistart;
opts.card_thresh   = varg.card_thresh;
opts.tol_liouville = varg.tol_liouville;
opts.alpha_sep     = varg.alpha_sep;

%% --- Verbose header ---
if varg.verbose >= 1
    rule = repmat('=', 1, 72);
    fprintf('\n%s\n', rule)
    fprintf('  findTrueFloquet_v2  --  spec v1.3: T0→A-→A→S pipeline\n')
    fprintf('%s\n', rule)
    fprintf('  Method     : %s\n', varg.method)
    if useAbsCriterion
        fprintf('  Criterion  : |sum(mu0)-tr(A0)| < thr  [absolute — tr(A0)~0]\n')
    else
        fprintf('  Criterion  : |sum(mu0)-tr(A0)| / |tr(A0)| < thr  [Floquet-Liouville]\n')
    end
    fprintf('%s\n', repmat('-', 1, 72))
    fprintf('  Problem\n')
    fprintf('    State dim  nx               = %d\n',   nx)
    fprintf('    Period     T                = %.6g\n', T)
    fprintf('    tr(A0)                      = %.6g%+.6gi\n', real(trA0), imag(trA0))
    fprintf('    hinit                       = %d\n',   varg.hinit)
    fprintf('    hmax                        = %d\n',   varg.hmax)
    fprintf('    thresholdResidual           = %.2e\n', varg.thresholdResidual)
    fprintf('%s\n', repmat('-', 1, 72))
    fprintf('  Stagnation\n')
    fprintf('    stagnationWindow            = %d\n',   varg.stagnationWindow)
    fprintf('    stagnationRatio             = %.2f   (min relative improvement)\n', varg.stagnationRatio)
    fprintf('%s\n', rule)
end
if varg.verbose == 1 || varg.verbose == 2
    if useAbsCriterion
        resLabel = 'traceResAbs';
    else
        resLabel = 'traceRes   ';
    end
    fprintf(' %4s | %3s | %11s | %10s | %8s | %9s | %s\n', ...
        'iter', 'h', resLabel, 'method', 'step(s)', 'total(s)', 'Note')
    fprintf('------|-----|-------------|------------|----------|-----------|------\n')
end

t_total = tic;

%% --- Pre-allocate history ---
hHist  = zeros(1, capacity);
reHist = zeros(1, capacity);

%% --- Initial computation (h = hinit) ---
h    = varg.hinit;
ev_m = o1.HmqNEig(h, T);
ev_l = exp(T * ev_m);

t_step = tic;
[mu0, lambda, quality_cur] = identifyFloquet(ev_m, nx, T, trA0, h, opts);
dt0 = toc(t_step);

resCur    = computeResidual(mu0, trA0, useAbsCriterion);
nIter     = 1;
hHist(1)  = h;
reHist(1) = resCur;

% Best-solution tracking
mu0_best     = mu0;
lambda_best  = lambda;
ev_m_best    = ev_m;
ev_l_best    = ev_l;
quality_best = quality_cur;
resBest      = resCur;
hBest        = h;

note = '';
if resCur <= varg.thresholdResidual, note = 'converged'; end

if varg.verbose == 1 || varg.verbose == 2
    fprintf(' %4d | %3d | %11.4e | %10s | %8.3f | %9.3f | %s\n', ...
        1, h, resCur, quality_cur.method_used, dt0, toc(t_total), note)
elseif varg.verbose >= 3
    fprintf('\n%s\n  iter 1  (h=%d)  method=%s\n%s\n', ...
        repmat('-', 1, 60), h, quality_cur.method_used, repmat('-', 1, 60))
    fprintf('  traceRes=%.4e | FM2=%d | low_conf=%d\n', ...
        resCur, quality_cur.fm2_suspected, quality_cur.low_confidence)
end

status = -1;
if resCur <= varg.thresholdResidual, status = 0; end

%% --- Iterative h-refinement ---
while status == -1 && h < varg.hmax
    t_step = tic;
    h      = h + 1;
    nIter  = nIter + 1;

    ev_m = o1.HmqNEig(h, T);
    ev_l = exp(T * ev_m);
    [mu0, lambda, quality_cur] = identifyFloquet(ev_m, nx, T, trA0, h, opts);

    resCur        = computeResidual(mu0, trA0, useAbsCriterion);
    hHist(nIter)  = h;
    reHist(nIter) = resCur;

    % Best-solution tracking
    if resCur < resBest
        mu0_best     = mu0;
        lambda_best  = lambda;
        ev_m_best    = ev_m;
        ev_l_best    = ev_l;
        quality_best = quality_cur;
        resBest      = resCur;
        hBest        = h;
    end

    note = '';

    % Convergence
    if resCur <= varg.thresholdResidual
        status = 0;
        note   = 'converged';
    end

    % Stagnation (only after filling window)
    if status == -1 && nIter >= varg.stagnationWindow
        win       = reHist(nIter - varg.stagnationWindow + 1 : nIter);
        relImprov = (win(1) - min(win)) / (win(1) + eps);
        if relImprov < varg.stagnationRatio
            status = 1;
            note   = sprintf('stagnated (%.1f%%)', relImprov * 100);
        end
    end

    if varg.verbose == 1 || varg.verbose == 2
        fprintf(' %4d | %3d | %11.4e | %10s | %8.3f | %9.3f | %s\n', ...
            nIter, h, resCur, quality_cur.method_used, toc(t_step), toc(t_total), note)
    elseif varg.verbose >= 3
        fprintf('\n%s\n  iter %d  (h=%d)  method=%s\n%s\n', ...
            repmat('-', 1, 60), nIter, h, quality_cur.method_used, repmat('-', 1, 60))
        fprintf('  traceRes=%.4e | step=%.3f s | total=%.3f s | %s\n', ...
            resCur, toc(t_step), toc(t_total), note)
    end
end

%% --- hmax exit ---
if status == -1
    status = 2;
end

%% --- Return best solution on non-converged exits ---
if status ~= 0
    mu0     = mu0_best;
    lambda  = lambda_best;
    ev_m    = ev_m_best;
    ev_l    = ev_l_best;
    quality = quality_best;
else
    quality = quality_cur;
end

%% --- Pack convergence fields into quality ---
hHist  = hHist(1:nIter);
reHist = reHist(1:nIter);

quality.status    = status;
quality.statusMsg = buildStatusMsg(status, h, varg.hmax, reHist(end), ...
    varg.thresholdResidual, nIter, varg.stagnationWindow, hBest, resBest);
quality.hhistory  = hHist;
quality.resTrace  = reHist;
quality.hBest     = hBest;
quality.resBest   = resBest;

%% --- Verbose footer ---
if varg.verbose >= 1
    rule = repmat('=', 1, 72);
    fprintf('\n  %s\n', quality.statusMsg)
    fprintf('%s\n\n', rule)
end

end  % findTrueFloquet_v2

%% =========================================================================
%% ORCHESTRATEUR INTERNE (spec §9)
%% =========================================================================

function [mu0, lambda, quality] = identifyFloquet(M_h, n, T, trA0, h, opts)
% identifyFloquet  Single-h clustering pipeline: T0 → A- → A → S (spec §9).
%   M_h    : N×1 complex — truncated Hill spectrum (N = n*(2h+1))
%   n      : system dimension
%   T      : period
%   trA0   : trace of zeroth Fourier coefficient of A(t)
%   h      : current truncation order
%   opts   : options struct

omega = 2*pi / T;

%% --- Initialize quality struct (all fields, spec §8) ---
quality.sigma           = NaN(n, 1);
quality.liouville_res   = Inf;
quality.multiplicity    = ones(n, 1);
quality.outlier_frac    = 0;
quality.method_used     = opts.method;
quality.fm2_suspected   = false;
quality.n_hat           = n;
quality.fallback_skipped = false;
quality.fallback_reason  = '';
quality.low_confidence  = false;

M_h = M_h(:);  % ensure column vector

%% --- Module 0: k̂ reconstruction and weights ---
[k_hat, weights] = computeKHat(M_h, omega);

%% --- Module 1: Tier 0 — band filter ---
if ismember(opts.method, {'auto', 'band'})
    [mu0_band, t0_ok] = bandFilter(M_h, n, omega, trA0, opts.tol_liouville);

    if t0_ok
        sigmas_band = estimateBandSigma(mu0_band, M_h, k_hat, omega, T);
        mu0         = liouvilleWLS(mu0_band, sigmas_band, ones(n,1), trA0);
        lambda      = exp(T * mu0);

        quality.method_used    = 'band';
        quality.sigma          = sigmas_band;
        quality.liouville_res  = abs(sum(mu0) - trA0);
        quality.outlier_frac   = 1 - n / numel(M_h);
        quality.multiplicity   = ones(n, 1);
        if quality.liouville_res > opts.tol_liouville
            quality.low_confidence = true;
        end
        return
    end

    if strcmp(opts.method, 'band')
        error('findTrueFloquet_v2:bandFailed', ...
            'Band filter failed at h=%d — FM-2 suspected. Use method=''auto''.', h)
    end

    quality.fm2_suspected = true;  % T0 failed → FM-2 suspected
end

%% --- Module 3: wrap dataset (AFTER k̂ is computed) ---
M_wrapped = wrapDataset(M_h, omega);

%% --- Module 4: DBSCAN seed (Rank A-) ---
[centroids, sigmas, labels, outlier_frac] = ...
    seedClustering(M_wrapped, weights, n, omega, T);

%% --- Module 5: FM-1 detection ---
[~, centroids, sigmas, multiplicities, n_hat, fallback_skipped, fallback_reason] = ...
    detectFusion(M_wrapped, labels, centroids, sigmas, n, h, outlier_frac, omega, T, opts);

quality.n_hat            = n_hat;
quality.fallback_skipped = fallback_skipped;
quality.fallback_reason  = fallback_reason;
if fallback_skipped
    quality.low_confidence = true;
end

%% --- Module 6: WLS Liouville projection (post-seed) ---
centroids = liouvilleWLS(centroids, sigmas, multiplicities, trA0);

%% --- HDBSCAN-only exit ---
if strcmp(opts.method, 'hdbscan')
    mu0 = centroids;
    quality.sigma         = sigmas;
    quality.liouville_res = abs(sum(mu0) - trA0);
    quality.outlier_frac  = outlier_frac;
    quality.multiplicity  = multiplicities;
    quality.method_used   = 'hdbscan';
    lambda = exp(T * mu0);
    if quality.liouville_res > opts.tol_liouville
        quality.low_confidence = true;
    end
    return
end

%% --- Module 7: GMM von Mises-Gaussian (Rank A) ---
if ismember(opts.method, {'auto', 'gmm_vm'})
    [mu0_gmm, sigmas_gmm] = gmmVonMises(M_wrapped, weights, n, T, trA0, centroids);

    if strcmp(opts.method, 'gmm_vm')
        mu0 = mu0_gmm;
        quality.sigma         = sigmas_gmm;
        quality.liouville_res = abs(sum(mu0) - trA0);
        quality.outlier_frac  = outlier_frac;
        quality.multiplicity  = multiplicities;
        quality.method_used   = 'gmm_vm';
        lambda = exp(T * mu0);
        if quality.liouville_res > opts.tol_liouville
            quality.low_confidence = true;
        end
        return
    end

    % Update seed for peigne fit
    centroids = mu0_gmm;
    sigmas    = sigmas_gmm;
end

%% --- Module 8: Parametric comb fit, multi-start (Rank S) ---
if ismember(opts.method, {'auto', 'peigne'})
    [mu0, quality_peigne] = peigFitMultistart(M_h, k_hat, weights, n, omega, T, ...
        trA0, centroids, sigmas, multiplicities, opts);

    quality.sigma         = quality_peigne.sigma;
    quality.liouville_res = quality_peigne.liouville_res;
    quality.outlier_frac  = quality_peigne.outlier_frac;
    quality.multiplicity  = multiplicities;
    quality.method_used   = 'peigne';
    lambda = exp(T * mu0);
end

%% --- FM-3 check ---
if quality.liouville_res > opts.tol_liouville
    quality.low_confidence = true;
end

end  % identifyFloquet

%% =========================================================================
%% MODULE 0 — k̂ RECONSTRUCTION AND WEIGHTS
%% =========================================================================

function [k_hat, weights] = computeKHat(M_h, omega)
% computeKHat  Reconstruct harmonic index k̂ from Im(ν) BEFORE wrapping.
%   k̂(ν) = round(Im(ν)/ω).  Accuracy degrades near |k̂| ~ h (trunc. boundary).
%   Weights: w(ν) = 1/(1+|k̂(ν)|) — depreciate poorly-converged aliases.

k_hat   = round(imag(M_h) ./ omega);
weights = 1 ./ (1 + abs(k_hat));

end

%% =========================================================================
%% MODULE 1 — TIER 0: BAND FILTER
%% =========================================================================

function [mu0_band, converged] = bandFilter(M_h, n, omega, trA0, tol_liouville)
% bandFilter  Select ν with |Im(ν)| < ω/2 and check cardinality + Liouville.
%   In M_∞, the k=0 band contains exactly n true Floquet exponents.
%   Two independent convergence checks: cardinality AND Liouville on Re.

band = M_h(abs(imag(M_h)) < omega/2);

if numel(band) ~= n
    mu0_band  = [];
    converged = false;
    return
end

residual = abs(sum(real(band)) - real(trA0));
if residual > tol_liouville
    mu0_band  = [];
    converged = false;
    return
end

mu0_band  = band;
converged = true;

end

%% =========================================================================
%% MODULE 1 — BAND SIGMA ESTIMATION
%% =========================================================================

function sigmas = estimateBandSigma(mu0_band, M_h, k_hat, omega, ~)
% estimateBandSigma  Intra-cluster dispersion for band filter path.
%   Folds all M_h back via k̂: nu_res = Re(ν) + j*wrap(Im(ν) - k̂·ω).
%   Assigns each nu_res to nearest band exponent and computes std(d_cyl).

n      = numel(mu0_band);
N      = numel(M_h);
nu_res = real(M_h) + 1j .* wrapToHalf(imag(M_h) - k_hat .* omega, omega);

% Pairwise cylindrical distances: N×n
D_ni = zeros(N, n);
for i = 1:n
    D_ni(:,i) = dCyl(nu_res, mu0_band(i), omega);
end

[~, assign] = min(D_ni, [], 2);

sigmas = zeros(n, 1);
for i = 1:n
    idx_i = assign == i;
    if ~any(idx_i)
        sigmas(i) = 1e-8;
    else
        sigmas(i) = max(std(D_ni(idx_i, i)), 1e-8);
    end
end

end

%% =========================================================================
%% MODULE 2 — CYLINDRICAL GEOMETRY UTILITIES
%% =========================================================================

function out = wrapToHalf(theta, omega)
% wrapToHalf  Wrap to [-omega/2, omega/2). Vectorized.
out = theta - omega .* floor(theta ./ omega + 0.5);
end

function d = dCirc(a, b, omega)
% dCirc  Circular distance on S^1 of period omega. Vectorized.
delta = abs(a - b);
d     = min(delta, omega - delta);
end

function d = dCyl(a, b, omega)
% dCyl  Cylindrical metric: Euclidean on Re, circular on Im. Vectorized.
%   a, b can be vectors (elementwise) or one scalar + one vector.
d = sqrt((real(a) - real(b)).^2 + dCirc(imag(a), imag(b), omega).^2);
end

function theta_mean = circularMean(thetas, weights, T)
% circularMean  Weighted Fréchet mean on S^1 via complex exponential (CON-004).
%   thetas in [-omega/2, omega/2), omega = 2pi/T.
%   Returns arithmetic mean = 0 if distribution is uniform.
if sum(weights) < 1e-14
    theta_mean = 0;
    return
end
z = sum(weights(:) .* exp(1j * thetas(:) * T)) / sum(weights);
if abs(z) < 1e-14
    theta_mean = 0;
else
    theta_mean = angle(z) / T;
end
end

function M_wrapped = wrapDataset(M_h, omega)
% wrapDataset  Apply Im → [-omega/2, omega/2) wrapping to M_h.
%   MUST be called AFTER computeKHat — k̂ requires unwrapped Im(ν).
M_wrapped = real(M_h) + 1j .* wrapToHalf(imag(M_h), omega);
end

%% =========================================================================
%% MODULE 4 — DBSCAN SEED (Rank A-)
%% =========================================================================

function [centroids, sigmas, labels, outlier_frac] = ...
    seedClustering(M_wrapped, weights, n, omega, T)
% seedClustering  DBSCAN with cylindrical metric, forced to exactly n clusters.
%   Approximation of HDBSCAN (not available natively in MATLAB).
%   Auto-epsilon from k-NN distance median. Forces exactly n clusters
%   by merging surplus and splitting deficit.

N      = numel(M_wrapped);
minPts = max(3, floor(0.3 * N / n));

%% Precomputed N×N cylindrical distance matrix (O(N²), acceptable for N≤500)
re_w = real(M_wrapped(:));
im_w = imag(M_wrapped(:));
dRe  = abs(bsxfun(@minus, re_w, re_w'));       % N×N
dIm  = abs(bsxfun(@minus, im_w, im_w'));
dIm  = min(dIm, omega - dIm);
D    = sqrt(dRe.^2 + dIm.^2);

% k̂-weight adjustment: depreciate pairs involving high-|k̂| points
w_mat = sqrt(max(bsxfun(@times, weights(:), weights(:)'), 1e-12));
D     = D ./ w_mat;
D(1:N+1:end) = 0;  % zero diagonal

%% Auto epsilon: median k-NN distance
k_nn  = min(minPts, N - 1);
Dsort = sort(D, 2);
eps_auto = max(median(Dsort(:, k_nn + 1)), 1e-8);

%% DBSCAN (requires Statistics and Machine Learning Toolbox, R2019a+)
if exist('dbscan', 'builtin') || exist('dbscan', 'file')
    labels = dbscan(D, eps_auto, minPts, 'Distance', 'precomputed');
else
    % Fallback: k-means cylindrical with Gonzalez init
    warning('findTrueFloquet_v2:noDBSCAN', ...
        'dbscan not found — using cylindrical k-means as A- fallback.')
    seeds  = gonzalezInit(M_wrapped, n, omega, T);
    labels = kmeansCyl(M_wrapped, n, seeds, omega, T);
    outlier_frac = 0;
    [centroids, sigmas] = clusterStats(M_wrapped, labels, n, omega, T);
    return
end

%% Force exactly n clusters
labels = forceNClusters(labels, M_wrapped, D, n, omega, T);

%% Compute centroids and intra-cluster dispersions
n_cl         = max(labels);
n_cl         = max(n_cl, n);
outlier_frac = sum(labels == 0) / N;  % DBSCAN uses 0 for noise in MATLAB

[centroids, sigmas] = clusterStats(M_wrapped, labels, n_cl, omega, T);

end

% -------------------------------------------------------------------------

function [centroids, sigmas] = clusterStats(M_wrapped, labels, n_cl, omega, T)
% clusterStats  Compute centroids (circularMean on Im) and sigma (d_cyl std).
centroids = zeros(n_cl, 1);
sigmas    = zeros(n_cl, 1);
for i = 1:n_cl
    idx_i = labels == i;
    if ~any(idx_i)
        centroids(i) = 0;
        sigmas(i)    = 1e-8;
        continue
    end
    pts_i     = M_wrapped(idx_i);
    w_i       = ones(sum(idx_i), 1);  % uniform weights at seed stage
    c_re      = mean(real(pts_i));
    c_im      = circularMean(imag(pts_i), w_i, T);
    centroids(i) = c_re + 1j * c_im;
    d_i = dCyl(pts_i, centroids(i), omega);
    sigmas(i) = max(std(d_i), 1e-8);
end
end

% -------------------------------------------------------------------------

function labels = forceNClusters(labels_in, M_wrapped, D, n_target, omega, T)
% forceNClusters  Adjust DBSCAN output to exactly n_target clusters.
%   Noise (label 0) → assigned to nearest cluster.
%   Surplus clusters → merge smallest to nearest.
%   Deficit clusters → split largest.

labels = labels_in;

% Step 1: Assign noise points (label 0) to nearest cluster
noise_idx = find(labels == 0);
if ~isempty(noise_idx) && max(labels) > 0
    for k = noise_idx'
        d_to_cls = zeros(max(labels), 1);
        for c = 1:max(labels)
            cl_idx = find(labels == c);
            if isempty(cl_idx)
                d_to_cls(c) = Inf;
            else
                d_to_cls(c) = mean(D(k, cl_idx));
            end
        end
        [~, best_c] = min(d_to_cls);
        labels(k) = best_c;
    end
end
if max(labels(labels > 0), [], 'all', 'omitnan') == 0 || isempty(labels(labels > 0))
    % No clusters found at all — use Gonzalez k-means as full fallback
    seeds  = gonzalezInit(M_wrapped, n_target, omega, T);
    labels = kmeansCyl(M_wrapped, n_target, seeds, omega, T);
    return
end

% Compact labels (remove gaps)
labels = renumberLabels(labels);
n_cl   = max(labels);

% Step 2: Too many clusters → merge smallest to nearest
while n_cl > n_target
    sizes  = arrayfun(@(c) sum(labels == c), 1:n_cl);
    [~, smallest] = min(sizes);
    mask_s = labels == smallest;
    best_other = 0;
    best_dist  = Inf;
    for c = 1:n_cl
        if c == smallest, continue; end
        d_sc  = mean(mean(D(mask_s, labels == c)));
        if d_sc < best_dist
            best_dist  = d_sc;
            best_other = c;
        end
    end
    labels(labels == smallest) = best_other;
    labels = renumberLabels(labels);
    n_cl   = max(labels);
end

% Step 3: Too few clusters → split largest
while n_cl < n_target
    sizes  = arrayfun(@(c) sum(labels == c), 1:n_cl);
    [~, largest] = max(sizes);
    idx_l  = find(labels == largest);
    pts_l  = M_wrapped(idx_l);
    seeds  = gonzalezInit(pts_l, 2, omega, T);
    sub_lb = kmeansCyl(pts_l, 2, seeds, omega, T);
    labels(idx_l(sub_lb == 2)) = n_cl + 1;
    n_cl   = max(labels);
end

end

% -------------------------------------------------------------------------

function labels = renumberLabels(labels)
% renumberLabels  Compact positive labels to 1..k without gaps.
unique_labels = unique(labels(labels > 0));
new_labels    = labels;
for i = 1:numel(unique_labels)
    new_labels(labels == unique_labels(i)) = i;
end
labels = new_labels;
end

%% =========================================================================
%% MODULE 5 — FM-1: FUSION DETECTION
%% =========================================================================

function [labels, centroids, sigmas, multiplicities, n_hat, fallback_skipped, fallback_reason] = ...
    detectFusion(M_wrapped, labels_in, centroids_in, sigmas_in, n, h, outlier_frac, omega, T, opts)
% detectFusion  Detect and attempt to split fused clusters (FM-1, spec §6).
%   Pre-test: n_hat = sum(p_est) must equal n.
%   Split method: Gonzalez init + cylindrical k-means.

expected_card   = 2*h + 1;
card_correction = 1 / max(1 - outlier_frac, 0.1);
n_c             = max(labels_in);
if isempty(n_c) || n_c == 0
    n_c = n;
end

labels     = labels_in;
centroids  = centroids_in;
sigmas     = sigmas_in;
multiplicities = ones(n_c, 1);

% Effective cardinalities and p_est
card_eff = zeros(n_c, 1);
for i = 1:n_c
    card_eff(i) = sum(labels == i) * card_correction;
end
p_est = max(1, round(card_eff ./ expected_card));

% Pre-test: coherence check
n_hat = sum(p_est);
if n_hat ~= n
    fallback_skipped = true;
    fallback_reason  = sprintf('n_hat=%d ~= n=%d | outlier_frac=%.3f | h=%d → augmenter h', ...
        n_hat, n, outlier_frac, h);
    return
end

fallback_skipped = false;
fallback_reason  = '';

% Process fused clusters (largest first)
[~, order] = sort(card_eff, 'descend');
n_current  = n_c;

for ii = 1:numel(order)
    i = order(ii);
    if n_current >= n, break; end

    p = p_est(i);
    if p < 2, continue; end

    pts_i = M_wrapped(labels == i);
    if numel(pts_i) < 2*p
        multiplicities(i) = p;
        continue
    end

    % Gonzalez init + cylindrical k-means split into p sub-clusters
    seeds  = gonzalezInit(pts_i, p, omega, T);
    sub_lb = kmeansCyl(pts_i, p, seeds, omega, T);

    % Sub-centroids
    sub_c = zeros(p, 1);
    for s = 1:p
        pts_s = pts_i(sub_lb == s);
        if isempty(pts_s)
            sub_c(s) = centroids(i);
        else
            sub_c(s) = mean(real(pts_s)) + 1j * circularMean(imag(pts_s), ones(numel(pts_s),1), T);
        end
    end

    % Minimum separation between sub-centroids
    sep_min = Inf;
    for a = 1:p
        for b = a+1:p
            sep_min = min(sep_min, dCyl(sub_c(a), sub_c(b), omega));
        end
    end

    if sep_min > opts.alpha_sep * sigmas(i)
        % Accept split
        orig_idx = find(labels == i);
        for s = 1:p
            cl  = (s == 1) * i + (s > 1) * (n_current + s - 1);
            in_s = orig_idx(sub_lb == s);
            labels(in_s) = cl;
            pts_s = M_wrapped(in_s);
            if isempty(pts_s)
                centroids(cl) = sub_c(s);
                sigmas(cl)    = 1e-8;
            else
                centroids(cl) = mean(real(pts_s)) + 1j * circularMean(imag(pts_s), ones(numel(pts_s),1), T);
                d_s = dCyl(pts_s, centroids(cl), omega);
                sigmas(cl) = max(std(d_s), 1e-8);
            end
        end
        n_current    = n_current + (p - 1);
        multiplicities = [multiplicities; ones(p-1, 1)]; %#ok<AGROW>
    else
        multiplicities(i) = p;
    end
end

% Trim to n_current
multiplicities = multiplicities(1:n_current);

end

%% =========================================================================
%% GONZALEZ INIT AND CYLINDRICAL K-MEANS
%% =========================================================================

function seeds = gonzalezInit(points, p, omega, T)
% gonzalezInit  Gonzalez (1985) initialization: p maximally-separated seeds.
%   O(p·N). Robust against near-degenerate configurations.

points = points(:);
N = numel(points);

% First seed: farthest from cloud centroid
c_re = mean(real(points));
c_im = circularMean(imag(points), ones(N,1), T);
centroid_g = c_re + 1j * c_im;

d_to_g    = dCyl(points, centroid_g, omega);
[~, idx1] = max(d_to_g);
seeds     = points(idx1);

for k = 2:p
    d_nearest = Inf(N, 1);
    for s = 1:numel(seeds)
        d_s = dCyl(points, seeds(s), omega);
        d_nearest = min(d_nearest, d_s);
    end
    [~, next_idx] = max(d_nearest);
    seeds(end+1) = points(next_idx); %#ok<AGROW>
end

seeds = seeds(:);

end

% -------------------------------------------------------------------------

function labels = kmeansCyl(points, k, seeds, omega, T)
% kmeansCyl  K-means with cylindrical metric and circular mean on Im.
%   Reinitializes empty clusters from farthest remaining point.

MAX_ITER  = 100;
points    = points(:);
N         = numel(points);
centroids = seeds(:);
labels    = ones(N, 1);

for iter = 1:MAX_ITER
    labels_prev = labels;

    % Assignment: N×k distance matrix
    D_nk = zeros(N, k);
    for c = 1:k
        D_nk(:, c) = dCyl(points, centroids(c), omega);
    end
    [~, labels] = min(D_nk, [], 2);

    % Update
    for c = 1:k
        idx_c = labels == c;
        if ~any(idx_c)
            % Reinitialize: farthest point from any centroid
            min_d = min(D_nk, [], 2);
            [~, far_idx] = max(min_d);
            centroids(c) = points(far_idx);
        else
            pts_c    = points(idx_c);
            w_c      = ones(sum(idx_c), 1);
            centroids(c) = mean(real(pts_c)) + 1j * circularMean(imag(pts_c), w_c, T);
        end
    end

    if isequal(labels, labels_prev), break; end
end

end

%% =========================================================================
%% MODULE 6 — WLS LIOUVILLE PROJECTION (REQ-004)
%% =========================================================================

function mu0_corr = liouvilleWLS(mu0_raw, sigmas, multiplicities, trA0)
% liouvilleWLS  Variance-weighted projection of Re(mu0) onto Liouville constraint.
%   Distributes residual proportionally to sigma_i^2 * multiplicity_i (spec §6).
%   Im(mu0) is unchanged — no branch ambiguity on Re.
%
%   Guard: if sigma2_total < 1e-14 (perfectly converged clusters),
%          use uniform correction epsilon/n instead.

mu0_raw       = mu0_raw(:);
sigmas        = sigmas(:);
multiplicities = multiplicities(:);
n             = numel(mu0_raw);

epsilon      = sum(real(mu0_raw)) - real(trA0);
sigma2_total = sum(sigmas.^2 .* multiplicities);

if sigma2_total < 1e-14
    correction = (epsilon / n) * ones(n, 1);
else
    correction = sigmas.^2 .* multiplicities .* epsilon ./ sigma2_total;
end

mu0_corr = mu0_raw - correction;  % only Re changes (correction is real)

end

%% =========================================================================
%% MODULE 7 — GMM VON MISES-GAUSSIAN (Rank A)
%% =========================================================================

function [mu0, sigmas] = gmmVonMises(M_wrapped, weights, n, T, trA0, centroids_init)
% gmmVonMises  EM on cylindrical model: Gaussian(Re) × von Mises(Im) + noise.
%   Geometrically correct on cylinder ℝ×S¹ (Rank A, audit §8).
%   Numerically stable via scaled Bessel: besseli(0,κ,1) = exp(-κ)*I₀(κ).
%   WLS Liouville projection applied post-EM.

MAX_ITER = 200;
TOL_LOGLIK = 1e-6;
omega   = 2*pi / T;
N       = numel(M_wrapped);
K       = n + 1;  % n signal + 1 noise

re_data = real(M_wrapped(:));
im_data = imag(M_wrapped(:));
w       = weights(:);
W_total = sum(w);

% --- Initialization ---
pi_mix   = ones(K, 1) / K;
m_R      = real(centroids_init(:));
sigma_R  = max(0.1 * omega * ones(n,1), 1e-6);
theta    = imag(centroids_init(:));
kappa    = ones(n, 1);

range_Re    = max(max(re_data) - min(re_data), 1e-6);
noise_dens  = 1 / (omega * range_Re);  % uniform on cylinder

log_lik_prev = -Inf;
r = zeros(N, K);

for iter = 1:MAX_ITER
    %% == E-step ==
    for i = 1:n
        % Gaussian on Re — log form then exp for stability
        z_R     = (re_data - m_R(i)) ./ sigma_R(i);
        log_gR  = -0.5 * z_R.^2 - log(sigma_R(i)) - 0.5*log(2*pi);

        % von Mises on Im — numerically stable via cos-1 form (spec §8 note)
        %   VM = exp(κ·cos(θ-θ₀)) / (2π·I₀(κ))
        %      = exp(κ·(cos(θ-θ₀)-1)) / (2π·besseli(0,κ,1))
        %   since besseli(0,κ,1) = exp(-κ)·I₀(κ)
        cos_arg  = (im_data - theta(i)) * T;
        I0_sc    = besseli(0, kappa(i), 1);  % scaled: exp(-κ)·I₀(κ)
        log_gIm  = kappa(i) * (cos(cos_arg) - 1) - log(2*pi) - log(max(I0_sc, 1e-300));

        r(:, i) = pi_mix(i) * exp(log_gR + log_gIm) .* w;
    end
    r(:, K) = pi_mix(K) * noise_dens .* w;

    r_sum = sum(r, 2);
    log_lik = sum(log(max(r_sum, 1e-300)));  % weighted log-likelihood
    r_sum = max(r_sum, 1e-300);
    r     = r ./ r_sum;

    %% == M-step ==
    for i = 1:n
        Ni         = max(sum(r(:,i)), 1e-10);
        pi_mix(i)  = Ni / W_total;
        m_R(i)     = sum(r(:,i) .* re_data) / Ni;
        var_R      = sum(r(:,i) .* (re_data - m_R(i)).^2) / Ni;
        sigma_R(i) = max(sqrt(var_R), 1e-6);

        % Circular mean for Im
        z_i      = sum(r(:,i) .* exp(1j * im_data * T)) / Ni;
        theta(i) = angle(z_i) / T;
        R_bar_i  = min(abs(z_i), 1 - 1e-6);
        kappa(i) = updateKappa(R_bar_i);
    end
    pi_mix(K) = max(0, 1 - sum(pi_mix(1:n)));

    if abs(log_lik - log_lik_prev) < TOL_LOGLIK, break; end
    log_lik_prev = log_lik;
end

%% --- Extract results ---
mu0_raw = m_R + 1j * theta;

% Combined Re + Im uncertainty per component (spec §8)
sigmas = sqrt(sigma_R.^2 + (1 ./ max(kappa, 1e-3)) * (1/T)^2);

% WLS Liouville projection
mu0 = liouvilleWLS(mu0_raw, sigmas, ones(n,1), trA0);

end

% -------------------------------------------------------------------------

function kappa = updateKappa(R_bar)
% updateKappa  Solve I₁(κ)/I₀(κ) = R_bar for concentration κ.
%   Uses Mardia & Jupp (2000) p.86 asymptotic for R_bar > 0.85
%   to avoid besseli() overflow and Newton-Raphson divergence.

if R_bar < 1e-6
    kappa = 0;
    return
end
if R_bar > 0.9999
    kappa = 1e6;
    return
end

% Asymptotic approximation for dense clusters (κ → ∞)
if R_bar > 0.85
    kappa = R_bar * (2 - R_bar^2) / (1 - R_bar^2);
    return
end

% Newton-Raphson: solve A(κ) = I₁(κ)/I₀(κ) = R_bar
kappa = 1.0;
for it = 1:20
    A  = besseli(1, kappa) / besseli(0, kappa);
    dA = 1 - A^2 - A / kappa;
    if abs(dA) < 1e-12, break; end
    kappa = max(kappa - (A - R_bar) / dA, 1e-6);
    if abs(A - R_bar) < 1e-10, break; end
end

end

%% =========================================================================
%% MODULE 8 — PARAMETRIC COMB FIT (Rank S)
%% =========================================================================

function [mu0, cost_total, sigmas, inlier_mask] = ...
    peigFitSingle(M_h, k_hat, weights, n, omega, T, trA0, mu0_init, delta)
% peigFitSingle  One run of the k̂-aware parametric comb fit (spec §8).
%   Folds each ν back via k̂: nu_res = Re(ν) + j*wrap(Im(ν) - k̂·ω).
%   This is the structural advantage over generic methods (spec §4.4).
%   WLS Liouville projection applied post-convergence.

MAX_ITER  = 100;
TOL_COST  = 1e-8;
N         = numel(M_h);
mu0       = mu0_init(:);
cost_prev = Inf;

% Folded residuals (fixed, used at every iteration)
nu_res_all = real(M_h) + 1j .* wrapToHalf(imag(M_h) - k_hat .* omega, omega);

assignment  = zeros(N, 1);
inlier_mask = true(N, 1);

for iter = 1:MAX_ITER
    %% -- Assignment --
    % N×n cylindrical distance matrix (vectorized)
    D_ni = zeros(N, n);
    for i = 1:n
        D_ni(:, i) = dCyl(nu_res_all, mu0(i), omega);
    end

    [best_dist, best_i] = min(D_ni, [], 2);
    inlier_mask = best_dist <= delta;
    assignment  = best_i;

    cost_total = sum(weights(inlier_mask) .* best_dist(inlier_mask).^2);

    %% -- Update centroids --
    for i = 1:n
        in_i = find(inlier_mask & assignment == i);
        if isempty(in_i)
            % Fallback: nearest point in folded spectrum
            [~, fb_idx] = min(D_ni(:, i));
            mu0(i) = nu_res_all(fb_idx);
        else
            res_i  = nu_res_all(in_i);
            w_i    = weights(in_i);
            mu0(i) = sum(w_i .* real(res_i)) / sum(w_i) + ...
                     1j * circularMean(imag(res_i), w_i, T);
        end
    end

    if abs(cost_total - cost_prev) < TOL_COST * (1 + abs(cost_prev)), break; end
    cost_prev = cost_total;
end

%% --- Final intra-cluster dispersions ---
sigmas = zeros(n, 1);
for i = 1:n
    in_i = find(inlier_mask & assignment == i);
    if isempty(in_i)
        sigmas(i) = 1e-8;
    else
        res_i   = nu_res_all(in_i);
        w_i     = weights(in_i);
        d_i     = dCyl(res_i, mu0(i), omega);
        sigmas(i) = max(sqrt(sum(w_i .* d_i.^2) / sum(w_i)), 1e-8);
    end
end

%% --- WLS Liouville projection ---
mu0 = liouvilleWLS(mu0, sigmas, ones(n,1), trA0);

end

% -------------------------------------------------------------------------

function [mu0_best, quality] = peigFitMultistart(M_h, k_hat, weights, n, omega, T, ...
    trA0, centroids, sigmas_seed, multiplicities, opts)
% peigFitMultistart  Multi-start comb fit (GUD-004): M runs from perturbed seeds.
%   Run 1: exact seed. Runs 2..n_multistart: perturbations bounded by sigma_i.
%   Retains the run minimizing total weighted cost.

n_runs    = opts.n_multistart;
delta     = max(3 * mean(sigmas_seed), 1e-8);
mu0_seed  = centroids(:);
best_cost = Inf;
mu0_best  = mu0_seed;
best_sigma = sigmas_seed;
best_inlier = true(numel(M_h), 1);

for run = 1:n_runs
    if run == 1
        mu0_init = mu0_seed;
    else
        mu0_init = mu0_seed + ...
            complex(randn(n,1) .* sigmas_seed, randn(n,1) .* sigmas_seed / T);
    end

    [mu0_run, cost_run, sigmas_run, inlier_run] = peigFitSingle(M_h, k_hat, weights, ...
        n, omega, T, trA0, mu0_init, delta);

    if cost_run < best_cost
        best_cost   = cost_run;
        mu0_best    = mu0_run;
        best_sigma  = sigmas_run;
        best_inlier = inlier_run;
    end
end

quality.sigma         = best_sigma;
quality.liouville_res = abs(sum(mu0_best) - trA0);
quality.outlier_frac  = sum(~best_inlier) / numel(M_h);
quality.multiplicity  = multiplicities;

end

%% =========================================================================
%% OUTER LOOP UTILITIES (from findTrueFloquet)
%% =========================================================================

function res = computeResidual(mu0, trA0, useAbsCriterion)
% computeResidual  Floquet-Liouville trace residual (relative or absolute).
traceErr = abs(sum(mu0) - trA0);
if useAbsCriterion
    res = traceErr;
else
    res = traceErr / abs(trA0);
end
end

% -------------------------------------------------------------------------

function msg = buildStatusMsg(status, h, hmax, resLast, thr, nIter, stagWin, hBest, resBest)
% buildStatusMsg  One-line exit summary.
switch status
    case 0
        msg = sprintf('Converged  h=%d | traceRes=%.4e (< %.2e)', h, resLast, thr);
    case 1
        msg = sprintf( ...
            'Stagnated at h=%d (%d iters). Best: h=%d, traceRes=%.4e. Try larger hmax.', ...
            h, nIter, hBest, resBest);
        if nIter < stagWin
            msg = [msg, sprintf(' [window=%d not yet full]', stagWin)];
        end
    case 2
        msg = sprintf( ...
            'hmax reached (h=%d). Best: h=%d, traceRes=%.4e >= %.2e — increase hmax.', ...
            hmax, hBest, resBest, thr);
    otherwise
        msg = sprintf('Unknown status %d at h=%d.', status, h);
end
end
