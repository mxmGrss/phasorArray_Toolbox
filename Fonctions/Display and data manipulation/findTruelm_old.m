function [trueMu, trueL, stats] = findTruelm(ev, T, nx, varg)
% FINDTRUELM  Extract true Floquet exponents from a Hill/HSS eigenvalue spectrum.
%
%   [trueMu, trueL, stats] = findTruelm(ev, T, nx) identifies the nx true
%   Floquet exponents from the (2h+1)*nx eigenvalues of the Hill matrix.
%   Each true exponent mu_k generates a "comb" of aliases: mu_k + j*m*omega
%   for m in {-h,...,h}. The solvers exploit this structure to cluster the
%   spectrum and recover one representative per comb.
%
%   SYNTAX:
%       [trueMu, trueL, stats] = findTruelm(ev, T, nx)
%       [trueMu, trueL, stats] = findTruelm(ev, T, nx, 'method', 'gap', ...)
%
%   INPUTS:
%       ev   - (:,1) double  Full Hill/HSS eigenvalue vector.
%       T    - (1,1) double  Fundamental period (s or rad depending on convention).
%       nx   - (1,1) int     Number of true modes (size of the original A matrix).
%
%   NAME-VALUE:
%       method   (char)   Clustering solver. Default: 'kde'.
%                  'kde'    Cylindrical KDE + greedy peak picking (original).
%                  'gap'    Gap detection on sorted imaginary parts (no toolbox).
%                  'dbscan' DBSCAN + K-means fallback (Stats Toolbox required).
%                  'kmeans' K-means with outlier trimming (Stats Toolbox required).
%       trA0     (1,1) double  dc harmonic trace tr(A_0) for Liouville validation.
%                              Default: [] (validation skipped, stats.liouvilleRes = NaN).
%       bandwidth (1,1) double KDE bandwidth. Default: 1e-3. ['kde' only]
%       epsilon  (1,1) double  DBSCAN neighbourhood radius. Default: 0 (auto).
%       minPts   (1,1) int     DBSCAN minimum neighbours. Default: 0 (auto).
%       outlierSigma (1,1) double  Intra-cluster outlier rejection threshold
%                                  (multiples of cluster std). Default: 3.
%                                  ['gap', 'kmeans' only]
%
%   OUTPUTS:
%       trueMu - (nx,1) double  True Floquet exponents (wrapped to Brillouin zone).
%       trueL  - (nx,1) double  True Floquet multipliers exp(T * trueMu).
%       stats  - struct         Diagnostics — all fields always present:
%           .method        : Solver used.
%           .confidence    : Normalised peak density (or cluster fill) per mode.
%           .rmse          : Intra-cluster spread of real parts.
%           .cluster_size  : Number of eigenvalues assigned to each mode.
%           .outlier_ratio : Fraction of spectrum marked as outliers.
%           .liouvilleRes  : |sum(trueMu) - trA0| / |trA0|  (NaN if trA0 not given).
%
%   SOLVER GUIDE:
%       'kde'    — safest default; works without toolbox; bandwidth-sensitive.
%       'gap'    — fastest, no parameters, best for clean well-separated combs.
%       'dbscan' — robust to outliers; auto-epsilon; falls back to kmeans if
%                  cluster count != nx.
%       'kmeans' — deterministic with Replicates=5; outliers trimmed post-hoc.
%
%   See also: findTrueFloquet, HmqNEig.

arguments
    ev   (:,1) double
    T    (1,1) double
    nx   (1,1) {mustBeInteger, mustBePositive}
    varg.method      {mustBeMember(varg.method, {'kde','gap','dbscan','kmeans'})} = 'kde'
    varg.trA0                                       = []
    varg.bandwidth   (1,1) double                   = 1e-3
    varg.epsilon     (1,1) double                   = 0
    varg.minPts      (1,1) {mustBeInteger, mustBeNonnegative} = 0
    varg.outlierSigma (1,1) double                  = 3
end

omega     = 2*pi / T;
imag_frac = wrapToHalf(imag(ev), omega);   % map all imaginary parts to [-omega/2, omega/2)

%% --- Dispatch ---

switch varg.method
    case 'kde'
        [trueMu, cs] = solveByCKDE(ev, imag_frac, omega, nx, varg.bandwidth);
    case 'gap'
        [trueMu, cs] = solveByGap(ev, imag_frac, omega, nx, varg.outlierSigma);
    case 'dbscan'
        [trueMu, cs] = solveByDBSCAN(ev, imag_frac, omega, nx, varg.epsilon, varg.minPts);
    case 'kmeans'
        [trueMu, cs] = solveByKMeans(ev, imag_frac, nx, varg.outlierSigma);
end

trueL = exp(T * trueMu);

%% --- Liouville validation (optional) ---

if ~isempty(varg.trA0)
    trA0 = varg.trA0;
    if abs(trA0) > 1e-12 * nx
        liouvilleRes = abs(sum(trueMu) - trA0) / abs(trA0);
    else
        liouvilleRes = abs(sum(trueMu) - trA0);   % absolute fallback: nilpotent dc
    end
else
    liouvilleRes = NaN;
end

%% --- Pack stats ---

stats.method        = varg.method;
stats.confidence    = cs.confidence;
stats.rmse          = cs.rmse;
stats.cluster_size  = cs.cluster_size;
stats.outlier_ratio = cs.outlier_ratio;
stats.liouvilleRes  = liouvilleRes;

end

%% =========================================================================
%% SOLVER 1 — Cylindrical KDE + greedy peak picking  (original algorithm)
%% =========================================================================

function [trueMu, cs] = solveByCKDE(ev, imag_frac, omega, nx, bw)
% solveByCKDE  Cylindrical KDE density over the HSS spectrum, greedy nx-peak extraction.
%   The cylindrical distance wraps the imaginary axis modulo omega so that
%   alias combs appear as co-located density peaks.

N            = length(ev);
pts_per_mode = N / nx;

% --- KDE density at each eigenvalue ---
density = zeros(N, 1);
for i = 1:N
    dx       = real(ev(i)) - real(ev);
    dy_wrap  = wrapToHalf(imag(ev(i)) - imag(ev), omega);
    dist_sq  = dx.^2 + dy_wrap.^2;
    density(i) = sum(exp(-dist_sq / (2 * bw^2)));
end

% --- Greedy peak extraction ---
trueMu    = zeros(nx, 1);
ev_pool   = ev;
dens_pool = density;
frac_pool = imag_frac;

cs.confidence   = zeros(nx, 1);
cs.rmse         = zeros(nx, 1);
cs.cluster_size = zeros(nx, 1);

for k = 1:nx
    if isempty(ev_pool)
        trueMu(k) = trueMu(k-1);
        continue
    end
    [max_dens, max_idx] = max(dens_pool);
    peak_ev  = ev_pool(max_idx);
    trueMu(k) = real(peak_ev) + 1i * frac_pool(max_idx);

    cs.confidence(k) = max_dens / pts_per_mode;

    dx       = real(ev_pool) - real(peak_ev);
    dy_wrap  = wrapToHalf(imag(ev_pool) - imag(peak_ev), omega);
    dist_sq  = dx.^2 + dy_wrap.^2;

    is_inlier  = dist_sq <= (5*bw)^2;
    cs.rmse(k)         = sqrt(mean(dist_sq(is_inlier)));
    cs.cluster_size(k) = sum(is_inlier);

    ev_pool   = ev_pool(~is_inlier);
    dens_pool = dens_pool(~is_inlier);
    frac_pool = frac_pool(~is_inlier);
end

cs.outlier_ratio = length(ev_pool) / N;
end

%% =========================================================================
%% SOLVER 2 — Gap detection on sorted imaginary parts  (no toolbox)
%% =========================================================================

function [trueMu, cs] = solveByGap(ev, imag_frac, omega, nx, outlierSigma)
% solveByGap  Sort imaginary parts on the circle [-omega/2, omega/2), find
%   the nx largest circular gaps to define cluster boundaries, then average
%   real parts within each cluster.  Intra-cluster outliers (real part)
%   are removed before the centroid computation.

N = length(ev);
[frac_sorted, sort_idx] = sort(imag_frac);

% Circular gaps: diff between consecutive points + wrap-around gap
gaps     = [diff(frac_sorted); frac_sorted(1) - frac_sorted(end) + omega];
[~, gap_order] = sort(gaps, 'descend');
boundary_starts = sort(gap_order(1:nx));   % sorted indices where a gap (= cluster start) begins

% Assign cluster labels (1..nx)
labels = zeros(N, 1);
for k = 1:nx
    start_idx = boundary_starts(k) + 1;            % first point after gap k
    if k < nx
        end_idx = boundary_starts(k+1);             % last point before gap k+1
        range   = start_idx : end_idx;
    else
        range   = [start_idx:N, 1:boundary_starts(1)];   % wrap-around cluster
    end
    labels(sort_idx(range)) = k;
end

% Per-cluster centroid with real-part outlier rejection
trueMu          = zeros(nx, 1);
cs.confidence   = zeros(nx, 1);
cs.rmse         = zeros(nx, 1);
cs.cluster_size = zeros(nx, 1);
outlier_count   = 0;

for k = 1:nx
    mask  = labels == k;
    ev_k  = ev(mask);
    fr_k  = imag_frac(mask);

    if isempty(ev_k), continue, end

    % Outlier rejection on real part (truncation noise)
    re_k   = real(ev_k);
    sigma_k = std(re_k) + eps;
    inlier  = abs(re_k - mean(re_k)) <= outlierSigma * sigma_k;

    if ~any(inlier), inlier = true(size(ev_k)); end   % safety: keep all if all rejected

    ev_in = ev_k(inlier);
    fr_in = fr_k(inlier);
    outlier_count           = outlier_count + sum(~inlier);

    trueMu(k)              = mean(real(ev_in)) + 1i * mean(fr_in);
    cs.cluster_size(k)     = sum(inlier);
    cs.rmse(k)             = std(real(ev_in));
    cs.confidence(k)       = cs.cluster_size(k) / (N / nx);
end

cs.outlier_ratio = outlier_count / N;
end

%% =========================================================================
%% SOLVER 3 — DBSCAN + K-means fallback  (Stats Toolbox)
%% =========================================================================

function [trueMu, cs] = solveByDBSCAN(ev, imag_frac, omega, nx, epsilon, minPts)
% solveByDBSCAN  DBSCAN on the 1-D imaginary-fraction axis.  Outliers (label -1)
%   are not forced into any cluster.  If the number of clusters found differs
%   from nx, falls back to K-means on non-outlier points, or on all points.

assert(~isempty(which('dbscan')), ...
    'findTruelm:noToolbox', 'method=''dbscan'' requires the Statistics and Machine Learning Toolbox (R2019a+).');

N = length(ev);

if epsilon <= 0, epsilon = omega / (4 * nx); end
if minPts  <= 0, minPts  = max(2, floor(N / (nx * 4))); end

labels       = dbscan(imag_frac, epsilon, minPts);
valid_labels = unique(labels(labels > 0));
K_found      = numel(valid_labels);

if K_found ~= nx
    warning('findTruelm:dbscanMismatch', ...
        'DBSCAN found %d clusters (expected %d, eps=%.3e, minPts=%d). Falling back to K-means.', ...
        K_found, nx, epsilon, minPts)
    non_out = labels > 0;
    if sum(non_out) >= nx
        km_labels              = kmeans(imag_frac(non_out), nx, 'Replicates', 5);
        labels(non_out)        = km_labels;
        labels(~non_out)       = -1;
    else
        labels       = kmeans(imag_frac, nx, 'Replicates', 5);
    end
    valid_labels = unique(labels(labels > 0));
end

trueMu          = zeros(nx, 1);
cs.confidence   = zeros(nx, 1);
cs.rmse         = zeros(nx, 1);
cs.cluster_size = zeros(nx, 1);

for k = 1:nx
    mask  = labels == valid_labels(k);
    ev_k  = ev(mask);
    fr_k  = imag_frac(mask);
    trueMu(k)          = mean(real(ev_k)) + 1i * mean(fr_k);
    cs.cluster_size(k) = sum(mask);
    cs.rmse(k)         = std(real(ev_k));
    cs.confidence(k)   = cs.cluster_size(k) / (N / nx);
end

cs.outlier_ratio = sum(labels == -1) / N;
end

%% =========================================================================
%% SOLVER 4 — K-means with intra-cluster outlier trimming  (Stats Toolbox)
%% =========================================================================

function [trueMu, cs] = solveByKMeans(ev, imag_frac, nx, outlierSigma)
% solveByKMeans  K-means (K=nx) on imaginary fractions, then real-part outlier
%   rejection within each cluster before centroid computation.

assert(~isempty(which('kmeans')), ...
    'findTruelm:noToolbox', 'method=''kmeans'' requires the Statistics and Machine Learning Toolbox.');

N = length(ev);
labels = kmeans(imag_frac, nx, 'Replicates', 5);

trueMu          = zeros(nx, 1);
cs.confidence   = zeros(nx, 1);
cs.rmse         = zeros(nx, 1);
cs.cluster_size = zeros(nx, 1);
outlier_count   = 0;

for k = 1:nx
    mask  = labels == k;
    ev_k  = ev(mask);
    fr_k  = imag_frac(mask);

    if isempty(ev_k), continue, end

    re_k    = real(ev_k);
    sigma_k = std(re_k) + eps;
    inlier  = abs(re_k - mean(re_k)) <= outlierSigma * sigma_k;

    if ~any(inlier), inlier = true(size(ev_k)); end

    ev_in  = ev_k(inlier);
    fr_in  = fr_k(inlier);
    outlier_count = outlier_count + sum(~inlier);

    trueMu(k)          = mean(real(ev_in)) + 1i * mean(fr_in);
    cs.cluster_size(k) = sum(inlier);
    cs.rmse(k)         = std(real(ev_in));
    cs.confidence(k)   = cs.cluster_size(k) / (N / nx);
end

cs.outlier_ratio = outlier_count / N;
end

%% =========================================================================
%% HELPER
%% =========================================================================

function y = wrapToHalf(x, omega)
% wrapToHalf  Map x to the half-open interval [-omega/2, omega/2).
y = mod(x + omega/2, omega) - omega/2;
end
