function [trueMu, trueL, stats] = findTruelmV3(ev, T, nx, varg)
% FINDTRUELMV3  Extract true Floquet exponents from a Hill/HSS spectrum.
%   Single coherent solver. Fixes the two known failure modes of the WIP:
%     (1) circular geometry: ALL distances use the cylindrical (toroidal-imag)
%         metric -- no flat Euclidean DBSCAN/kmeans that tear the Brillouin edge.
%     (2) coincident multiplicity (repeated Floquet exponents): detected by
%         cluster MASS and emitted by REPLICATION, never by k-means splitting
%         a coincident cloud (which returns one copy + parasites).
%
%   ev  : (:,1) full HSS spectrum, (2h+1)*nx values.
%   T   : fundamental period (omega = 2*pi/T).
%   nx  : descriptor dimension = total Floquet multiplicity (sum of mults).
%
%   stats.liouvilleRes uses Sum(mu) = tr(A0) as the ultimate failure detector.

arguments
    ev   (:,1) double
    T    (1,1) double
    nx   (1,1) {mustBeInteger, mustBePositive}
    varg.trA0    = []
    varg.verbose (1,1) logical = false
end

omega   = 2*pi / T;
N       = numel(ev);
rho_ref = N / nx;                 % points per simple mode = (2h+1)
reEV    = real(ev);
imF     = wrapToHalf(imag(ev), omega);

% ---- scales for the cylindrical metric -------------------------------------
% imag is circular over omega; real is the damping axis. Scale real by the
% robust spread of the DENSE core only (parasites must not set the scale).
% ponytail: density via O(N^2) gaussian sum; N is (2h+1)*nx, a few hundred. Fine.
bw = max(1e-3, omega/200);        % sharp kernel to resolve comb cores
dens = zeros(N,1);
for i = 1:N
    dx = reEV(i) - reEV;
    dy = wrapToHalf(imF(i) - imF, omega);
    dens(i) = sum(exp(-(dx.^2 + dy.^2) / (2*bw^2)));
end
coreReal = reEV(dens > median(dens));
s_r = max(iqr(coreReal), 10*bw);  % robust real-axis scale, never zero
s_i = omega;                      % imag scale = one Brillouin period

cdist2 = @(ar,ai,br,bi) ((ar-br)/s_r).^2 + (wrapToHalf(ai-bi,omega)/s_i).^2;

% ---- greedy distinct-peak extraction (circular non-max suppression) --------
R_merge = 0.35;                   % scaled units: < inter-mode gap, > comb spread
tau     = 0.15;                   % a real comb core exceeds tau*max(density);
                                  % isolated parasites fall below and are not peaks
densMax = max(dens);
pool    = true(N,1);
peaksR  = []; peaksI = []; peakDens = [];
while any(pool)
    idxPool = find(pool);
    [dp, j] = max(dens(idxPool));
    if dp < tau * densMax, break; end   % remaining pool = outliers/parasites
    p = idxPool(j);
    peaksR(end+1,1) = reEV(p);  peaksI(end+1,1) = imF(p);  peakDens(end+1,1) = dens(p); %#ok<AGROW>
    near = pool & (cdist2(reEV, imF, reEV(p), imF(p)) <= R_merge^2);
    pool(near) = false;
end

% ---- Voronoi assignment of ALL points to nearest peak (circular) -----------
nP  = numel(peaksR);
D2  = zeros(N, nP);
for k = 1:nP
    D2(:,k) = cdist2(reEV, imF, peaksR(k), peaksI(k));
end
[~, assign] = min(D2, [], 2);
mass = accumarray(assign, 1, [nP 1]);

% ---- multiplicity by mass, keep the nx-budget densest modes ----------------
[~, order] = sort(peakDens, 'descend');        % trust dense peaks first
mult   = zeros(nP,1);
budget = nx;
for k = order(:)'
    if budget <= 0, break; end
    m = max(1, round(mass(k)/rho_ref));
    m = min(m, budget);
    mult(k) = m;
    budget  = budget - m;
end
% if rounding left budget unfilled, top up the densest kept modes
ki = 1;
keptOrder = order(mult(order) > 0);
while budget > 0 && ~isempty(keptOrder)
    kk = keptOrder(mod(ki-1, numel(keptOrder))+1);
    mult(kk) = mult(kk) + 1; budget = budget - 1; ki = ki + 1;
end

% ---- refine each kept peak to its cluster's density-weighted core ----------
keep = find(mult > 0);
trueMu = zeros(nx,1);
conf = zeros(nx,1); rmseV = zeros(nx,1); csz = zeros(nx,1);
w = 0;
for k = keep(:)'
    m   = assign == k;
    evk = ev(m); dk = dens(m);
    [~, c] = max(dk);                          % density peak = ground-truth core
    muk = reEV(find(m,1)); %#ok<NASGU>
    evkR = real(evk); evkI = wrapToHalf(imag(evk), omega);
    muVal = evkR(c) + 1i*evkI(c);
    for r = 1:mult(k)
        w = w + 1;
        trueMu(w) = muVal;
        conf(w)   = mass(k)/rho_ref/mult(k);
        rmseV(w)  = std(evkR);
        csz(w)    = mass(k);
    end
end

trueL = exp(T * trueMu);

% ---- Liouville validation --------------------------------------------------
if ~isempty(varg.trA0)
    if abs(varg.trA0) > 1e-12*nx
        liou = abs(sum(trueMu) - varg.trA0) / abs(varg.trA0);
    else
        liou = abs(sum(trueMu) - varg.trA0);
    end
else
    liou = NaN;
end

stats = struct('method','v3-cylindrical-replicate', 'confidence',conf, ...
    'rmse',rmseV, 'cluster_size',csz, 'outlier_ratio', NaN, ...
    'liouvilleRes', liou, 'nPeaks', nP, 'mult', {mult(keep)});

if varg.verbose
    fprintf('[v3] %d distinct peaks, masses=%s, mult=%s, Liouville=%.2e\n', ...
        nP, mat2str(mass'), mat2str(mult(keep)'), liou);
end
end

function y = wrapToHalf(x, omega)
y = mod(x + omega/2, omega) - omega/2;
end
