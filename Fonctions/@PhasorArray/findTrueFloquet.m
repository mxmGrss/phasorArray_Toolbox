function [true_mu, true_l, stats, ev_m, ev_l] = findTrueFloquet(o1, T, nvp)
% FINDTRUEFLOQUET  Recover true Floquet exponents and multipliers (V3 Engine)
%
%   [true_mu, true_l, stats, ev_m, ev_l] = findTrueFloquet(o1, T, Name, Value)
%
%   Replaces the legacy k-means solver with the V3 topological solver.

arguments
    o1
    T    (1,1) double = 2*pi
    nvp.hinit             = []
    nvp.hmax              (1,1) {mustBeInteger, mustBePositive} = 200
    nvp.thresholdResidual (1,1) double = 1e-3
    nvp.trA0              = []
    nvp.verbose           = false
    
    % Ignored legacy arguments for drop-in compatibility
    nvp.stagnationWindow  = []
    nvp.stagnationRatio   = []
    nvp.method            = []
end

if isa(o1, 'PhasorArray')
    nx = o1.size(1);
    if isempty(nvp.hinit), nvp.hinit = o1.h; end
    if isempty(nvp.trA0), nvp.trA0 = trace(o1.phas(0)); end
else
    nx = size(o1, 1);
    if isempty(nvp.hinit), nvp.hinit = 8; end
end

omega = 2*pi/T;
prevMu = [];
hist = [];
h = nvp.hinit;
change = inf;

while h <= nvp.hmax
    ev_m = HmqNEig(o1, h, T); ev_m = ev_m(:);
    ev_l = exp(T * ev_m);
    [true_mu, true_l, st] = findTruelm(ev_m, T, nx, 'trA0', nvp.trA0);
    
    if isempty(prevMu)
        change = inf;
    else
        change = maxMatchedChange(prevMu, true_mu, omega);   % stability between h
    end
    
    hist(end+1,:) = [h, change, st.liouvilleRes]; %#ok<AGROW>
    if nvp.verbose
        fprintf('  h=%3d | dMu=%.2e | liou=%.2e\n', h, change, st.liouvilleRes);
    end
    
    if change < nvp.thresholdResidual
        stats = struct('h',h,'change',change,'liouvilleRes',st.liouvilleRes, ...
                       'converged',true,'history',hist);
        return
    end
    
    prevMu = true_mu;
    h = ceil(h * 1.5);          % geometric growth keeps the eig count small
end

stats = struct('h',h,'change',change,'liouvilleRes',st.liouvilleRes, ...
               'converged',false,'history',hist);
warning('findTrueFloquet:noConverge','did not converge by hmax=%d (dMu=%.2e)', nvp.hmax, change);
end

function c = maxMatchedChange(a, b, omega)
% greedy unique match of b onto a, max circular distance between matched pairs
a = a(:); b = b(:); claimed = false(numel(a),1); c = 0;
for q = 1:numel(b)
    dr = abs(real(a)-real(b(q)));
    di = abs(mod(imag(a)-imag(b(q))+omega/2, omega) - omega/2);
    d  = hypot(dr,di); d(claimed) = inf;
    [dm, im] = min(d);
    claimed(im) = true; c = max(c, dm);
end
end
