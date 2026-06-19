function [mu, info] = findTrueFloquetV3(A, nx, varg)
% FINDTRUEFLOQUETV3  h-refinement driver: grow the HSS truncation until the
%   recovered Floquet exponents stop moving. findTruelmV3 is cheap (~ms); the
%   cost is HmqNEig (the eig). So we simply raise h and re-solve until the
%   spectrum is good enough -- robustness at large h is what carries this.
%
%   [mu, info] = findTrueFloquetV3(A, nx)
%   A   : PhasorArray (periodic matrix). nx : number of Floquet exponents.
%
%   Convergence = recovered mu stable between two successive h (max matched
%   change < tol). Liouville residual is reported as a secondary check (but is
%   unreliable for modes exactly on the +/- omega/2 edge, so it is not the gate).

arguments
    A
    nx   (1,1) {mustBeInteger, mustBePositive}
    varg.T     (1,1) double = 2*pi
    varg.h0    (1,1) double = 8
    varg.hmax  (1,1) double = 200
    varg.tol   (1,1) double = 1e-3
    varg.trA0  = []
    varg.verbose (1,1) logical = false
end

omega = 2*pi/varg.T;
prevMu = [];
hist = [];
h = varg.h0;
while h <= varg.hmax
    E  = HmqNEig(A, h); E = E(:);
    [mu, ~, st] = findTruelmV3(E, varg.T, nx, 'trA0', varg.trA0);
    if isempty(prevMu)
        change = inf;
    else
        change = maxMatchedChange(prevMu, mu, omega);   % stability between h
    end
    hist(end+1,:) = [h, change, st.liouvilleRes]; %#ok<AGROW>
    if varg.verbose
        fprintf('  h=%3d | dMu=%.2e | liou=%.2e\n', h, change, st.liouvilleRes);
    end
    if change < varg.tol
        info = struct('h',h,'change',change,'liouvilleRes',st.liouvilleRes, ...
                      'converged',true,'history',hist);
        return
    end
    prevMu = mu;
    h = ceil(h * 1.5);          % geometric growth keeps the eig count small
end
info = struct('h',h,'change',change,'liouvilleRes',st.liouvilleRes, ...
              'converged',false,'history',hist);
warning('findTrueFloquetV3:noConverge','did not converge by hmax=%d (dMu=%.2e)', varg.hmax, change);
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
