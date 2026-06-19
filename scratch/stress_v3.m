function stress_v3()
% STRESS_V3  Synthetic-spectrum breaker hunt for findTruelmV3.
%   Builds Hill/HSS combs with controlled omega, multiplicities and truncation
%   parasites, then probes the geometric edge cases the toolbox tests can't
%   reach (scale, near-collision, Brillouin edge, real-only separation).

rng(0);                                  % reproducible parasites
pass = 0; fail = 0;

% --- baseline (matches the real A) ---------------------------------------
report('baseline omega=1', [2.6465; -2.0164-0.1855i; -2.0164+0.1855i], [1 1 1], 1, 15);

% --- 1. SCALE INVARIANCE: same modes, different omega --------------------
for om = [0.1, 2, 10, 2*pi]
    mu = [2.6465; -2.0164-0.1855i; -2.0164+0.1855i];   % Re fixed, Im will be wrapped by om
    report(sprintf('scale omega=%.3g', om), mu, [1 1 1], om, 15);
end

% --- 2. NEAR-COLLISION: two distinct modes, shrinking frac gap ----------
for beta = [0.30 0.20 0.12 0.08 0.05 0.02]
    mu = [-1.0+1i*beta; -1.0-1i*beta; 2.0];            % conj-like pair, gap=2*beta
    report(sprintf('near-collision gap=%.2f', 2*beta), mu, [1 1 1], 1, 25);
end

% --- 3. BRILLOUIN-EDGE-CENTERED mode (frac = omega/2) -------------------
om = 1;
report('edge-centered frac=om/2', [-0.5+1i*(om/2); 1.5+0i; 0.3+1i*0.1], [1 1 1], om, 20);

% --- 4. CONJ PAIR STRADDLING THE EDGE (frac just inside +/- om/2) -------
report('pair straddling edge', [-1+1i*(om/2-0.02); -1-1i*(om/2-0.02); 2.0], [1 1 1], om, 25);

% --- 5. MODES DIFFERING ONLY IN Re (same frac=0) -----------------------
for dRe = [2.0 0.5 0.2 0.1]
    report(sprintf('real-only sep dRe=%.1f', dRe), [-1+0i; -1-dRe; 3.0], [1 1 1], 1, 20);
end

% --- 6. DEGENERATE at the edge (mult 2 on the Brillouin edge) ----------
report('degenerate-2 on edge', [-0.5+1i*(om/2); 1.5+0i], [2 1], om, 20);

% --- 7. WRONG nx (graceful?) -------------------------------------------
report_raw('wrong nx: ask 4 of 3', [2.6465; -2.0164-0.1855i; -2.0164+0.1855i], [1 1 1], 1, 15, 4);
report_raw('wrong nx: ask 2 of 3', [2.6465; -2.0164-0.1855i; -2.0164+0.1855i], [1 1 1], 1, 15, 2);

fprintf('\n==== %d pass / %d fail ====\n', pass, fail);

    function report(name, mu, mult, om, h)
        report_raw(name, mu, mult, om, h, sum(mult));
    end

    function report_raw(name, mu, mult, om, h, nx_ask)
        ev = synth(mu, mult, om, h);
        T  = 2*pi/om;
        trA0 = sum(mu(:).*mult(:));        % Liouville reference for the synth
        ok = true; info = '';
        try
            [rec,~,st] = findTruelmV3(ev, T, nx_ask, 'trA0', trA0);
            % build the expanded truth list (each mode repeated its multiplicity)
            truth = [];
            for i=1:numel(mu), truth = [truth; repmat(mu(i),mult(i),1)]; end %#ok<AGROW>
            % greedy UNIQUE matching: each recovered claims its nearest unclaimed truth
            recS = rec(:); claimed = false(numel(truth),1); nMatch = 0;
            for q = 1:numel(recS)
                dr = abs(real(truth)-real(recS(q)));
                di = abs(mod(imag(truth)-imag(recS(q))+om/2,om)-om/2);
                d  = hypot(dr,di); d(claimed) = inf;
                [dm, im] = min(d);
                if dm < 2e-2, claimed(im)=true; nMatch=nMatch+1; end
            end
            if nx_ask == sum(mult)
                ok = (nMatch == numel(truth));
            else
                ok = (numel(rec)==nx_ask);   % wrong-nx: just no crash, right count
            end
            info = sprintf('nPk=%d match=%d/%d liou=%.1e rec=%s', st.nPeaks, ...
                           nMatch, numel(truth), st.liouvilleRes, mat2str(round(sort(rec),3).'));
        catch ME
            ok = false; info = ['ERROR ' regexprep(ME.message,'\s+',' ')];
        end
        if ok, pass=pass+1; tag='PASS'; else, fail=fail+1; tag='**FAIL**'; end
        fprintf('%-28s | %-8s | %s\n', name, tag, info);
    end
end

function ev = synth(mu, mult, omega, h)
% Build the HSS comb: each mode mu_k repeated mult_k times, aliases m=-h..h.
% Edge aliases (|m| >= h-1) get truncation drift to mimic real Hill pollution.
ev = [];
drift = 0.15;
for k = 1:numel(mu)
    for c = 1:mult(k)
        for m = -h:h
            v = mu(k) + 1i*m*omega;
            if abs(m) >= h-1
                v = v + drift*(randn + 1i*randn*omega);  % comet-tail parasite
            end
            ev(end+1,1) = v; %#ok<AGROW>
        end
    end
end
end
