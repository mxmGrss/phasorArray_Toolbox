function E = eig(o1, nPt, nvp)
    % EIG Compute the continuous spectrum of the PhasorArray point-by-point.
    %
    %   E = EIG(o1) evaluates the PhasorArray at points in time over one period
    %   and computes the eigenvalues of the instantaneous matrix A(t).
    %   The continuous spectrum of A(t) asymptotically bounds the spectrum 
    %   of the Toeplitz block matrix T_tb(A) (Szegő's theorem).
    %
    %   Inputs:
    %     o1  - (PhasorArray) The input PhasorArray.
    %     nPt - (integer, optional) Number of points for the time evaluation.
    %            - Default: `2^nextpow2(2*o1.h + 1)`
    %
    %   Name-Value Options:
    %     'Track' - (logical) If true, tracks eigenvalues across time steps 
    %               to form continuous trajectories.
    %               Default is false.
    %
    %   Outputs:
    %     E   - (N x 1 x nPt) The eigenvalues of A(t) evaluated point-by-point.
    arguments
        o1 PhasorArray
        nPt = []
        nvp.Track (1,1) logical = false
    end
    if isempty(nPt)
        if o1.h == 0
            nPt = 2;
        else
            nPt = 2^nextpow2(2*o1.h + 1);
        end
    end
    
    th = linspace(0, 2*pi, nPt + 1);
    th = th(1:end-1);
    
    A_t = o1.evalp(th);
    N = size(A_t, 1);
    
    % Capability check rather than a release number: the previous gate named
    % R2022a, and being wrong by one release calls a function that is not there.
    if isempty(which('pageeig'))
        E = zeros(N, 1, nPt);
        for k = 1:nPt
            E(:, 1, k) = eig(A_t(:, :, k));
        end
    else
        E = pageeig(A_t);
    end
    
    if nvp.Track && nPt > 1
        % Track eigenvalues using a predictor and optimal linear assignment
        for k = 2:nPt
            E_curr = E(:, 1, k);
            
            if k == 2
                E_pred = E(:, 1, k-1);
            else
                % Linear predictor (momentum) to naturally resolve crossings
                E_pred = 2 * E(:, 1, k-1) - E(:, 1, k-2);
            end
            
            % Distance matrix between predicted and current eigenvalues.
            % Using values instead of eigenvectors eliminates numeric blowups
            % at exceptional points (bifurcations).
            D = abs(E_pred - E_curr.');
            
            % Global optimal assignment minimizes the total sum of distances,
            % preventing greedy "stealing" that causes massive jumps.
            p = matchpairs(D, 1e6);
            p = sortrows(p, 1);
            idx = p(:, 2);
            
            E(:, 1, k) = E_curr(idx);
        end
    end
end
