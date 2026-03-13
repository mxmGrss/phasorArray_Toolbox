function [trueMu, trueL, stats] = findTruelm(ev, T, nx)
% FINDTRUELM Extraction of Floquet exponents using cylindrical KDE.
%
%   [trueMu, trueL, stats] = findTruelm(ev, T, nx) identifies the nx true
%   Floquet exponents (mu) and multipliers (l) from a pool of eigenvalues 
%   retrieved from a Hill State-Space (HSS) spectrum. It uses a Cylindrical 
%   Kernel Density Estimation (KDE) to cluster points belonging to the same 
%   Floquet mode.
%
%   SYNTAX:
%       [trueMu, trueL, stats] = findTruelm(ev, T, nx)
%
%   INPUTS:
%       ev      - Eigenvalues vector (e.g., from HmqNEig)
%       T       - Fundamental period of the system
%       nx      - Theoretical number of modes (size of the original A matrix)
%
%   OUTPUTS:
%       trueMu  - Identified Floquet exponents (within the Brillouin zone)
%       trueL   - Identified Floquet multipliers (exp(T * mu))
%       stats   - Struct containing diagnostic metrics:
%           - confidence   : Normalized peak density per mode
%           - rmse         : Spread of the identified cluster
%           - cluster_size : Number of points assigned to the mode
%           - outlier_ratio: Percentage of spectrum pollution
%
%   PhasorArray Toolbox
%   Author: Maxime GROSSO
%   Date: March 2026

    arguments
        ev (:,1) double
        T (1,1) double
        nx (1,1) {mustBeInteger, mustBePositive}
    end

    % 1. Fundamental parameters
    Omega = 2*pi / T;
    N = length(ev);
    bw = 1e-3; % Bandwidth for KDE
    
    % Theoretical points per mode (if no truncation/pollution)
    pts_per_mode = N / nx; 
    
    % Initialize statistics structure
    stats.confidence = zeros(nx, 1);
    stats.rmse       = zeros(nx, 1);
    stats.cluster_size = zeros(nx, 1);
    
    % 2. KDE Density Evaluation with Cylindrical Distance
    density = zeros(N, 1);
    for i = 1:N
        dx = real(ev(i)) - real(ev);
        dy = imag(ev(i)) - imag(ev);
        % Wrapping imaginary part to account for [-Omega/2, Omega/2] periodicity
        dy_wrap = mod(dy + Omega/2, Omega) - Omega/2;
        dist_sq = dx.^2 + dy_wrap.^2;
        density(i) = sum(exp(-dist_sq / (2 * bw^2)));
    end
    
    % 3. Peak Extraction (Identifying the nx clusters)
    trueMu = zeros(nx, 1);
    ev_pool = ev;
    dens_pool = density;
    
    for k = 1:nx
        if isempty(ev_pool)
            trueMu(k) = trueMu(k-1);
            continue;
        end
        
        % Find the densest peak
        [max_dens, max_idx] = max(dens_pool);
        peak_ev = ev_pool(max_idx);
        
        % Wrap to Brillouin Zone
        peak_mu = real(peak_ev) + 1i * (mod(imag(peak_ev) + Omega/2, Omega) - Omega/2);
        trueMu(k) = peak_mu;
        
        % Record confidence (Normalized density)
        stats.confidence(k) = max_dens / pts_per_mode;
        
        % Compute distances for this cluster
        dx = real(ev_pool) - real(peak_ev);
        dy = imag(ev_pool) - imag(peak_ev);
        dy_wrap = mod(dy + Omega/2, Omega) - Omega/2;
        dist_sq = dx.^2 + dy_wrap.^2;
        
        % Separate Inliers (cluster members) from Outliers
        is_outlier = dist_sq > (5 * bw)^2;
        is_inlier  = ~is_outlier;
        
        % Cluster dispersion metrics
        cluster_dist_sq = dist_sq(is_inlier);
        stats.cluster_size(k) = length(cluster_dist_sq);
        stats.rmse(k) = sqrt(mean(cluster_dist_sq)); % Ideal: close to 0
        
        % Update pool (keep only outliers for next iteration)
        ev_pool = ev_pool(is_outlier);
        dens_pool = dens_pool(is_outlier);
    end
    
    % 4. Global statistics
    stats.outlier_ratio = length(ev_pool) / N; % Percentage of spectral pollution
    
    % 5. Return multipliers
    trueL = exp(T * trueMu);
end