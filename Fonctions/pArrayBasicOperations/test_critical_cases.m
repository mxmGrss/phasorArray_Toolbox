%% Test des cas critiques pour shift# Tracé rapide de l'évolution
fprintf('\nVisualization des premières secondes:\n');
for i = 1:min(10, length(t))
    if mod(i-1, round(fs/10)) == 0  % Affiche tous les 0.1s
        % Version corrigée
        dtheta_val = dtheta(min(i, length(dtheta)));
        if dtheta_val > 0
            direction = 'forward';
        else
            direction = 'backward';
        end
        fprintf('  t=%.1fs: theta=%.2f rad, omega=%.2f rad/s, %s\n', ...
                t(i), theta(i), omega(i), direction);
    end
end % Script rapide pour tester uniquement les cas problématiques

clear; clc;
fprintf('=== Test des cas critiques shift2pi ===\n\n');

%% Configuration
fs = 1000;  % Fréquence d'échantillonnage

%% Test Case 6: Vitesse angulaire variable avec inversions initiales
fprintf('Test du cas critique: omega = 2*pi*(5*t + 2*cos(2*pi*0.3*t))\n');

t = (0:1/fs:4)';
omega = 2*pi*(5*t + 2*cos(2*pi*0.3*t));
theta = cumtrapz(t, omega);

% Analyse de la monotonie
dtheta = diff(theta);
non_monotonic_indices = find(dtheta <= 0);
first_positive_run = find(dtheta > 0, 1, 'first');

fprintf('Analyse des données:\n');
fprintf('  - Échantillons non-monotones: %d/%d (%.1f%%)\n', ...
        length(non_monotonic_indices), length(dtheta), ...
        100*length(non_monotonic_indices)/length(dtheta));

if ~isempty(non_monotonic_indices)
    fprintf('  - Dernier mouvement arrière à t=%.3fs (index %d)\n', ...
            t(non_monotonic_indices(end)), non_monotonic_indices(end));
end

if ~isempty(first_positive_run)
    fprintf('  - Première dérivée positive à t=%.3fs (index %d)\n', ...
            t(first_positive_run), first_positive_run);
end

% Tracé rapide de l'évolution
fprintf('\nVisualization des premières secondes:\n');
for i = 1:min(10, length(t))
    if mod(i-1, round(fs/10)) == 0  % Affiche tous les 0.1s
        fprintf('  t=%.1fs: θ=%.2f rad, ω=%.2f rad/s, dθ/dt=%s\n', ...
                t(i), theta(i), omega(i), ...
                char(dtheta(min(i, length(dtheta))) > 0)*'↗' + ...
                char(dtheta(min(i, length(dtheta))) <= 0)*'↘');
    end
end

% Test de la fonction
U = sin(2*pi*15*t);
fprintf('\nTest de shift2pi...\n');

try
    tic;
    [UHat, istart, i_first_rev, thetaHat, timeHat] = shift2pi(theta, t, U);
    elapsed = toc;
    
    fprintf('Résultats:\n');
    fprintf('  ✓ Exécution réussie en %.2f ms\n', elapsed*1000);
    fprintf('  ✓ Région monotone démarre à l''index: %d (t=%.3fs)\n', istart, t(istart));
    fprintf('  ✓ Données conservées: %d/%d (%.1f%%)\n', ...
            length(theta)-istart+1, length(theta), ...
            100*(length(theta)-istart+1)/length(theta));
    
    if ~isempty(i_first_rev)
        fprintf('  ✓ Première révolution à l''index: %d\n', i_first_rev);
    else
        fprintf('  ⚠ Aucune révolution complète détectée\n');
    end
    
    % Vérification de la monotonie de la région sélectionnée
    theta_selected = theta(istart:end);
    dtheta_selected = diff(theta_selected);
    if all(dtheta_selected > 0)
        fprintf('  ✓ Région sélectionnée strictement monotone\n');
    else
        bad_points = sum(dtheta_selected <= 0);
        fprintf('  ✗ Région sélectionnée contient %d points non-monotones!\n', bad_points);
    end
    
    % Plage angulaire finale
    if istart <= length(theta)
        final_range = theta(end) - theta(istart);
        fprintf('  ✓ Plage angulaire (région monotone): %.2f rad (%.1f tours)\n', ...
                final_range, final_range/(2*pi));
    end
    
catch ME
    fprintf('  ✗ Erreur: %s\n', ME.message);
end

fprintf('\n=== Test terminé ===\n');
