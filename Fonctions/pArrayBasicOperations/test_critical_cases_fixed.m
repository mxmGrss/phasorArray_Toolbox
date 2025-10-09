%% Test des cas critiques pour shift2pi
% Script rapide pour tester uniquement les cas problématiques

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

% Visualisation des premières secondes
fprintf('\nVisualization des premières secondes:\n');
for i = 1:min(10, length(t))
    if mod(i-1, round(fs/10)) == 0  % Affiche tous les 0.1s
        dtheta_val = dtheta(min(i, length(dtheta)));
        if dtheta_val > 0
            direction = 'forward';
        else
            direction = 'backward';
        end
        fprintf('  t=%.1fs: theta=%.2f rad, omega=%.2f rad/s, %s\n', ...
                t(i), theta(i), omega(i), direction);
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

%% Test documentaire: Cas "cloche de positivité" (limitation connue)
fprintf('\n=== Test documentaire: Limitation "cloche de positivité" ===\n');
fprintf('Cas pathologique: courte_negative → longue_positive → courte_negative\n');

t_bell = (0:1/fs:6)';
% Crée une "cloche" : négative au début et à la fin, positive au milieu
omega_bell = 2*pi*(10*exp(-((t_bell-3)/1.5).^2) - 1);  % Gaussienne centrée - offset
theta_bell = cumtrapz(t_bell, omega_bell);

% Analyse de cette "cloche"
dtheta_bell = diff(theta_bell);
positive_regions = dtheta_bell > 0;

% Trouve les transitions
transitions = diff([false; positive_regions; false]);
starts = find(transitions == 1);
ends = find(transitions == -1) - 1;

fprintf('Analyse de la "cloche de positivité":\n');
if ~isempty(starts) && ~isempty(ends)
    for i = 1:min(length(starts), length(ends))
        region_length = ends(i) - starts(i) + 1;
        region_duration = t_bell(ends(i)) - t_bell(starts(i));
        region_range = theta_bell(ends(i)) - theta_bell(starts(i));
        
        fprintf('  Région monotone %d: indices %d→%d (%d points, %.2fs, %.1f rad)\n', ...
                i, starts(i), ends(i), region_length, region_duration, region_range);
    end
    
    % Trouve la plus longue région (ce que l'algorithme ne fait PAS)
    [~, best_idx] = max(ends - starts);
    best_start = starts(best_idx);
    best_range = theta_bell(ends(best_idx)) - theta_bell(best_start);
    
    fprintf('  → Région optimale théorique: index %d (%.1f rad)\n', best_start, best_range);
end

% Test avec l'algorithme actuel
fprintf('\nTest algorithme actuel sur la "cloche":\n');
try
    U_bell = sin(2*pi*5*t_bell);
    [~, istart_bell, ~] = shift2pi(theta_bell, t_bell, U_bell);
    
    fprintf('  Algorithme sélectionne: index %d\n', istart_bell);
    
    if ~isempty(starts) && istart_bell ~= best_start
        fprintf('  ⚠ Note: Ce n''est pas la région optimale théorique\n');
        fprintf('  ✓ Mais c''est le comportement attendu (algorithme simple et rapide)\n');
        fprintf('  💡 Pour ce cas, l''ingénieur devrait pré-segmenter les données\n');
    else
        fprintf('  ✓ Par chance, l''algorithme a trouvé une bonne région\n');
    end
    
catch ME
    fprintf('  ✗ Algorithme échoue sur ce cas pathologique: %s\n', ME.message);
    fprintf('  💡 Ceci confirme le besoin de pré-traitement pour les cas extrêmes\n');
end

fprintf('\n=== Recommandations d''usage ===\n');
fprintf('✓ Cas simples: Algorithme automatique suffisant\n');
fprintf('⚠ Cas complexes: Pré-analyse et segmentation recommandées\n');
fprintf('💡 Critères de pré-traitement suggérés:\n');
fprintf('  - Visualiser ω(t) et θ(t) avant traitement\n');
fprintf('  - Identifier manuellement les régions d''intérêt\n');
fprintf('  - Segmenter les données si plusieurs régions monotones distinctes\n');
fprintf('  - Utiliser shift2pi sur chaque segment individuellement\n');

fprintf('\n=== Test terminé ===\n');
