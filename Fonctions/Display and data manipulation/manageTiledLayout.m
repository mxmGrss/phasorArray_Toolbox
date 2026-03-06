function T = manageTiledLayout(parent, nx, ny, Tag, varg)
    arguments
        parent = []          
        nx (1,1) double = 1
        ny (1,1) double = 1
        Tag string = "plotTimePhasor"
        varg.ishold = []     
    end

    % 1. Trouver l'ancre
    if isempty(parent)
        targetAx = gca;
    elseif isa(parent, 'matlab.graphics.axis.Axes')
        targetAx = parent;
    else
        targetAx = get(parent, 'CurrentAxes');
        if isempty(targetAx), targetAx = axes('Parent', parent); end
    end

    holdState = varg.ishold;
    if isempty(holdState), holdState = ishold(targetAx); end

    % 2. RECHERCHE DU LAYOUT EXISTANT
    T_existing = [];
    curr = targetAx;
    while ~isempty(curr) && ~isa(curr, 'matlab.ui.Figure')
        if isa(curr, 'matlab.graphics.layout.TiledChartLayout') && strcmp(curr.Tag, Tag)
            T_existing = curr;
            break;
        end
        curr = curr.Parent;
    end

    if isempty(T_existing) && isa(targetAx.Parent, 'matlab.graphics.layout.TiledChartLayout')
        container = targetAx.Parent;
        % Sécurité : on vérifie si targetAx a une propriété Tile
        if isprop(targetAx.Layout, 'Tile')
            tileIdx = targetAx.Layout.Tile;
            siblings = findobj(container.Children, 'flat', 'Tag', Tag, 'Type', 'tiledchartlayout');
            for i = 1:numel(siblings)
                if isequal(siblings(i).Layout.Tile, tileIdx)
                    T_existing = siblings(i); break;
                end
            end
        end
    end

    % 3. LOGIQUE DE DÉCISION
    if holdState && ~isempty(T_existing)
        if all(T_existing.GridSize == [nx, ny])
            T = T_existing;
            return; 
        else
            warning('Dimensions mismatch. Replacing layout.');
        end
    end

    % 4. NETTOYAGE ET RÉCUPÉRATION DE L'INDEX DE TUILE
    if ~isempty(T_existing)
        container = T_existing.Parent;
        % RÉCUPÉRATION SÉCURISÉE DE LA TUILE
        tileIdx = [];
        if isprop(T_existing.Layout, 'Tile')
            tileIdx = T_existing.Layout.Tile;
        end
        delete(T_existing);
    else
        container = targetAx.Parent;
        tileIdx = [];
        if isprop(targetAx.Layout, 'Tile')
            tileIdx = targetAx.Layout.Tile;
        end
        if isvalid(targetAx), delete(targetAx); end
    end

    % 5. CRÉATION
    T = tiledlayout(container, nx, ny, ...
        'TileSpacing', 'compact', 'Padding', 'compact', 'Tag', Tag);
    
    % On n'applique Tile que si on a un index valide (cas imbriqué)
    if ~isempty(tileIdx)
        T.Layout.Tile = tileIdx;
    end
end