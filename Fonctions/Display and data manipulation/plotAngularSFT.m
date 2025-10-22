function  plotAngularSFT(phasorStruct,plotTAPRI,optarg)
arguments
    phasorStruct
    plotTAPRI {mustBeNumericOrLogical} = [true true false false false]
    optarg.plotOmega logical = false
    optarg.plotDebut logical = false
    optarg.xAxes {mustBeMember(optarg.xAxes,{'time','phase','revolution'})} = 'time'
    optarg.orientation {mustBeMember(optarg.orientation,{'ver','hor'})} = 'hor'
    optarg.plotIndexes = 1:numel(phasorStruct)
    optarg.Hm2plot = []
    optarg.lang ='fr';
end
switch optarg.lang
    case 'fr'
        titleStr = 'Phaseurs de ';
        legendStr = 'Phaseur ';
    otherwise  
        titleStr = 'Phasors of';
        legendStr = 'Phasor ';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Affichage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%normalisation de plotTAPRI qui peut contenir des scalaire ou des str
ploto=zeros(size(plotTAPRI));
for plotTAPRIi = 1:numel(plotTAPRI)
    scale_i=plotTAPRI(plotTAPRIi);
    if strcmp(scale_i,"linear")  ||  double(scale_i)==1
        ploto(plotTAPRIi)=1;
    elseif strcmp(scale_i,"log") ||  double(scale_i)==2
        ploto(plotTAPRIi)=2;
    else
        ploto(plotTAPRIi)=0;
    end
end

phasorStruct=phasorStruct(optarg.plotIndexes);
if isempty(optarg.Hm2plot)
    optarg.Hm2plot = {phasorStruct.harmonics};
end

if isnumeric(optarg.Hm2plot)
    optarg.Hm2plot = repmat({optarg.Hm2plot},1,numel(phasorStruct));
end

if iscell(optarg.Hm2plot)
    assert(all(cellfun(@isnumeric,optarg.Hm2plot)),'Hm2plot must be a cell array of numeric arrays')
    assert(numel(phasorStruct)==numel(optarg.Hm2plot),'Hm2plot must have the same number of elements as phasorStruct')
    assert(all(cellfun(@(x) all(x==round(x)),optarg.Hm2plot)),'Hm2plot must be an array of integers')
    %check that all elements of Hm2plot are in the harmonics of the corresponding phasorStruct
    assert(all(arrayfun(@(i) all(ismember(optarg.Hm2plot{i},phasorStruct(i).harmonics)),1:numel(phasorStruct))),'All elements of Hm2plot must be in the harmonics of the corresponding phasorStruct')
end

plotTAPRI=ploto;
if sum(plotTAPRI)>0

    if ishold
        holdvar='on';
    else
        holdvar='off';
    end

    nplot_hor=numel(phasorStruct);
    nplot_vert=sum((plotTAPRI>0))+(optarg.plotOmega>0);

    ff = gcf; % Get the current figure handle
    set(ff, 'Visible', 'off'); % Make the current figure invisible
    pause(1e-1) % Pause for a short time to allow the figure to become invisible
    TT = ff.Children; % Get the children of the current figure

    % Function to create a new tiled layout
    createTiledLayout = @(nplot_vert, nplot_hor) tiledlayout(ff,nplot_vert, nplot_hor, 'TileSpacing', 'tight', 'Padding', 'tight');

    if ~ishold % Check if the hold state is off
        clf(ff); % Clear the current figure
        ff.Tag = ''; % Reset the tag of the figure
        T = createTiledLayout(nplot_vert, nplot_hor); % Create a new tiled layout
    else
        changeLineStyles(ff, {'-', '--', '-.', ':'});
        try
            if strcmp(ff.Tag,'vertical') && isa(TT, 'matlab.graphics.layout.TiledChartLayout')
                transposeTiledLayout(TT);
                TT = ff.Children; % Get the children of the current figure
            end
            if isempty(TT) || ~isa(TT, 'matlab.graphics.layout.TiledChartLayout') || any(TT.GridSize - [nplot_vert, nplot_hor])
                %include a detailed warning witch show wich condition is not met
                if isempty(TT)
                    warning('No tiled layout exists in figure %d.', ff.Number);
                elseif ~isa(TT, 'matlab.graphics.layout.TiledChartLayout')
                    warning('Invalid tiled layout exists in figure %d.', ff.Number);
                elseif any(TT.GridSize - [nplot_vert, nplot_hor])
                    warning('Tiled layout in figure %d has dimensions [%d, %d], but [%d, %d] is required.', ff.Number, TT.GridSize(1), TT.GridSize(2), nplot_vert, nplot_hor);
                end
                
                clf; % Clear the current figure if no valid tiled layout exists
                ff.Tag = ''; % Reset the tag of the figure
                T = createTiledLayout(nplot_vert, nplot_hor); % Create a new tiled layout
                %display a warning message including figure number and hold state
                warning('Hold state is off. Creating a new tiled layout in figure %d. Hold state may be lost.', ff.Number);
                
                    
            else
                T = TT; % Use the existing tiled layout
            end
        catch e
            clf; % Clear the current figure if an error occurs
            ff.Tag = ''; % Reset the tag of the figure
            warning(e.message); % Display the error message
            warning('Error occurred while creating tiled layout. Creating a new tiled layout. Hold state may be lost.'); % Display a warning message
            T = createTiledLayout(nplot_vert, nplot_hor); % Create a new tiled layout
        end
    end

    IDX   = phasorStruct(1).IDX;
    time  = phasorStruct(1).time;
    theta = phasorStruct(1).theta;
    omega = phasorStruct(1).omega;

    if optarg.plotDebut
        IDXp=ones(size(IDX))*2;
    else
        IDXp=IDX;
    end
    switch optarg.xAxes
        case 'time'
            xAx = time(IDXp>1);
            xlabelVar = 'time (s)';
        case 'phase'
            xAx = theta(IDXp>1);
            xlabelVar = 'angle (rad)';
        case 'revolution'
            xAx = theta(IDXp>1)/2/pi;
            xlabelVar = 'revolution';
    end

    for ii = 1:numel(phasorStruct)
        clear temp_sig_plot

        Signal_ii=phasorStruct(ii).signal;
        if ~isrow(Signal_ii)
            Signal_ii=Signal_ii.';
        end

        phasor_v=phasorStruct(ii).phasors;

        titleVar = phasorStruct(ii).name;

        temp_sig_plot.time=xAx;

        %get the indexes of optarg.Hm2plot{ii} in phasorStruct(ii).harmonics
        [~,indPh]=ismember(optarg.Hm2plot{ii},phasorStruct(ii).harmonics);

        temp_sig_plot.sig{1}=Signal_ii(IDXp>1);
        temp_sig_plot.sig{2}=abs(phasor_v(indPh,IDXp>1));
        temp_sig_plot.sig{3}=angle(phasor_v(indPh,IDXp>1));
        temp_sig_plot.sig{4}=real(phasor_v(indPh,IDXp>1));
        temp_sig_plot.sig{5}=imag(phasor_v(indPh,IDXp>1));
        temp_sig_plot.title={titleVar;titleVar+" abs phasor (dB)";titleVar+" angle phasor (dB)";titleVar+" real part phasor";titleVar+" imag part phasor"};

        temp_sig_plot.title={titleVar;...
            strcat(titleStr, " ",titleVar ," (dB)");...
             strcat(titleVar, " angle " ,legendStr, "(dB)");...
             strcat(titleVar, " real part " ,legendStr);...
             strcat(titleVar, " imag part " ,legendStr)};

        iter_i_real=0;
        % ax=gobjects(nplot_vert,nplot_hor);
        if logical(optarg.plotOmega)
            % ax(1,ii)=subplot(nplot_vert,nplot_hor,ii);
            nexttile(T,ii)
            switch optarg.plotOmega
                case 1
                    plot(xAx,omega)
                    title('\omega(t) (rad/s)')
                    xlabel('time (sec)')
                    grid on
                case 2
                    plot(xAx,omega/2/pi)
                    title('f(t) (Hz)')
                    grid on
                    xlabel('time (sec)')
            end
            grid minor
        end

        for iter_i=1:numel(plotTAPRI)
            if logical(plotTAPRI(iter_i))

                iter_i_real=iter_i_real+1;
                ligne_i=iter_i_real+(optarg.plotOmega>0);

                nexttile(T,ii+(ligne_i-1)*nplot_hor);

                if iter_i>1;set(gca,'ColorOrderIndex',1);end

                hold(holdvar);
                plot(temp_sig_plot.time,temp_sig_plot.sig{iter_i})

                xlabel('time (sec)')
                title(temp_sig_plot.title{iter_i})

                if plotTAPRI(iter_i)==2; set(gca,'YScale','log'); end
                
                hold off
                grid on
                grid minor

            end
        end
        legend(legendStr+string(optarg.Hm2plot{ii}))
        legend('Location','best')

    end
        sgtitle('')


    aaa = findobj(gcf,'Type','axes','Visible','on');
    try
        linkaxes(aaa,'x')
    catch ME
        disp(ME.message)
    end


    if strcmp(optarg.orientation,'ver')
        transposeTiledLayout(T)
    end
    setBottomXLabel(ff.Children, xlabelVar) % Set the x-label of the bottom axes
    
    pause(1e-2) % Pause for a short time to allow the figure to become visible
    set(ff,'Visible','on')

    dcm = datacursormode(gcf);
    set(dcm, 'UpdateFcn', @myupdatefcn)
            

end


% Define a custom update function that displays the legend name
function output_txt = myupdatefcn(~, event)
    % Extract the data index and the x and y values
    idx = event.DataIndex;
    x = event.Position(1);
    y = event.Position(2);
    % Extract the legend name for the corresponding series
    leg = event.Target.DisplayName;
    % Create a cell array with custom display text
    output_txt = {[leg], ['time:', num2str(x)], ['value:', num2str(y)]};
end

    function transposeTiledLayout(TL)
        % change the tag of the figure to 'vertical' or 'horizontal', if empty to 'vertical'
        if isempty(ff.Tag)
            ff.Tag = 'vertical';
        else
            if strcmp(ff.Tag, 'vertical')
                ff.Tag = 'horizontal';
            else
                ff.Tag = 'vertical';
            end
        end

        % Get the current grid size
        [nplot_vert2, nplot_hor2] = deal(TL.GridSize(1), TL.GridSize(2));

        % Create a new figure and a new tiled layout with transposed dimensions
        f2 = figure('Visible', 'off');
        TL2 = tiledlayout(f2, nplot_hor2, nplot_vert2, 'TileSpacing', TL.TileSpacing, 'Padding', TL.Padding);

        % Loop through each child and copy it to the new tiled layout
        for axe_i = numel(TL.Children):-1:2
            Axe_i = TL.Children(axe_i);
            
            %check if Axe_i is an axes, if not, skip
            if ~isa(Axe_i,'matlab.graphics.axis.Axes')
                continue                
            end

            %if axe_i == 2
                %axe_i = [2 1]; % We know that the legend is the first child and its associated axe is the second
            %end

            % Copy the child to the new tiled layout
            copyobj(TL.Children(axe_i), TL2);
            npos = Axe_i.Layout.Tile;

            [col, row] = ind2sub([ nplot_hor2,nplot_vert2], npos);
            newIndex = sub2ind([nplot_vert2,nplot_hor2], row, col);
            [npos row col newIndex];

            % Set the new position for the copied child
            TL2.Children(1).Layout.Tile = newIndex; % New children are always in position 1
            try
            legend(TL2.Children(1),Axe_i.Legend.String{:});
            legend('Location','best')
            end
        end

        % Clear the original figure and copy the new tiled layout back to it
        clf(ff)
        copyobj(TL2, ff);
        close(f2)

        % Link axes if possible
        try
            linkaxes(findall(ff, 'type', 'axes'), 'x')
        catch
        end


        
    end
    function setBottomXLabel(tiledLayoutHandle, xlabelStr)
        % Get all the axes in the tiled layout
        axesHandles = findall(tiledLayoutHandle, 'type', 'axes');
        % Get the number of rows in the tiled layout
        numTiles = tiledLayoutHandle.GridSize(1)*tiledLayoutHandle.GridSize(2);
        nplot_vert2 = tiledLayoutHandle.GridSize(1);
        nplot_hor2 = tiledLayoutHandle.GridSize(2);
        % Loop through each axis and delete the x-label
        for i = 1:(numTiles)
            pos = axesHandles(i).Layout.Tile;
            [col, row] = ind2sub([ nplot_hor2,nplot_vert2], pos);
            if row == nplot_vert2
                axesHandles(i).XLabel.String = xlabelStr;
                % set(axesHandles(i),'Xticklabel',get(axesHandles(i),'Xtick'))
            else
                axesHandles(i).XLabel.String = '';
                set(axesHandles(i),'Xticklabel',[])
            end
        end
    end
    function changeLineStyles(fig, lineStyles)
        % Ensure '-' is the first element of lineStyles
        if ~strcmp(lineStyles{1}, '-')
            idx = find(strcmp(lineStyles, '-'));
            if isempty(idx)
                lineStyles = [{'-'}, lineStyles];
            else
                % Move '-' to the first position
                lineStyles = [{'-'}, lineStyles(1:idx-1), lineStyles(idx+1:end)];
            end
        end

        % Get all axes in the figure
        axesHandles = findall(fig, 'type', 'axes');

        % Loop through each axes handle
        for i = 1:length(axesHandles)
            % Get all line objects in the current axes
            lineHandles = findall(axesHandles(i), 'type', 'line');

            % Set the LineStyle of each line or delete if index is maximal
            for j = 1:length(lineHandles)
                currentStyle = lineHandles(j).LineStyle;
                % Find the index of the current style in the lineStyles array
                idx = find(strcmp(lineStyles, currentStyle));
                if ~isempty(idx)
                    if idx == length(lineStyles)
                        % Delete the line if the index is maximal
                        delete(lineHandles(j));
                    else
                        % Otherwise, set to the next style
                        lineHandles(j).LineStyle = lineStyles{idx + 1};
                    end
                end
            end
        end
    end

end