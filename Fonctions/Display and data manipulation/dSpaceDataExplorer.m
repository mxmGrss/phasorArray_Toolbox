classdef dSpaceDataExplorer < handle
    properties
        Xdata
        Ydata
        fig
        mainPanel
        listbox
        axesHandle
        angularSpeedDropdown
        phaseDropdown
        sftFactorField
        harmonicsField
        signalTree
        signalNames
        clearedSignalNames
        flattenedTree
        openedFiles = cell.empty(0,3)
        fileDropdown
        refHarmonicField
        maxOrderField
        treeSig
        selectedFile = 1
        LeftGrid
    end

    methods
        function obj = dSpaceDataExplorer(s,optarg)
            %DSPACE DATA EXPLORER Explore and visualize dSpace data
            %   dSpaceDataExplorer(s) takes a struct variable s and displays its fields
            %   recursively in a hierarchical manner, allowing the user to select and
            %   visualize signals based on their names.

            % Extract Xdata and Ydata
            arguments
                s = [];
                optarg.speedSigName = 'motor_speed_res_rads';
                optarg.phaseSigName = 'motor_angle_res_rad';
                optarg.sftFactor = 1;
                optarg.folder = []
            end
            if isempty(s)
                %open file explorer, load the file, and use it as s
                
                [filename, pathname] = uigetfile('*.mat', 'Select a .mat file');
                if isequal(filename,0)
                    error('User selected Cancel');
                else
                    s = load(fullfile(pathname, filename));
                    %extract the first variable in the struct
                    ff = fieldnames(s);
                    s = s.(ff{1});

                    %initialize the openedFiles cell array
                    obj.openedFiles{end+1,1} = {filename};
                    obj.openedFiles{end,2} = {s};
                    obj.openedFiles{end,3} = {};
                    % Initialize storage for opened files (names and data)
                    % openedFiles{1} will be a cell array of names
                    % openedFiles{2} will be a cell array of corresponding data structs
                    % openedFiles{3} will be a cell array of tree data

                end
            else
                %if s is not empty, use it as the data
                %get the variable name of s in base workspace

                inputname(1)
                %initialize the openedFiles cell array
                obj.openedFiles = {{inputname(1)}, {s}, {}};
            end

            obj.Xdata = s.X;
            obj.Ydata = s.Y;


            reFreshData(obj);


            % Create the figure
            obj.fig = uifigure('Name', 'dSpace Data Explorer', 'NumberTitle', 'off', 'Resize', 'on', ...
                'Position', [20 20 1200 600]);
            MainGrid = uigridlayout(obj.fig,[2,2]);
            MainGrid.RowHeight = {'1x','fit'};

            CtrlPanel = uipanel(MainGrid);
            CtrlPanel.Layout.Row = 2;
            CtrlPanel.Layout.Column = [1 2];
            CtrlGridLayout = uigridlayout(CtrlPanel,[1,1]);
            CtrlGridLayout.RowHeight = {'1x'};
            CtrlGridLayout.ColumnWidth = {'1x'};

            LeftPanel = uipanel(MainGrid);
            LeftPanel.Layout.Row = 1;
            LeftPanel.Layout.Column = 1;
            obj.LeftGrid = uigridlayout(LeftPanel,[2,2]);
            obj.LeftGrid.RowHeight = {'fit', '1x'};


            RightPanel = uipanel(MainGrid);
            RightPanel.Layout.Row = 1;
            RightPanel.Layout.Column = 2;
            RightGrid = uigridlayout(RightPanel,[2,1]);
            RightGrid.RowHeight = {'fit','1x'};

            ActionPanel = uipanel(RightGrid);
            ActionPanel.Layout.Row = 1;
            ActionPanel.Layout.Column = 1;

            FileSelPanel = uipanel(obj.LeftGrid);
            FileSelPanel.Layout.Row = 1;
            FileSelPanel.Layout.Column = [1 2];

            FileSelGrid = uigridlayout(FileSelPanel,[1 3]);

            % Label for file selection dropdown
            labelSel = uilabel(FileSelGrid, 'Text', 'Select File:');
            labelSel.Layout.Column = 1;
            labelSel.Layout.Row = 1;


            % Dropdown to select already opened files
            obj.fileDropdown = uidropdown(FileSelGrid, 'Items', obj.openedFiles{1}, ...
                'ValueChangedFcn', @(src, evt) fileDropdownCallback(obj, src, evt));
            obj.fileDropdown.Layout.Column = 2;
            obj.fileDropdown.Layout.Row = 1;

            % Button to open a new file
            buttonSelFile = uibutton(FileSelGrid, 'Text', 'Open File', ...
                'ButtonPushedFcn', @(src, evt) openFileCallback(obj));
            buttonSelFile.Layout.Column = 3;
            buttonSelFile.Layout.Row = 1;


            % Create the listbox to display the signal names
            obj.listbox = uilistbox(obj.LeftGrid, 'Items', obj.flattenedTree(:,1), ...
                'ItemsData', obj.flattenedTree(:,2), ...
                'Multiselect', 'on', ...
                'ValueChangedFcn', @(src, evt) listboxCallback(obj));
            obj.listbox.Layout.Row = 2;
            obj.listbox.Layout.Column = 2;

            %create a tree on the left
            obj.treeSig = uitree(obj.LeftGrid,'FontSize',12,'Multiselect','on');
            obj.treeSig.Layout.Row = 2;
            obj.treeSig.Layout.Column = 1;

            %populate the tree with the signal tree
            populateTreeSig(obj);


            % Create the axes for plotting
            obj.axesHandle = uiaxes( RightGrid, ...
                'Units', 'normalized');
            obj.axesHandle.Layout.Row = 2;
            obj.axesHandle.Layout.Column = 1;



            set(obj.listbox, 'Items', obj.flattenedTree(:,1));
            set(obj.listbox, 'ItemsData', obj.flattenedTree(:,2));

            % Create UI GRID for the action panel
            ActionGrid = uigridlayout(ActionPanel,[3,4]);
            ActionGrid.RowHeight = {'1x','1x','1x'};
            ActionGrid.ColumnWidth = {'1x','1x','1x','1x'};
            ActionGrid.Padding = [10 10 10 10];


            % Create dropdowns for angular speed and phase
            angularSpeedLabel = uilabel(ActionGrid, 'Text', 'Select Angular Speed:');
            angularSpeedLabel.Layout.Row = 1;
            angularSpeedLabel.Layout.Column = 1;

            obj.angularSpeedDropdown = uidropdown(ActionGrid, 'Items', obj.clearedSignalNames(:,1), ...
                'ItemsData', obj.clearedSignalNames(:,2));
            obj.angularSpeedDropdown.Layout.Row = 1;
            obj.angularSpeedDropdown.Layout.Column = 2;

            phaseLabel = uilabel(ActionGrid, 'Text', 'Select Phase:');
            phaseLabel.Layout.Row = 1;
            phaseLabel.Layout.Column = 3;

            obj.phaseDropdown = uidropdown(ActionGrid, 'Items', obj.clearedSignalNames(:,1),'ItemsData', obj.clearedSignalNames(:,2));
            obj.phaseDropdown.Layout.Row = 1;
            obj.phaseDropdown.Layout.Column = 4;

            % Create numeric field for SFT factor
            sftFactorLabel = uilabel(ActionGrid, 'Text', 'SFT Factor:');
            sftFactorLabel.Layout.Row = 2;
            sftFactorLabel.Layout.Column = 1;

            obj.sftFactorField = uieditfield(ActionGrid, 'numeric', 'Value', optarg.sftFactor);
            obj.sftFactorField.Layout.Row = 2;
            obj.sftFactorField.Layout.Column = 2;

            % Create field for harmonics
            harmonicsLabel = uilabel(ActionGrid, 'Text', 'Harmonics:');
            harmonicsLabel.Layout.Row = 2;
            harmonicsLabel.Layout.Column = 3;

            obj.harmonicsField = uieditfield(ActionGrid, 'text', 'Value', '0 1 2 3 4 5');
            obj.harmonicsField.Layout.Row = 2;
            obj.harmonicsField.Layout.Column = 4;

            % Create button to compute SFT
            computeSFTButton = uibutton(CtrlGridLayout, 'Text', 'Compute SFT', ...
                'ButtonPushedFcn', @(src, evt) computeSFTCallback(obj));
            computeSFTButton.Layout.Row = 1;
            computeSFTButton.Layout.Column = 1;


            % Set the default values for angular speed and phase dropdowns
            try
            angularSpeedIdx = find(strcmp(obj.clearedSignalNames, optarg.speedSigName));
            phaseIdx = find(strcmp(obj.clearedSignalNames, optarg.phaseSigName));
            if isempty(angularSpeedIdx);angularSpeedIdx=1;end
            if isempty(phaseIdx);phaseIdx=1;end
            set(obj.angularSpeedDropdown, 'Value', angularSpeedIdx);
            set(obj.phaseDropdown, 'Value', phaseIdx);
            catch
            end

            % Set the default value for SFT factor
            set(obj.sftFactorField, 'Value', optarg.sftFactor);

            % Create text fields for ref harmonic and max order
            refHarmonicLabel = uilabel(ActionGrid, 'Text', 'Ref Harmonic:');
            refHarmonicLabel.Layout.Row = 3;
            refHarmonicLabel.Layout.Column = 1;

            obj.refHarmonicField = uieditfield(ActionGrid, 'numeric', 'Value', 0);
            obj.refHarmonicField.Layout.Row = 3;
            obj.refHarmonicField.Layout.Column = 2;

            maxOrderLabel = uilabel(ActionGrid, 'Text', 'Max Order:');
            maxOrderLabel.Layout.Row = 3;
            maxOrderLabel.Layout.Column = 3;

            obj.maxOrderField = uieditfield(ActionGrid, 'numeric', 'Value', 10);
            obj.maxOrderField.Layout.Row = 3;
            obj.maxOrderField.Layout.Column = 4;

            % Button to compute SFT with 0:sftFactor:maxOrder
            computeTHDButton = uibutton(CtrlGridLayout, 'Text', 'Compute THD', ...
                'ButtonPushedFcn', @(src, evt) newComputeCallback(obj));
            computeTHDButton.Layout.Row = 1;
            computeTHDButton.Layout.Column = 2;

            computeApproxButton = uibutton(CtrlGridLayout, 'Text', 'Compute Approx hmq', ...
                'ButtonPushedFcn', @(src, evt) computerApproxHmq(obj));
            computeApproxButton.Layout.Row = 1;
            computeApproxButton.Layout.Column = 3;


            %button in ctrl panel, to close all figure except the main one
            closeAllButton = uibutton(CtrlGridLayout, 'Text', 'Close All Figures', ...
                'ButtonPushedFcn', @(src, evt) closeAllFigures(obj));
            closeAllButton.Layout.Row = 1;
            closeAllButton.Layout.Column = 5;

            %button to close the main figure
            closeButton = uibutton(CtrlGridLayout, 'Text', 'Close Explorer', ...
                'ButtonPushedFcn', @(src, evt) close(obj.fig));
            closeButton.Layout.Row = 1;
            closeButton.Layout.Column = 6;

            function closeAllFigures(obj)
                % Close all figures except the main one
                figs = findall(0, 'Type', 'figure');
                for i = 1:length(figs)
                    if figs(i) ~= obj.fig
                        close(figs(i));
                    end
                end
            end





            % Nested callback to open a new file
            function openFileCallback(obj)
                [filename, pathname] = uigetfile('*.mat', 'Select a .mat file');
                if ~isequal(filename, 0)
                    newdata = load(fullfile(pathname, filename));
                    fn = fieldnames(newdata);
                    newdata = newdata.(fn{1});
                    % Append to openedFiles
                    obj.openedFiles{1}{end+1} = filename;
                    obj.openedFiles{2}{end+1} = newdata;
                    obj.openedFiles{3}{end+1} = {};

                    set(obj.fileDropdown, 'Items', obj.openedFiles{1});
                    set(obj.fileDropdown, 'ItemsData', 1: numel(obj.openedFiles{1}));

                    set(obj.fileDropdown, 'Value', numel(obj.openedFiles{1}));
                    obj.selectedFile = numel(obj.openedFiles{1});
                    rebuildApp(obj, newdata);
                end
            end

            % Nested callback to rebuild the app when a file is selected
            function fileDropdownCallback(obj, src, ~)
                idx = get(src, 'Value');
                rebuildApp(obj, obj.openedFiles{2}{idx});
                obj.selectedFile = idx;
            end

            % Helper to re-init figure data
            function rebuildApp(obj, newStruct)
                disph('rebuilding...')
                obj.Xdata = newStruct.X;
                obj.Ydata = newStruct.Y;

                %save current listbox selection
                idx        = get(obj.listbox, 'Value');


                % Clear current selection and refresh listbox
                %set(obj.listbox, 'Value', []);
                %set(obj.listbox, 'String', {});

                idxx = cell2mat(idx);
                idxx = idxx(idxx ~= 0);
                if ~isempty(idxx)
                    selName    = obj.signalNames(idxx);
                else
                    selName = [];
                end
                reFreshData(obj);
                updateUI(obj);


                %restore selection if name still exists, else select first item
                if ~isempty(selName) %&& any(strcmp(get(obj.listbox, 'String'), selName))
                    %iterate through selName, and reselect the index in the listbox if it exists
                    selectionArray = {};
                    for i = 1:length(selName)
                        nIdx = find(strcmp(obj.signalNames, selName{i}));
                        if ~isempty(nIdx)
                            selectionArray = [selectionArray, nIdx];
                        end
                    end
                    set(obj.listbox, 'Value', selectionArray);
                else
                    set(obj.listbox, 'ValueIndex', 1);
                end
                listboxCallback(obj);
            end

            function reFreshData(obj)
                % Refresh SignalTree, signalNames, clearedSignalNames, and flattenedTree
                %get the signal names
                obj.signalNames = {obj.Ydata.Name};

                %parse the signal names to create a hierarchical tree
                obj.signalTree = obj.parseName();


                obj.clearedSignalNames = getFruits(obj,'Signals');

                %flatten the tree to display in the listbox
                obj.flattenedTree = obj.flattenTree();
            end

            function updateUI(obj)

                % Display the signal names in the listbox
                set(obj.listbox, 'Items', obj.flattenedTree(:,1));
                set(obj.listbox, 'ItemsData', obj.flattenedTree(:,2));

                % Set the default values for angular speed and phase dropdowns
                angularSpeedIdx = find(strcmp(obj.clearedSignalNames, optarg.speedSigName));
                phaseIdx = find(strcmp(obj.clearedSignalNames, optarg.phaseSigName));
                set(obj.angularSpeedDropdown, 'Value', angularSpeedIdx);
                set(obj.phaseDropdown, 'Value', phaseIdx);

                populateTreeSig(obj);
            end

        end
        function newComputeCallback(obj,opt)
            if nargin <2
                opt =1
            end


            sftFactor = obj.sftFactorField.Value;
            refHarm = (obj.refHarmonicField.Value);
            maxHarm = (obj.maxOrderField.Value);

            if maxHarm <= 10
                warning('The maximum order is set to 10 or less. Consider increasing the maximum order for better results.');
            end


            fitleredFlattenedTreeIndex = obj.flattenedTree([obj.flattenedTree{:,2}] ~= 0,2);

            % Get the selected angular speed and phase
            angularSpeedIdx = fitleredFlattenedTreeIndex{get(obj.angularSpeedDropdown, 'Value')};
            phaseIdx = fitleredFlattenedTreeIndex{get(obj.phaseDropdown, 'Value')};

            % Get the corresponding Xdata for angular speed and phase
            angularSpeedStruct = obj.Ydata(angularSpeedIdx);
            phaseStruct = obj.Ydata(phaseIdx);

            % Get the selected signals from the listbox
            idx = get(obj.listbox, 'Value');
            realIndexes = cell2mat(idx);
            realIndexes = realIndexes(realIndexes ~= 0);

            % Prepare the signals and names for angularsft
            signals = {obj.Ydata(realIndexes).Data};
            names = {obj.Ydata(realIndexes).Name};
            names = cellfun(@clearName,names,'UniformOutput',false);

            % Find the corresponding Xdata struct based on the raster field
            raster = phaseStruct.Raster;

            xStruct = obj.Xdata(strcmp({obj.Xdata.Raster}, raster));
            time = xStruct.Data;

            % Check that all signals and angularSpeed and phase have the same raster
            if ~all(strcmp({obj.Ydata(realIndexes).Raster}, angularSpeedStruct.Raster)) || ~all(strcmp({obj.Ydata(realIndexes).Raster}, phaseStruct.Raster))
                errordlg('All selected signals must have the same raster as the angular speed and phase signals.', 'Error');
                return;
            end


            if opt == 1
                harmonics = 0:sftFactor:(maxHarm*sftFactor);

                % Unwrap the angular signal and apply the factor
                theta = unwrap(phaseStruct.Data);
                omega = angularSpeedStruct.Data;

            else
                harmonics = 0:maxHarm;
                % Unwrap the angular signal and apply the factor
                theta = unwrap(phaseStruct.Data)*sftFactor;
                omega = angularSpeedStruct.Data*sftFactor;

            end

            % Compute the SFT
            [~, ~, ~, ~, phasorStruct] = angularsft(theta, time, omega, signals, harmonics, names, [0 0 0 0 0]);

            if ~isempty(findall(0,'Type','Figure','Name','THD'))
                f = findall(0,'Type','Figure','Name','THD');
                figure(f)
            else
                f = figure("Name","THD");
                figure(f)
            end
            for ii = 1:numel(phasorStruct)

                %square abs phasors field
                phasorStruct(ii).abs2 = abs(phasorStruct(ii).phasors).^2;

                phasorStruct(ii).refHarm = phasorStruct(ii).abs2(refHarm+1,:);
                phasorStruct(ii).otherHarm = sum(phasorStruct(ii).abs2,1) - phasorStruct(ii).refHarm;

                phasorStruct(ii).THD = sqrt((phasorStruct(ii).otherHarm))./sqrt((phasorStruct(ii).refHarm));

                hold on
                plot(phasorStruct(ii).time,phasorStruct(ii).THD,'DisplayName',phasorStruct(ii).name)

            end
            hold off
            legend('show','Interpreter','None')
            xlabel('Time')
            ylabel('THD')
            %switch to semilogy
            set(gca,'YScale','log')

            %title with filename and sft param
            get(obj.fileDropdown, 'Value')
            title(sprintf("THD for %s with sftFactor = %d, refHarm = %d, maxHarm = %d",get(obj.fileDropdown, 'Value'),sftFactor,refHarm,maxHarm),'Interpreter','none')

            % No plotting for now; we'll complete this callback later
            % Example call: angularsft(theta, time, omega, signals, harmonics, names, ...)
        end

        function listboxCallback(obj)
            % Get the selected items
            idx = get(obj.listbox, 'Value');
            if all(cell2mat(idx) == 0)
                return
            end

            idx = cell2mat(idx);
            idx = idx(idx ~= 0);

            selectedItems = obj.signalNames(idx);
            selectedItems = cellfun(@clearName,selectedItems,'UniformOutput',false);
            % Clear the current plot
            cla(obj.axesHandle);

            % Hold the plot for multiple signals
            hold(obj.axesHandle, 'on');

            % Loop through each selected item and plot the corresponding signal
            for i = 1:length(idx)
                realIndex = idx(i);
                % Find the corresponding Ydata struct
                yStruct = obj.Ydata(realIndex);

                % Find the corresponding Xdata struct based on the raster field
                raster = yStruct.Raster;
                xStruct = obj.Xdata(strcmp({obj.Xdata.Raster}, raster));

                % Plot the data
                if ~isempty(xStruct)
                    time = xStruct.Data;
                    signal = yStruct.Data;
                    plot(obj.axesHandle, time, signal, 'DisplayName', selectedItems{i});
                else
                    errordlg('No matching Xdata found for the selected signal.', 'Error');
                end
            end

            % Add legend and labels
            legend(obj.axesHandle, 'show','Interpreter','None');
            xlabel(obj.axesHandle, 'Time');
            ylabel(obj.axesHandle, 'Signal');
            hold(obj.axesHandle, 'off');
        end

        function computeSFTCallback(obj,opt)
            if nargin <2
                opt =1;
            end
            fitleredFlattenedTreeIndex = obj.flattenedTree([obj.flattenedTree{:,2}] ~= 0,2);

            % Get the selected angular speed and phase
            angularSpeedIdx = fitleredFlattenedTreeIndex{get(obj.angularSpeedDropdown, 'Value')};
            phaseIdx = fitleredFlattenedTreeIndex{get(obj.phaseDropdown, 'Value')};

            % Get the corresponding Xdata for angular speed and phase
            angularSpeedStruct = obj.Ydata(angularSpeedIdx);
            phaseStruct = obj.Ydata(phaseIdx);

            % Get the SFT factor
            sftFactor = obj.sftFactorField.Value;

            % Get the harmonics
            harmonicsStr = get(obj.harmonicsField, 'Value');
            harmonics = str2num(harmonicsStr); %#ok<ST2NM>

            % Get the selected signals from the listbox
            idx = get(obj.listbox, 'Value');
            realIndexes = cell2mat(idx);
            selectedYdata = obj.Ydata(realIndexes);

            % Prepare the signals and names for angularsft
            signals = {selectedYdata.Data};
            names = {selectedYdata.Name};


            if opt == 1
                harmonics = harmonics*sftFactor;

                % Unwrap the angular signal and apply the factor
                theta = unwrap(phaseStruct.Data);
                omega = angularSpeedStruct.Data;

            else
                harmonics = harmonics;
                % Unwrap the angular signal and apply the factor
                theta = unwrap(phaseStruct.Data)*sftFactor;
                omega = angularSpeedStruct.Data*sftFactor;

            end

            % Find the corresponding Xdata struct based on the raster field
            raster = phaseStruct.Raster;
            xStruct = obj.Xdata(strcmp({obj.Xdata.Raster}, raster));
            time = xStruct.Data;

            % Check that all signals and angularSpeed and phase have the same raster
            if ~all(strcmp({selectedYdata.Raster}, angularSpeedStruct.Raster)) || ~all(strcmp({selectedYdata.Raster}, phaseStruct.Raster))
                errordlg('All selected signals must have the same raster as the angular speed and phase signals.', 'Error');
                return;
            end

            % Compute the SFT
            [phasor_cell, theta, omega, IDX, phasorStruct] = angularsft(theta, time, omega, signals, harmonics, names, [0 0 0 0 0],'method','angle');

            % Store the results in the UserData property of the figure for later use
            set(obj.fig, 'UserData', phasorStruct);

            % Display the results in a new figure
            fig = figure('Name', 'SFT Results', 'NumberTitle', 'off', 'Resize', 'on');
            plotAngularSFT(phasorStruct, 'orientation', 'hor', 'plotDebut', 1, 'xAxes', 'time');
        end



        function Tree = parseName(obj)
            %use recursiveParse to parse the signal names
            Tree = struct();
            for i = 1:length(obj.signalNames)
                parts = strsplit(obj.signalNames{i}, '.');
                Tree = recursiveParse(Tree,parts,i);
            end

            obj.signalTree = Tree;

            function parTree = recursiveParse(tree,parsedName,ogIndex)
                if numel(parsedName) == 1
                    if ~isfield(tree, 'Signals')
                        tree.Signals = {};
                    end
                    tree.Signals = [tree.Signals, {parsedName{1};ogIndex}];
                    parTree = tree;
                    return
                end
                part = parsedName{1};
                if ~isfield(tree, part)
                    tree.(part) = struct();
                end
                subtree = tree.(part);
                newsubtree = recursiveParse(subtree,parsedName(2:end),ogIndex);
                tree.(part) = newsubtree;
                parTree = tree;
            end
        end
        function fruits = getFruits(obj,fruitFieldName)
            %find all the fields name 'fruitFieldName' in obj.signalTree and store them in a cell array
            obj.signalNames = {};
            fruits = getFruitsHelper(obj.signalTree, fruitFieldName);

            function fruits = getFruitsHelper(tree, fruitFieldName)
                fruits = {};
                fields = fieldnames(tree);
                for i = 1:length(fields)
                    field = fields{i};
                    if strcmp(field, fruitFieldName)
                        fruits = [fruits; tree.(field)'];
                    else
                        subFruits = getFruitsHelper(tree.(field), fruitFieldName);
                        fruits = [fruits; subFruits];
                    end
                end
            end

            obj.signalNames = {obj.Ydata.Name};

        end

        function flattenedTree = flattenTree(obj, tree, prefix)
            % Flatten the hierarchical tree to display in the listbox
            if nargin < 3
                tree = obj.signalTree;
                prefix = '';
            end
            flattenedTree = {};
            fields = fieldnames(tree);
            for i = 1:length(fields)
                field = fields{i};
                if strcmp(field, 'Signals')
                    % Add the signals at this level
                    for j = 1:size(tree.Signals,2)
                        %find index of the signal in Ydata
                        idx = find(strcmp(obj.clearedSignalNames, [tree.Signals{1,j}]));
                        flattenedTree{end+1, 1} = sprintf('[%3i] %s%s', (idx), prefix, tree.Signals{1,j});
                        flattenedTree{end, 2} = tree.Signals{2,j}; % Signal flag
                    end
                else
                    % Add the field and recurse into the subtree
                    flattenedTree{end+1, 1} = ['      ' prefix field];
                    flattenedTree{end, 2} = 0; % Folder flag
                    subTree = obj.flattenTree(tree.(field), [prefix '|    ']);
                    flattenedTree = [flattenedTree; subTree];
                end
            end
        end
        function computerApproxHmq(obj)
            fitleredFlattenedTreeIndex = obj.flattenedTree([obj.flattenedTree{:,2}] ~= 0,2);

            % Get the selected angular speed and phase
            angularSpeedIdx = fitleredFlattenedTreeIndex{get(obj.angularSpeedDropdown, 'Value')};
            phaseIdx = fitleredFlattenedTreeIndex{get(obj.phaseDropdown, 'Value')};

            % Get the corresponding Xdata for angular speed and phase
            angularSpeedStruct = obj.Ydata(angularSpeedIdx);
            phaseStruct = obj.Ydata(phaseIdx);

            % get time
            time = obj.Xdata(strcmp({obj.Xdata.Raster}, phaseStruct.Raster)).Data;

            selIdx = [];

            %prompt the user to select the signal that contains derivative of the omega (dedicated uifigure with dropdown), and create a button to select Ok, another to select "no signal available"
            % Create a new figure
            domFig = uifigure('Name', 'Select Angular Acceleration Signal', 'NumberTitle', 'off', 'Resize', 'on', ...
                'Position', [20 20 400 200]);
            g = uigridlayout(domFig,[2,2]);
            g.RowHeight = {'fit','fit'};
            g.ColumnWidth = {'1x','1x'};

            uilabel(g,'Text','Select Angular Acceleration Signal:', ...
                'Layout', matlab.ui.layout.GridLayoutOptions('Row',1,'Column',1));

            accelDropdown = uidropdown(g, ...
                'Items', obj.clearedSignalNames(:,1), ...
                'ItemsData', obj.clearedSignalNames(:,2), ...
                'Layout', matlab.ui.layout.GridLayoutOptions('Row',1,'Column',2));

            uButtonOk = uibutton(g,'Text','OK',...
                'Layout',  matlab.ui.layout.GridLayoutOptions('Row',2,'Column',1),...
                'ButtonPushedFcn', @(src,evt) getSelectedSignal(src,evt));

            uButtonNo = uibutton(g,'Text','No signal available',...
                'Layout', matlab.ui.layout.GridLayoutOptions('Row',2,'Column',2),...
                'ButtonPushedFcn', @(src,evt) close(domFig));

            %nested callback to close figure and get selected signal
            function  getSelectedSignal(src,evt)
                selIdx = accelDropdown.Value;
                close(domFig);
            end

            %wait for figure to be closed
            waitfor(domFig)

            %get the selected signal
            if isempty(selIdx)
                %get derivative of speed
                dOmega = gradient(angularSpeedStruct.Data,time);

            else
                %get the selected signal
                disp('Using signal %d (index %d) as the angular acceleration signal',obj.signalNames(selIdx),selIdx)
                dOmega = obj.Ydata(selIdx).Data;

            end



            figure()
            tiledlayout("vert")
            nexttile
            plot(time,dOmega)
            title('Derivative of Angular Speed')
            xlabel('Time')
            ylabel('dOmega/dt')

            nexttile
            plot(time,angularSpeedStruct.Data)
            title('Angular Speed')
            xlabel('Time')
            ylabel('Omega')

            nexttile
            plot(time,dOmega./angularSpeedStruct.Data)
            title('Approximate hmq dω/dt /ω')
            xlabel('Time')
            ylabel('hmq')



            % Get the SFT factor
            sftFactor = obj.sftFactorField.Value;

            % Get the harmonics
            harmonicsStr = get(obj.harmonicsField, 'Value');
            harmonics = str2num(harmonicsStr); %#ok<ST2NM>

            % signal are dom, om, and dom/om
            % Prepare the signals and names for angularsft
            signals = {dOmega,angularSpeedStruct.Data, dOmega./angularSpeedStruct.Data};
            names = {'dω/dt','ω','dω/dt /ω'};

            % Unwrap the angular signal and apply the factor
            theta = unwrap(phaseStruct.Data) ;
            omega = angularSpeedStruct.Data ;

            % Compute the SFT
            [phasor_cell, theta, omega, IDX, phasorStruct] = angularsft(theta, time, omega, signals, harmonics*sftFactor, names, [0 0 0 0 0]);

            figure()
            plotAngularSFT(phasorStruct, 'orientation', 'hor', 'plotDebut', 1, 'xAxes', 'time');


            %compute ω(t-T(t)) and T(t)
            [omega_a_tmT,istartmooving,ifirstrev,~,timeHat,~,~,IDX] = shift2pi(theta,time,omega);

            %theta(i)-2pi = theta(IDX(i))
            %T(i) = time(i)-timeHat(i)

            T = time(istartmooving:end)-timeHat;
            tmooving = time(istartmooving:end);
            figure()
            tiledlayout("vert")
            nexttile
            plot(tmooving,T/sftFactor)
            title('pseudo period T(t)')

            %compute the derivative of T
            dT = gradient(T/sftFactor,tmooving);
            nexttile
            hold on
            rap_om = omega(istartmooving:end)./omega_a_tmT;
            plot(tmooving,(1-rap_om)/sftFactor,"DisplayName","1-ω(t)/ω(t-T(t))")
            plot(tmooving,dT,"DisplayName","dT/dt via grad")
            xlabel('Time')
            title('derivative of T(t)')
            legend('show','Interpreter','None')

            nexttile
            plot(tmooving,rap_om)
            title('ω(t)/ω(t-T(t))')
            xlabel('Time')


            %get selected signals
            idx = get(obj.listbox, 'Value');
            realIndexes = cell2mat(idx);

            extractedSignals = obj.Ydata(realIndexes);

            figure()
            mainT = tiledlayout(1,2);
            leftT = tiledlayout(mainT,"vert");
            leftT.Layout.Tile = 1;
            rightT = tiledlayout(mainT,"vert");
            rightT.Layout.Tile = 2;
            %loop through the signals and compute the hmq
            for i = 1:length(realIndexes)
                signal = extractedSignals(i).Data;
                name = clearName(extractedSignals(i).Name);

                %[omega_a_tmT,istartmooving,ifirstrev,~,timeHat,~,~,IDX] = shift2pi(theta,time,omega)
                [shiftedSignal,~,~,~,~,~,~,~] = shift2pi(theta,time,signal);

                neglectedTerm = shiftedSignal.*omega_a_tmT.*(1-rap_om);
                neglectedTerm(isnan(neglectedTerm))=0;


                shiftedTheta = theta(istartmooving:end)-2*pi;
                neglectedTerm = neglectedTerm.*exp(-1i*shiftedTheta*harmonics);


                nexttile(leftT)
                plot(tmooving,(neglectedTerm))
                legend("Harmonic "+harmonics)
                title(sprintf('Neglected term for %s',name),'Interpreter','none')
                xlabel('Time')
                nexttile(rightT)
                plot(tmooving,abs(cumtrapz(tmooving,neglectedTerm)))
                title(sprintf('Integral of neglected term for %s',name),'Interpreter','none')
                xlabel('Time')


            end


        end

        function populateTreeSig(obj)
            try
            dataStruct = obj.signalTree;
            try
            treeStyle = uistyle("FontWeight","bold");
            endStyle = uistyle("FontAngle","italic","FontColor","blue","Icon","signalico.png");
            catch
            end
            
            try
            %look intro openedFiles{3} if it contains uiTree
                if ~isempty(obj.openedFiles{3}{obj.selectedFile})
                    if isa(obj.openedFiles{3}{obj.selectedFile},'matlab.ui.container.Tree')
                        %make the tree visible
                    obj.openedFiles{3}{obj.selectedFile}.Visible = 'on';
                    %iterate through obj.openedFiles{3} and make them invisible
                    for ii = 1:length(obj.openedFiles{3})
                        if ii ~= obj.selectedFile
                            obj.openedFiles{3}{ii}.Visible = 'off';
                        end
                    end
                    obj.treeSig = obj.openedFiles{3}{obj.selectedFile};
                    return
                    
                    else
                    error('not a tree')
                    end
                end
            catch e
                warning(e.message)
                warning('rebuilding tree from scratch')
            end

            
            %create a tree on the left
            obj.treeSig = uitree(obj.LeftGrid,'FontSize',12,'Multiselect','on');
            obj.treeSig.Layout.Row = 2;
            obj.treeSig.Layout.Column = 1;

            % Recursively populate the tree
            buildNodes(obj.treeSig, dataStruct);
            % Store the tree in the openedFiles cell array
            obj.openedFiles{3}{obj.selectedFile} = obj.treeSig;

            %expand the tree
            expand(obj.treeSig,'all');
            catch E
                warning(E.message)
                warning('Failed to build or populate the tree')
            end


            % Nested function to handle recursion
            function buildNodes(parentNode, s)
                % If s is not a struct, just make a leaf node
                if ~isstruct(s)
                    uitreenode(parentNode, 'Text', (s));
                    disp('found neither struct nor cell, creating a leaf with name %s',s)
                    return;
                end

                % Otherwise, iterate through fields
                fnames = fieldnames(s);
                for i = 1:numel(fnames)
                    thisField = fnames{i};

                    if strcmp(thisField,'Signals')

                        % field is a Signal field, create leaves directly
                        subData = s.(thisField);
                        for cIdx = 1:size(subData,2)
                            uitreenode(parentNode, 'Text',['(' mat2str(subData{2,cIdx}) ') | ' (subData{1,cIdx})],'NodeData',subData{2,cIdx});
                            try
                            %set style
                            addStyle(obj.treeSig,endStyle,'node',parentNode.Children(end));
                            catch
                            end
                        end
                        %continue to next field
                        continue
                    end

                    childNode = uitreenode(parentNode, 'Text', thisField);
                    try
                    addStyle(obj.treeSig,treeStyle,'node',parentNode.Children(end));
                    catch
                    end

                    subData = s.(thisField);
                    % If nested struct, recurse
                    if isstruct(subData)
                        buildNodes(childNode, subData);
                    elseif iscell(subData)
                        % If it's a cell array, create leaves
                        disp('haha')
                        for cIdx = 1:size(subData,2)
                            uitreenode(childNode, 'Text', (subData{1,cIdx}),'NodeData',subData{2,cIdx});
                            %set style
                            addStyle(obj.treeSig,endStyle,'node',parentNode.Children(end));
                        end
                    else
                        % Not struct or cell => direct leaf
                        uitreenode(childNode, 'Text', mat2str(subData));
                    end
                end
            end
        end
    end
end

function clrName = clearName(name)
%break on dot and keep last part
parts = strsplit(name,'.');
clrName = parts{end};
end