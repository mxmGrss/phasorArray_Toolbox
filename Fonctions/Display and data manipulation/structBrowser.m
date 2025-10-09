function [fields,values] = structBrowser(s)
    %STRUCTBROWSER Explore a struct variable recursively and display its fields
    %   structBrowser(s) takes a struct variable s and displays its fields
    %   recursively in a hierarchical manner.

    % Create the figure
    fig = figure('Name', 'Struct Browser', 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', 'Resize', 'on');

    % Create the main panel
    mainPanel = uipanel('Parent', fig, 'Units', 'normalized', 'Position', [0 0 1 1]);

    % Create the listbox to display the struct fields
    listbox = uicontrol('Parent', mainPanel, 'Style', 'listbox', 'Units', 'normalized', 'Position', [0 0 1 0.9], 'FontName', 'Courier New', 'Callback', @(src, event) listboxCallback(src, s));

    % Get the struct fields recursively
    [fields, values] = getStructFields(s);

    % Create the export button
    exportButton = uicontrol('Parent', mainPanel, 'Style', 'pushbutton', 'String', 'Export', 'Units', 'normalized', 'Position', [0 0.9 1 0.1], 'Callback', @(src, event) exportToText(fields));

    % Display the fields in the listbox
    set(listbox, 'String', fields);
    set(listbox, 'UserData', values); % Store the values in the listbox's UserData property
end

function [fields, values] = getStructFields(s, prefix, isAlreadyMultiple)
    if nargin < 2
        prefix = '';
    end
    if nargin < 3
        isAlreadyMultiple = false;
    end

    fields = {};
    values = {};
    if isstruct(s)
        fn = fieldnames(s);
        for i = 1:length(fn)
            field = fn{i};
            isMultiple = isAlreadyMultiple;
            value = {s.(field)}; % Store all values in a cell array
            fieldType = class(s(1).(field)); % Get the type from the first element
            if isstruct(value{1})
                if numel(value{1}) > 1
                    isMultiple = true;
                    sizeStr = mat2str(size(value{1}));
                    fields{end+1,1} = sprintf('%s%s: %s (%s)', prefix, field, sizeStr, fieldType);
                else
                    fields{end+1,1} = sprintf('%s%s: (%s)', prefix, field, fieldType);
                end
                [subFields, subValues] = getStructFields(value{1}, [prefix '  |  '], isMultiple);
                values{end+1,1} = fieldnames(value{1});
                fields = [fields; subFields];
                values = [values; subValues];
            else
                sizeStr = mat2str(size(value{1}));
                if ~isMultiple
                    if isnumeric(value{1})
                        valueStr = num2str(value{1});
                        if numel(value{1}) > 3
                            valueStr = [num2str(value{1}(1:3)) '...'];
                        end
                    elseif ischar(value{1})
                        valueStr = value{1};
                        if numel(value{1}) > 30
                            valueStr = [value{1}(1:30) '...'];
                        end
                    else
                        valueStr = '';
                    end
                else
                    valueStr = '';
                end
                fields{end+1,1} = sprintf('%s%s: %s (%s): %s', prefix, field, sizeStr, fieldType, valueStr);
                values{end+1,1} = value;
            end
        end
    end
end

function exportToText(fields)
    % Convert the hierarchical fields to a single string with line returns
    textString = strjoin(fields, '\n');

    % Create a new figure to display the text string
    fig = figure('Name', 'Exported Text', 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', 'Resize', 'on');
    uicontrol('Parent', fig, 'Style', 'edit', 'Max', 2, 'Min', 0, 'Units', 'normalized', 'Position', [0 0 1 1], 'FontName', 'Courier New', 'String', textString);
end

function listboxCallback(src, s)
    % Get the selected item
    idx = get(src, 'Value');
    items = get(src, 'String');
    selectedItem = items{idx};

    % Get the corresponding value from the UserData property
    values = get(src, 'UserData');
    value = values{idx};

    value

    % Display or plot the value based on its type
    if ischar(value{1})
        % Display the full string in a new figure
        fig = figure('Name', 'Field Content', 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', 'Resize', 'on');
        uicontrol('Parent', fig, 'Style', 'edit', 'Max', 2, 'Min', 0, 'Units', 'normalized', 'Position', [0 0 1 1], 'FontName', 'Courier New', 'String', value{1});
    elseif isnumeric(value{1}) && isvector(value{1})
        % Plot the numeric vector in a new figure
        fig = figure('Name', 'Field Content', 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', 'Resize', 'on');
        plot(value{1});
        title(selectedItem);
    elseif isstruct(value{1}) && numel(value{1}) > 1
        % Display a listbox with the names of the struct array
        names = arrayfun(@(x) x.Name, value{1}, 'UniformOutput', false);
        fig = figure('Name', 'Field Content', 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', 'Resize', 'on');
        uicontrol('Parent', fig, 'Style', 'listbox', 'Max', 2, 'Min', 0, 'Units', 'normalized', 'Position', [0 0 1 1], 'FontName', 'Courier New', 'String', names);
    else
        % Display the value in a new figure
        fig = figure('Name', 'Field Content', 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', 'Resize', 'on');
        uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0 0 1 1], 'FontName', 'Courier New', 'String', mat2str(value{1}));
    end
end