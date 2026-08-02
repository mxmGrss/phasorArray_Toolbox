% filepath: SyncSubfoldersToPath.m
% This script removes subfolders of the current directory that no longer
% exist from the MATLAB path, ensures the current folder is in the path,
% then re-adds subfolders with a single addpath(genpath(...)), and finally saves.

rootFolder = pwd;

% 1) Identify folders on the path that start with the current root folder
currentPathFolders = strsplit(path, pathsep);
foldersUnderRoot = currentPathFolders(startsWith(currentPathFolders, rootFolder));

% 2) Remove any folders that no longer exist on disk
totalFolders = numel(foldersUnderRoot);
prevLength = 0; % Keep track of how many characters were last printed
varPro = -00.01;
for i = 1:totalFolders
    % Drop folders that no longer exist, plus the excluded ones.
    if ~isfolder(foldersUnderRoot{i}) || isExcludedFolder(foldersUnderRoot{i}, rootFolder)
        rmpath(foldersUnderRoot{i});
    end

    % Compute progress and build a simple progress bar line
    progress = i / totalFolders * 100;
    if progress > varPro
        numDots = round(progress / 2);
        progressLine = sprintf('Scanning Path, progress: [%s%s] %.2f%%', ...
                               repmat('.', 1, numDots), ...
                               repmat(' ', 1, 50 - numDots), ...
                               progress);

        % Backspace to overwrite the previous line if needed
        if prevLength > 0
            fprintf(repmat('\b', 1, prevLength - 1));
        end

        % Print the current line and record its length
        fprintf(progressLine);
        prevLength = length(progressLine);
        varPro = varPro + 5;
        pause(0.01)
    end
end
fprintf('\n');

% 3) Ensure the current folder itself is on the path
if ~contains(path, rootFolder)
    addpath(rootFolder);
end

% 4) Recursively add all valid subfolders, minus the excluded ones.
allFolders = strsplit(genpath(rootFolder), pathsep);
keep = ~cellfun(@isempty, allFolders) & ...
       ~cellfun(@(f) isExcludedFolder(f, rootFolder), allFolders);
addpath(strjoin(allFolders(keep), pathsep));

% 5) Save the path. savepath cannot write pathdef.m without elevation and only
%    warns, so check its status instead of assuming it persisted.
if savepath ~= 0
    warning('PhasorArray:installToolbox:pathNotSaved', ...
        ['The MATLAB path was updated for this session but could NOT be saved ' ...
         '(pathdef.m is not writable). Re-run installToolbox in each new session, ' ...
         'or run MATLAB as administrator once to persist it.'])
end

fprintf('Path updated successfully.\nJump to the <a href="matlab:doc PhasorArray">Phasor Array Documentation</a> and <a href="matlab:doc PhasorSS">PhasorSS Documentation</a> \nOr see an exemple with <a href="matlab:open GettingStarted.mlx">GettingStarted</a>\nAll the basis are <a href="matlab:open BasicToolBox.mlx">here</a>\nFor Periodic State Space system jump to <a href="matlab:open Periodic_State_space_example.mlx">PeriodicStateSpace</a>\nFor LMIs exemple see <a href="matlab:open Exemple_Toolbox_LMI.mlx">Exemple_Toolbox_LMI</a>\n');
fprintf('To verify the installation, you can <a href="matlab:test_PhasorArray_basic; test_PhasorArray_advanced">Run All Global Tests</a>\n');
fprintf('To check optional features and solvers, run <a href="matlab:checkDependencies">checkDependencies</a>\n');

%% =========================================================================
function tf = isExcludedFolder(folder, rootFolder)
%ISEXCLUDEDFOLDER  True for hidden folders (.git, .claude...) and scratch/.
%   scratch/ is gitignored and machine-local: a file there whose name collides
%   with a real function shadows it, on that machine only.
tf = contains(folder, [filesep '.']) || ...
     strcmp(folder, fullfile(rootFolder, 'scratch')) || ...
     startsWith(folder, [fullfile(rootFolder, 'scratch') filesep]);
end