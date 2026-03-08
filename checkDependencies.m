function deps = checkDependencies(options)
%CHECKDEPENDENCIES Diagnose required/optional dependencies for PhasorArray Toolbox.
%
%   deps = CHECKDEPENDENCIES() returns a struct with environment and feature
%   availability for the toolbox.
%
%   deps = CHECKDEPENDENCIES("verbose", true) also prints a human-readable
%   summary to the command window.
%
%   OUTPUT STRUCT FIELDS
%       deps.timestamp
%       deps.matlab.version
%       deps.matlab.release
%       deps.matlab.isR2021bOrNewer
%       deps.matlab.isR2022aOrNewer
%       deps.optional.symbolic
%       deps.optional.controlLpvLtv
%       deps.optional.yalmip
%       deps.optional.solvers
%       deps.features
%       deps.status
%
%   EXAMPLE
%       deps = checkDependencies;
%       deps = checkDependencies("verbose", true);
%
%   See also: test_PhasorArray_basic, test_PhasorArray_advanced

arguments
    options.verbose (1,1) logical = true
end

% Timestamp
deps.timestamp = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");

% MATLAB info
deps.matlab.version = version;
deps.matlab.release = version('-release');
deps.matlab.isR2021bOrNewer = ~verLessThan('matlab', '9.11');
deps.matlab.isR2022aOrNewer = ~verLessThan('matlab', '9.12');

% Optional toolboxes / frameworks
deps.optional.symbolic = ~isempty(ver('symbolic'));

hasLpvss = exist('lpvss', 'file') == 2;
hasLtvss = exist('ltvss', 'file') == 2;
deps.optional.controlLpvLtv = hasLpvss && hasLtvss;

deps.optional.yalmip = exist('sdpvar', 'file') == 2;

% Solver discovery (relevant when YALMIP is present)
deps.optional.solvers.mosek = (exist('mosekopt', 'file') == 2) || (exist('mosekopt', 'file') == 3);
deps.optional.solvers.sdpt3 = exist('sqlp', 'file') == 2;
deps.optional.solvers.sedumi = exist('sedumi', 'file') == 2;
deps.optional.solvers.sdpa = exist('sdpam', 'file') == 2;
solverNames = fieldnames(deps.optional.solvers);
solverFlags = struct2array(deps.optional.solvers);
deps.optional.solvers.any = any(solverFlags);
deps.optional.solvers.availableList = solverNames(solverFlags);

% Feature availability matrix
deps.features.coreToolbox = deps.matlab.isR2021bOrNewer;
deps.features.fastTensorprodPath = deps.matlab.isR2022aOrNewer;
deps.features.symbolicPhasorArray = deps.optional.symbolic;
deps.features.phasorSS_toLPV_LTV = deps.optional.controlLpvLtv;
deps.features.lmiWithYalmip = deps.optional.yalmip;
deps.features.lmiSolve = deps.optional.yalmip && deps.optional.solvers.any;

% Global status
if deps.features.coreToolbox
    deps.status.core = "OK";
else
    deps.status.core = "BLOCKED";
end

missingOptional = strings(0,1);
if ~deps.optional.symbolic
    missingOptional(end+1) = "Symbolic Math Toolbox"; %#ok<AGROW>
end
if ~deps.optional.controlLpvLtv
    missingOptional(end+1) = "Control System Toolbox (lpvss/ltvss)"; %#ok<AGROW>
end
if ~deps.optional.yalmip
    missingOptional(end+1) = "YALMIP"; %#ok<AGROW>
elseif ~deps.optional.solvers.any
    missingOptional(end+1) = "SDP solver (MOSEK/SDPT3/SeDuMi/SDPA)"; %#ok<AGROW>
end

deps.status.missingOptional = missingOptional;

% Human-readable report
if options.verbose
    printReport(deps);
end

end

function printReport(deps)
fprintf('\n============================================================\n');
fprintf('  PhasorArray Toolbox - Dependency Report\n');
fprintf('============================================================\n');
fprintf('Timestamp         : %s\n', char(deps.timestamp));
fprintf('MATLAB            : %s (%s)\n', deps.matlab.version, deps.matlab.release);

fprintf('\n[Core]\n');
printFlag('MATLAB R2021b+ (required)', deps.matlab.isR2021bOrNewer);
printFlag('Fast path (tensorprod, R2022a+)', deps.matlab.isR2022aOrNewer);

fprintf('\n[Optional]\n');
printFlag('Symbolic Math Toolbox', deps.optional.symbolic);
printFlag('Control Toolbox lpvss/ltvss', deps.optional.controlLpvLtv);
printFlag('YALMIP', deps.optional.yalmip);
printFlag('SDP solver available', deps.optional.solvers.any);

if deps.optional.solvers.any
    fprintf('  - Solvers detected : %s\n', strjoin(cellstr(deps.optional.solvers.availableList), ', '));
end

fprintf('\n[Feature Availability]\n');
printFlag('Core toolbox usage', deps.features.coreToolbox);
printFlag('Symbolic workflows', deps.features.symbolicPhasorArray);
printFlag('PhasorSS export to lpvss/ltvss', deps.features.phasorSS_toLPV_LTV);
printFlag('LMI modeling (YALMIP)', deps.features.lmiWithYalmip);
printFlag('LMI solving (YALMIP + solver)', deps.features.lmiSolve);

if ~isempty(deps.status.missingOptional)
    fprintf('\n[Info] Missing optional components:\n');
    for k = 1:numel(deps.status.missingOptional)
        fprintf('  - %s\n', deps.status.missingOptional(k));
    end
end

if deps.matlab.isR2021bOrNewer
    fprintf('\nStatus: OK - Core toolbox can run on this setup.\n');
else
    fprintf('\nStatus: BLOCKED - MATLAB R2021b+ is required.\n');
end
fprintf('============================================================\n\n');
end

function printFlag(label, flag)
if flag
    fprintf('  [OK]   %s\n', label);
else
    fprintf('  [MISS] %s\n', label);
end
end
