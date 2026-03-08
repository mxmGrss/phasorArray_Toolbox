function advice = truncationAdvisor(problemType, options)
%TRUNCATIONADVISOR Recommend practical harmonic truncation settings.
%
%   advice = TRUNCATIONADVISOR(problemType) returns recommended truncation
%   settings and convergence checks for typical PhasorArray workflows.
%
%   advice = TRUNCATIONADVISOR(problemType, "profile", profile) selects a
%   preset profile:
%       - "quick"        : fast exploratory runs
%       - "balanced"     : default trade-off (recommended)
%       - "high-accuracy": conservative high-fidelity settings
%
%   Supported problemType values:
%       - "lyap"
%       - "riccati"
%       - "lmi"
%       - "simulation"
%       - "generic"
%
%   The function does not solve equations; it provides initialization and
%   stopping criteria that can be fed into existing solvers.
%
%   EXAMPLES
%       A = truncationAdvisor("lyap");
%       A = truncationAdvisor("lmi", "profile", "high-accuracy");
%       A = truncationAdvisor("riccati", "h0", 8, "hmax", 120);
%
%   See also: lyap, RicHarmonicKlein, checkDependencies

arguments
    problemType (1,1) string {mustBeMember(problemType,["lyap","riccati","lmi","simulation","generic"])} = "generic"
    options.profile (1,1) string {mustBeMember(options.profile,["quick","balanced","high-accuracy"])} = "balanced"
    options.h0 (1,1) double {mustBeInteger,mustBePositive} = 0
    options.hmax (1,1) double {mustBeInteger,mustBePositive} = 0
    options.targetResidual (1,1) double {mustBePositive} = NaN
    options.verbose (1,1) logical = true
end

% Baseline presets (heuristics)
preset = localPreset(problemType, options.profile);

% User overrides
if options.h0 > 0
    preset.h0 = options.h0;
end
if options.hmax > 0
    preset.hmax = options.hmax;
end
if ~isnan(options.targetResidual)
    preset.targetResidual = options.targetResidual;
end

% Derived LMI truncation recommendation
if problemType == "lmi"
    hLMI = max(2*preset.h0, preset.h0 + 4);
    hP = max(round(0.5*hLMI), 4);
    hA = max(round(0.5*hLMI), 4);
else
    hLMI = NaN;
    hP = NaN;
    hA = NaN;
end

% Convergence checklist
convergence = struct();
convergence.primaryMetric = "resnorm";
convergence.secondaryMetric = "delta-gain";
convergence.targetResidual = preset.targetResidual;
convergence.relativeGainTolerance = preset.gainTolerance;
convergence.maxHarmonicOrder = preset.hmax;
convergence.updateStep = preset.hStep;
convergence.stopRule = sprintf("Stop when residual < %.1e and gain delta < %.1e", ...
    preset.targetResidual, preset.gainTolerance);

% Recommended solver options
solverOptions = struct();
solverOptions.h = preset.h0;
solverOptions.autoUpdateh = true;
solverOptions.thresholdResidual = preset.targetResidual;
solverOptions.hmax = preset.hmax;
solverOptions.hStep = preset.hStep;
solverOptions.maxIter = preset.maxIter;

% Build output
advice = struct();
advice.problemType = problemType;
advice.profile = options.profile;
advice.recommended = struct(...
    "h0", preset.h0, ...
    "hmax", preset.hmax, ...
    "hStep", preset.hStep, ...
    "targetResidual", preset.targetResidual, ...
    "hLMI", hLMI, ...
    "hP", hP, ...
    "hA", hA);
advice.solverOptions = solverOptions;
advice.convergence = convergence;
advice.notes = preset.notes;

if options.verbose
    localPrint(advice);
end

end

function preset = localPreset(problemType, profile)
% Heuristics tuned for practical first runs on medium-size models.

switch profile
    case "quick"
        base = struct("h0",4,"hmax",24,"hStep",2,"targetResidual",1e-3,"gainTolerance",1e-2,"maxIter",30);
    case "balanced"
        base = struct("h0",6,"hmax",60,"hStep",2,"targetResidual",1e-5,"gainTolerance",1e-3,"maxIter",60);
    case "high-accuracy"
        base = struct("h0",10,"hmax",120,"hStep",2,"targetResidual",1e-7,"gainTolerance",1e-4,"maxIter",120);
end

switch problemType
    case "lyap"
        base.notes = [...
            "Use residual norm on Lyapunov equation as primary indicator.", ...
            "Increase h until residual stabilizes and symmetry residual remains low."];
    case "riccati"
        base.h0 = base.h0 + 2;
        base.notes = [...
            "Riccati generally needs higher h than Lyapunov.", ...
            "Track both residual and gain variation between iterations."];
    case "lmi"
        base.h0 = max(base.h0, 8);
        base.targetResidual = min(base.targetResidual, 1e-6);
        base.notes = [...
            "Use separate truncation orders (hLMI, hP, hA).", ...
            "Start moderate, then increase hLMI to check robustness of feasibility."];
    case "simulation"
        base.h0 = max(base.h0, 4);
        base.notes = [...
            "Match h with dominant harmonic content in excitation and plant.", ...
            "Validate in time-domain reconstruction after each increment."];
    otherwise
        base.notes = [...
            "Start balanced and escalate only if residual or gain variation is high.", ...
            "Prefer incremental h updates over large jumps."];
end

preset = base;
end

function localPrint(advice)
fprintf('\n============================================================\n');
fprintf('  Truncation Advisor\n');
fprintf('============================================================\n');
fprintf('Problem type      : %s\n', advice.problemType);
fprintf('Profile           : %s\n', advice.profile);
fprintf('h0 / hmax / hStep : %d / %d / %d\n', ...
    advice.recommended.h0, advice.recommended.hmax, advice.recommended.hStep);
fprintf('Target residual   : %.1e\n', advice.recommended.targetResidual);

if ~isnan(advice.recommended.hLMI)
    fprintf('LMI split orders  : hLMI=%d, hP=%d, hA=%d\n', ...
        advice.recommended.hLMI, advice.recommended.hP, advice.recommended.hA);
end

fprintf('Stop rule         : %s\n', advice.convergence.stopRule);
fprintf('Notes:\n');
for k = 1:numel(advice.notes)
    fprintf('  - %s\n', advice.notes(k));
end
fprintf('============================================================\n\n');
end
