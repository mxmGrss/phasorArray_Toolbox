function results = run_all_tests(mode)
%RUN_ALL_TESTS  Run the phasorArray test suite.
%
%   run_all_tests()            runs everything, the development suite.
%   run_all_tests("install")   runs the Install-tagged smoke set only.
%
%   The two answer different questions. The full suite is the regression net:
%   it is meant to fail when a change breaks a contract, and it exercises paths
%   a user never touches directly -- fallback kernels, symbolic payloads,
%   solver residuals. The install set answers "does this copy work here": one
%   check per layer, no optional toolbox, a couple of seconds. A red install
%   run means the installation is wrong; a red development run means the code is.
%
%   Tests carry the Install tag inside a
%       methods (Test, TestTags = {'Install'})
%   block. Adding one to the smoke set means moving it into that block.
%
%   See also installToolbox, checkDependencies.

arguments
    mode (1,1) string {mustBeMember(mode, ["all", "install"])} = "all"
end

projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'Fonctions')));
addpath(fullfile(projectRoot, 'tests'));

suite = matlab.unittest.TestSuite.fromFolder(fullfile(projectRoot, 'tests'));
if mode == "install"
    suite = suite.selectIf(matlab.unittest.selectors.HasTag('Install'));
    banner = 'PhasorArray Toolbox -- install check';
else
    banner = 'PhasorArray Toolbox -- full suite';
end

fprintf('\n=== %s: %d tests ===\n', banner, numel(suite));

runner = matlab.unittest.TestRunner.withTextOutput();
runner.addPlugin(matlab.unittest.plugins.DiagnosticsRecordingPlugin());
results = runner.run(suite);

fprintf('\n===================================================\n');
if all([results.Passed])
    fprintf('SUCCESS: %d/%d passed (%.1f s)\n', ...
        numel(results), numel(results), sum([results.Duration]));
    if mode == "install"
        fprintf('The installation is functional. Run run_all_tests() for the full suite.\n');
    end
else
    fprintf('FAILURE: %d of %d failed.\n', sum([results.Failed]), numel(results));
    fprintf('%s\n', results([results.Failed]).Name);
    if mode == "install"
        fprintf(['A failure here points at the installation: run installToolbox, then\n' ...
                 'checkDependencies("verbose", true).\n']);
    end
end
fprintf('===================================================\n');

if nargout == 0
    clear results
end
end
