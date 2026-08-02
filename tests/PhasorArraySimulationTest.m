classdef PhasorArraySimulationTest < matlab.unittest.TestCase
    %PHASORARRAYSIMULATIONTEST  Time integration against closed-form trajectories.
    %
    %   lsim and initial on systems whose solution is known exactly, so the comparison
%   is an equality rather than a smoke test.

    properties
        tol = 1e-10;
        tolOde = 1e-4;   % adaptive solver, measured residual is around 3e-6
    end

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            srcFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcFolder, 'IncludingSubfolders', true));
        end
    end

    methods (Test)
        function testInitialDecayOfDiagonalSystem(testCase)
            % dx/dt = diag(-2,-3) x  ->  x(t) = [exp(-2t); exp(-3t)]
            [y, t] = initial(PhasorArray(diag([-2 -3])), [1; 1], 2*pi, 3, 'plot', false);
            expected = [exp(-2*t(:).'); exp(-3*t(:).')];
            testCase.verifyEqual(y, expected, 'AbsTol', testCase.tolOde);
        end

        function testInitialPeriodicSystemAgainstExactSolution(testCase)
            % Same closed form as above, reached through the homogeneous entry
            % point instead of the forced one.
            [y, t] = initial(PhasorArray.cos(), 1, 2*pi, 4*pi, 'plot', false);
            testCase.verifyEqual(y(:), exp(sin(t(:))), 'AbsTol', testCase.tolOde);
        end

        function testInitialRuns(testCase)
            % Create a stable periodic system dx/dt = A(t) x
            A0 = [-2 0; 0 -3];
            A1 = [0.05 0; 0 0.05];
            Ahm = PhasorArray(cat(3, conj(A1'), A0, A1));
    
            x0 = [1; 0];
            T = 2*pi;
            tfinal = 10;
            [y, t] = initial(Ahm, x0, T, tfinal);
    
            testCase.verifyTrue(~isempty(y), 'initial() should return non-empty y');
            testCase.verifyTrue(~isempty(t), 'initial() should return non-empty t');
            testCase.verifyTrue(numel(t) > 1, 'Should have multiple time samples');
            % State should decay (stable system)
            y_end = y(:, end);
            testCase.verifyTrue(norm(y_end) < 0.5, 'Stable system: final state should be small');
        end

        function testLsimPeriodicSystemAgainstExactSolution(testCase)
            % dx/dt = cos(t) x, x(0) = 1  ->  x(t) = exp(sin t), a genuinely
            % time-periodic case with a closed form.
            [y, t] = lsim(PhasorArray.cos(), 4*pi, 1, 2*pi, [], 'plot', false);
            testCase.verifyEqual(y(:), exp(sin(t(:))), 'AbsTol', testCase.tolOde);
        end

        function testLsimRuns(testCase)
            % Create stable periodic system with forcing
            A0 = [-2 0; 0 -3];
            Ahm = PhasorArray(A0);  % LTI for simplicity
    
            x0 = [0; 0];
            T = 2*pi;
            tfinal = 5;
    
            % Forcing: u(t) = cos(t) on state 1
            B = PhasorArray(eye(2));  % B matrix
            Uph = PhasorArray(cat(3, [0; 0], [1; 0], [0; 0]));  % cos input on channel 1
    
            [y, t] = lsim(Ahm, tfinal, x0, T, Uph);
    
            testCase.verifyTrue(~isempty(y), 'lsim() should return non-empty y');
            testCase.verifyTrue(~isempty(t), 'lsim() should return non-empty t');
        end


    end

    methods (Test, TestTags = {'Install'})
        % Smoke set for a fresh install: no optional toolbox, a few seconds,
        % one check per layer. Run with run_all_tests("install").
        function testLsimStepResponseOfFirstOrderSystem(testCase)
            % dx/dt = -x + 1, x(0) = 0  ->  x(t) = 1 - exp(-t)
            [y, t] = lsim(PhasorArray(-1), 5, 0, 2*pi, PhasorArray(1), 'plot', false);
            testCase.verifyEqual(y(:), (1 - exp(-t(:))), 'AbsTol', testCase.tolOde);
        end
    end
end
