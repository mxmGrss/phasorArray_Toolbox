classdef PhasorArrayCompatibilityTest < matlab.unittest.TestCase
    %PHASORARRAYCOMPATIBILITYTEST  Branches a recent release never takes.
    %
    %   The toolbox supports R2021b upward, so several kernels carry a second
%   implementation the version gate never selects on a modern MATLAB. A mutation
%   campaign found breaking them turned no test red. Each path here is reachable
%   by direct call or name-value pair, so these run on any release.

    properties
        tol = 1e-10;
    end

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            srcFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcFolder, 'IncludingSubfolders', true));
        end
    end

    methods (Test)
        function testEvaluationMethodsAgreeComplex(testCase)
            % forceReal false selects the exponential basis; the real sin/cos
            % basis above never exercises it.
            A = PhasorArray(complex(randn(2,2,5), randn(2,2,5)));
            th = linspace(0, 2*pi, 33);
            th(end) = [];
            ref = PhasorArray2time(A, 2*pi, th, 'computationMethod', 1, 'forceReal', false);
            for method = [2 3]
                alt = PhasorArray2time(A, 2*pi, th, 'computationMethod', method, 'forceReal', false);
                testCase.verifyEqual(alt, ref, 'AbsTol', testCase.tol, ...
                    sprintf('complex computationMethod %d differs from 1', method));
            end
        end

        function testEvaluationMethodsAgreeReal(testCase)
            A = PhasorArray.random(3, 2, 4);
            th = linspace(0, 2*pi, 33);
            th(end) = [];
            ref = PhasorArray2time(A, 2*pi, th, 'computationMethod', 1);
            for method = [2 3]
                alt = PhasorArray2time(A, 2*pi, th, 'computationMethod', method);
                testCase.verifyEqual(alt, ref, 'AbsTol', testCase.tol, ...
                    sprintf('computationMethod %d differs from 1', method));
            end
        end


        function testFallbackConvolutionOnComplexOperands(testCase)
            % Conjugate symmetry is what the two branches index differently, so
            % a signal without it is the discriminating case.
            v = complex(randn(2,2,5), randn(2,2,5));
            w = complex(randn(2,2,5), randn(2,2,5));
            fast = PhasorArrayTimes(PhasorArray(v), PhasorArray(w));
            slow = PhasorArrayTimes2(v, w);
            if isa(fast, 'PhasorArray'), fast = fast.Value; end
            if isa(slow, 'PhasorArray'), slow = slow.Value; end
            testCase.verifyEqual(slow, fast, 'AbsTol', testCase.tol);
        end

        function testRealAndComplexBasisAgreeOnRealSignal(testCase)
            % A real signal may go through either basis and must land on the
            % same samples; only the round-off differs.
            A = PhasorArray.random(2, 2, 3);
            th = [0.31, 1.77, 3.9, 5.5];
            viaReal    = PhasorArray2time(A, 2*pi, th, 'forceReal', true);
            viaComplex = PhasorArray2time(A, 2*pi, th, 'forceReal', false);
            testCase.verifyEqual(real(viaComplex), viaReal, 'AbsTol', testCase.tol);
            testCase.verifyLessThan(max(abs(imag(viaComplex)), [], 'all'), testCase.tol, ...
                'A real signal must reconstruct with no imaginary part');
        end

        function testSparseToeplitzMatchesDense(testCase)
            A = PhasorArray.random(2, 2, 3);
            B = PhasorArray.random(2, 2, 3);
            h = 5;
            dense = A.TBmtimes(B, h);
            sparseM = A.spTBmtimes(B, h);
            testCase.verifyTrue(issparse(sparseM), 'spTBmtimes should stay sparse');
            testCase.verifyEqual(full(sparseM), dense, 'AbsTol', testCase.tol);
        end

    end

    methods (Test, TestTags = {'Install'})
        % Smoke set for a fresh install: no optional toolbox, a few seconds,
        % one check per layer. Run with run_all_tests("install").
        function testFallbackConvolutionMatchesTensorprod(testCase)
            for shape = {[2 2 2], [3 1 4], [4 3 1]}
                n = shape{1}(1); m = shape{1}(2); h = shape{1}(3);
                A = PhasorArray.random(n, m, h);
                B = PhasorArray.random(m, n, h);

                fast = PhasorArrayTimes(A, B);
                slow = PhasorArrayTimes2(A.Value, B.Value);
                if isa(fast, 'PhasorArray'), fast = fast.Value; end
                if isa(slow, 'PhasorArray'), slow = slow.Value; end

                testCase.verifyEqual(size(slow), size(fast), ...
                    sprintf('shape mismatch on %dx%dx%d', n, m, h));
                testCase.verifyEqual(slow, fast, 'AbsTol', testCase.tol, ...
                    sprintf('fallback convolution differs on %dx%dx%d', n, m, h));
            end
        end
    end
end
