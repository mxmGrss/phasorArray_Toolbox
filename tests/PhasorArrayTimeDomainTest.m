classdef PhasorArrayTimeDomainTest < matlab.unittest.TestCase
    %PHASORARRAYTIMEDOMAINTEST  The bijection between A(t) and its coefficients.
    %
    %   Everything defined by what the signal does over a period rather than by its
%   coefficient array: the trigonometric round trip with exact coefficients, the
%   comparison operators, abs, regularize, setHarmonic and the symmetry
%   projections.

    properties
        tol = 1e-10;
        tolOde = 1e-4;
        nBins = 7;   % 128 samples: above Nyquist for every harmonic used here
    end

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            srcFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcFolder, 'IncludingSubfolders', true));
        end
    end

    methods (Access = private)
        function c = coef(~, A, i, j, k)
            % Harmonic k of entry (i,j); slice k sits at k + h + 1.
            v = A.value;
            c = v(i, j, k + A.h + 1);
        end

        function A = build(testCase, func)
            A = PhasorArray.funcToPhasorArray(func, 2*pi, testCase.nBins);
        end

        function At = sample(~, func, th)
            cells = arrayfun(func, th, 'UniformOutput', false);
            At = cat(3, cells{:});
        end

    end

    methods (Test)
        function testAbsFoldsANegativeSignal(testCase)
            % |cos| has a corner, so the truncated series converges as 1/k^2 and
            % the tolerance is loose on purpose; the point is the folding.
            a = abs(PhasorArray.cos());
            th = [0.4, pi/2 + 0.3, pi, 3*pi/2 + 0.2];
            testCase.verifyEqual(squeeze(evalp(a, th)).', abs(cos(th)), 'AbsTol', 5e-3);
            testCase.verifyGreaterThan(real(a{1,1,0}), 0.5, ...
                '|cos| has a positive mean, unlike cos');
        end

        function testAbsOfStrictlyPositiveSignalIsIdentity(testCase)
            % 1.5 + cos never changes sign, so |A| = A and the check is exact.
            A = PhasorArray(1.5) + PhasorArray.cos();
            th = linspace(0, 2*pi, 65);
            testCase.verifyEqual(squeeze(evalp(abs(A), th)), squeeze(evalp(A, th)), ...
                'AbsTol', 1e-8, 'abs changed a positive signal');
        end

        function testComparisonExtraOutputsAgree(testCase)
            % r == (frac == 1), and both the indicator and its phasor form must
            % carry the same fraction: the DC term of a square wave is its duty.
            [r, frac, crenel, Cph] = ge(PhasorArray.cos(), PhasorArray(0));
            testCase.verifyEqual(r, all(crenel, 3), 'r disagrees with the indicator');
            testCase.verifyEqual(frac, mean(double(crenel), 3), 'AbsTol', testCase.tol);
            testCase.verifyEqual(real(Cph{1,1,0}), frac, 'AbsTol', 1e-8, ...
                'the DC term of the square wave is not the duty cycle');
        end

        function testComparisonFractionIsExact(testCase)
            % cos(th) < 1/2 holds on (pi/3, 5pi/3), a fraction 2/3 of the period.
            % The verdict is sampled, so the fraction carries the grid step.
            [r, frac] = lt(PhasorArray.cos(), PhasorArray(0.5));
            testCase.verifyFalse(r, 'the strict inequality does not hold everywhere');
            testCase.verifyEqual(frac, 2/3, 'AbsTol', 0.01, 'wrong fraction of the period');
        end

        function testComparisonHoldingEverywhereReportsTrue(testCase)
            [r, frac] = gt(PhasorArray.cos(), PhasorArray(-2));
            testCase.verifyTrue(r, 'cos(t) > -2 holds for every t');
            testCase.verifyEqual(frac, 1, 'AbsTol', testCase.tol);
        end

        function testComparisonRefusesAComplexDifference(testCase)
            % MATLAB's ordering silently compares real parts on complex input;
            % refusing is the whole point of the rewrite.
            A = PhasorArray(complex(randn(2,2,3), randn(2,2,3)));
            testCase.verifyError(@() lt(A, PhasorArray(0)), 'PhasorArray:compare:complexOrder');
        end

        function testComplexCoefficientsAreExact(testCase)
            f = @(th) 3*exp(1i*th) + 2i*exp(-2i*th) + (1 - 1i);
            A = testCase.build(f);

            testCase.verifyFalse(isrealp(A), 'A one-sided spectrum is not conjugate symmetric');
            testCase.verifyEqual(testCase.coef(A,1,1,+1),  3,      'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(A,1,1,-2), +2i,     'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(A,1,1, 0),  1 - 1i, 'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(A,1,1,-1),  0,      'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(A,1,1,+2),  0,      'AbsTol', testCase.tol);
        end

        function testComplexEvalpMatchesTheHandle(testCase)
            F = @(th) [exp(1i*th) + 2, 1i*exp(-3i*th); 0, exp(2i*th) - 1i];
            A = testCase.build(F);
            th = [0.137, 1.913, 4.221, 6.02];
            testCase.verifyEqual(evalp(A, th), testCase.sample(F, th), 'AbsTol', testCase.tol);
        end

        function testComplexProductShiftsHarmonics(testCase)
            % exp(i th) * exp(2i th) = exp(3i th): a single coefficient moves.
            A = testCase.build(@(th) exp(1i*th));
            B = testCase.build(@(th) exp(2i*th));
            P = A*B;
            testCase.verifyEqual(testCase.coef(P,1,1,+3), 1, 'AbsTol', testCase.tol);
            testCase.verifyEqual(max(abs(P.value), [], 'all'), 1, 'AbsTol', testCase.tol, ...
                'No other harmonic may be populated');
        end

        function testConjugateOfComplexSignal(testCase)
            % conj(f)(th) has coefficients conj(c_-k): conjugate and flip.
            f = @(th) 3*exp(1i*th) + 2i*exp(-2i*th);
            A = testCase.build(f);
            C = testCase.build(@(th) conj(f(th)));
            testCase.verifyEqual(C.value, conj(flip(A.value, 3)), 'AbsTol', testCase.tol);
        end

        function testConjugateSymmetryOfRealSignal(testCase)
            F = @(th) [1 + cos(th), sin(2*th); 2*sin(th), 3 - cos(3*th)];
            A = testCase.build(F);
            v = A.value;
            testCase.verifyEqual(v, conj(flip(v, 3)), 'AbsTol', testCase.tol, ...
                'A real periodic matrix must satisfy c_-k = conj(c_k)');
            testCase.verifyTrue(isrealp(A), 'isrealp should recognise the symmetry');
        end

        function testCtrretroIsTheRetroConjugateTranspose(testCase)
            % It returned A(t)' before, duplicating mctranspose. Ground truth is
            % evaluated at -t, so no algebraic restatement can hide the error.
            A = PhasorArray.random(2, 3, 3);
            th = [0.3, 1.7, 2.9, 4.8];
            expected = conj(permute(evalp(A, -th), [2 1 3]));
            testCase.verifyEqual(evalp(ctrretro(A), th), expected, 'AbsTol', testCase.tol);
        end

        function testEqualityIsExactOnIdenticalArrays(testCase)
            A = PhasorArray.random(2, 3, 2);
            testCase.verifyTrue(all(A == A, 'all'), 'A == A must hold entrywise');
            testCase.verifyFalse(all(A == A + PhasorArray.eye(2, 3), 'all'));
        end


        function testEvalpIsPeriodic(testCase)
            A = testCase.build(@(th) [cos(th), sin(2*th)]);
            th = [0.3, 1.7, 5.1];
            testCase.verifyEqual(evalp(A, th + 2*pi), evalp(A, th), 'AbsTol', testCase.tol);
        end

        function testEvalpMatchesTheHandleOffGrid(testCase)
            F = @(th) [1 + cos(th),  sin(2*th); 2*sin(th), 3 - cos(3*th)];
            A = testCase.build(F);

            % Deliberately off the sampling grid: reconstruction is exact for a
            % band-limited signal, not merely exact at the sample points.
            th = [0.137, 1.913, 2.5, 4.221, 6.02];
            testCase.verifyEqual(evalp(A, th), testCase.sample(F, th), 'AbsTol', testCase.tol);
        end

        function testForcedSymmetryOnComplexIsOptIn(testCase)
            % isReal defaults to isreal of the samples. Asking for the symmetry on
            % a complex signal still averages c_k with conj(c_-k), which is a
            % destructive request, not a silent default.
            f = @(th) 3*exp(1i*th);
            auto   = PhasorArray.funcToPhasorArray(f, 2*pi, testCase.nBins);
            forced = PhasorArray.funcToPhasorArray(f, 2*pi, testCase.nBins, 'isReal', true);

            testCase.verifyEqual(testCase.coef(auto,1,1,+1), 3, 'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(auto,1,1,-1), 0, 'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(forced,1,1,+1), 1.5, 'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(forced,1,1,-1), 1.5, 'AbsTol', testCase.tol);
        end

        function testIndependentDeclarationsMatch(testCase)
            f1 = @(th) 4*cos(2*th) + sin(th) - 1;
            f2 = @(th) -1 + sin(th) + 4*cos(2*th);      % same signal, written differently
            f3 = @(th) -1 + sin(th) + 4*cos(2*th + 0.1); % phase shifted: must differ

            A = testCase.build(f1);
            B = testCase.build(f2);
            C = testCase.build(f3);

            testCase.verifyEqual(A.value, B.value, 'AbsTol', testCase.tol);
            testCase.verifyGreaterThan(max(abs(A.value - C.value), [], 'all'), 1e-3, ...
                'A phase shift must move the coefficients');
        end

        function testMatrixCoefficientsAreExact(testCase)
            F = @(th) [1 + cos(th),  sin(2*th); ...
                       2*sin(th),    3 - cos(3*th)];
            A = testCase.build(F);

            testCase.verifyEqual(A.h, 3, 'Bandwidth is set by the 3rd harmonic entry');

            % (1,1) = 1 + cos(th)
            testCase.verifyEqual(testCase.coef(A,1,1, 0),  1,     'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(A,1,1,+1),  0.5,   'AbsTol', testCase.tol);
            % (1,2) = sin(2 th)
            testCase.verifyEqual(testCase.coef(A,1,2,+2), -0.5i,  'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(A,1,2,-2), +0.5i,  'AbsTol', testCase.tol);
            % (2,1) = 2 sin(th)
            testCase.verifyEqual(testCase.coef(A,2,1,+1), -1i,    'AbsTol', testCase.tol);
            % (2,2) = 3 - cos(3 th)
            testCase.verifyEqual(testCase.coef(A,2,2, 0),  3,     'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(A,2,2,+3), -0.5,   'AbsTol', testCase.tol);

            % Entries below their own bandwidth must be padded with exact zeros.
            testCase.verifyEqual(testCase.coef(A,1,1,+3), 0, 'AbsTol', testCase.tol);
        end

        function testMatrixProductAgainstTimeDomain(testCase)
            F = @(th) [1 + cos(th),  sin(2*th); 2*sin(th), 3 - cos(3*th)];
            G = @(th) [cos(2*th),    1 + sin(th); -sin(3*th), cos(th) + 2];

            A = testCase.build(F);
            B = testCase.build(G);
            P = A*B;
            testCase.verifyEqual(P.h, 6, 'Convolution adds the bandwidths');

            % The harmonic product must agree with the pointwise matrix product,
            % which is the whole point of the Toeplitz representation.
            th = linspace(0, 2*pi, 41);
            th(end) = [];
            testCase.verifyEqual(evalp(P, th), pagemtimes(evalp(A, th), evalp(B, th)), ...
                'AbsTol', 1e-11);

            % And with the FFT of the product taken in the time domain.
            Pdirect = testCase.build(@(th) F(th)*G(th));
            testCase.verifyEqual(evalp(P, th), evalp(Pdirect, th), 'AbsTol', 1e-11);
        end

        function testOversamplingDoesNotChangeCoefficients(testCase)
            % Any sampling above Nyquist gives the same coefficients; below it,
            % the 5th harmonic folds and the result is a different signal.
            f = @(th) cos(5*th) + sin(th);
            A = PhasorArray.funcToPhasorArray(f, 2*pi, 5);
            B = PhasorArray.funcToPhasorArray(f, 2*pi, 8);
            testCase.verifyEqual(A.h, 5);
            testCase.verifyEqual(A.value, B.value, 'AbsTol', testCase.tol);
        end

        function testParityOfPureSineAndCosine(testCase)
            C = testCase.build(@(th) cos(3*th));
            S = testCase.build(@(th) sin(3*th));

            % cos is even: coefficients real and symmetric.
            testCase.verifyEqual(testCase.coef(C,1,1,+3), 0.5, 'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(C,1,1,-3), 0.5, 'AbsTol', testCase.tol);

            % sin is odd: coefficients purely imaginary and antisymmetric.
            testCase.verifyEqual(testCase.coef(S,1,1,+3), -0.5i, 'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(S,1,1,-3), +0.5i, 'AbsTol', testCase.tol);
            testCase.verifyLessThan(abs(real(squeeze(S.value))), testCase.tol, ...
                'An odd real signal has no real part in its coefficients');
        end

        function testPhasorSymmetryProjectsAndIsIdempotent(testCase)
            A = PhasorArray(complex(randn(2,2,5), randn(2,2,5)));
            B = phasorSymmetry(A, "real");
            testCase.verifyTrue(isrealp(B), 'the projection did not produce a real signal');
            testCase.verifyEqual(phasorSymmetry(B, "real").value, B.value, 'AbsTol', testCase.tol, ...
                'a projection must be idempotent');
        end

        function testProductIsNotCommutative(testCase)
            F = @(th) [cos(th), 1; 0, sin(2*th)];
            G = @(th) [0, sin(th); cos(3*th), 1];
            A = testCase.build(F);
            B = testCase.build(G);
            AB = A*B;
            BA = B*A;
            testCase.verifyGreaterThan(max(abs(AB.value - BA.value), [], 'all'), 1e-3, ...
                'Matrix phasor product must keep the ordering of the time-domain product');
        end

        function testRegularizeDampsHigherHarmonicsMonotonically(testCase)
            A = PhasorArray.random(2, 2, 4);
            widths = [0.1, 0.3, 0.6];
            tops = zeros(size(widths));
            for k = 1:numel(widths)
                Rk = regularize(A, widths(k));
                tops(k) = norm(Rk{:,:,4});
            end
            testCase.verifyLessThan(tops(1), norm(A{:,:,4}), ...
                'mollifying must damp the top harmonic');
            testCase.verifyTrue(all(diff(tops) < 0), ...
                sprintf('damping is not monotone in the width: %s', mat2str(tops, 3)));
        end

        function testRegularizeLeavesTheMeanUntouched(testCase)
            % The mollifier has unit mass, so phi_hat(0) = 1 and the DC is fixed.
            A = PhasorArray.random(2, 2, 3);
            R = regularize(A, 0.3);
            testCase.verifyEqual(R{:,:,0}, A{:,:,0}, 'AbsTol', testCase.tol);
        end

        function testRegularizeTendsToTheIdentityAsWidthVanishes(testCase)
            A = PhasorArray.random(2, 2, 3);
            R = regularize(A, 1e-9);
            testCase.verifyEqual(R.value, A.value, 'AbsTol', 1e-9);
        end


        function testScalarProductMatchesProductToSum(testCase)
            % cos(th)*sin(2 th) = ( sin(3 th) + sin(th) ) / 2
            A = testCase.build(@(th) cos(th));
            B = testCase.build(@(th) sin(2*th));
            P = A*B;

            testCase.verifyEqual(testCase.coef(P,1,1,+3), -0.25i, 'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(P,1,1,-3), +0.25i, 'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(P,1,1,+1), -0.25i, 'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(P,1,1,-1), +0.25i, 'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(P,1,1, 0),  0,     'AbsTol', testCase.tol);
        end

        function testSetHarmonicAcceptsEveryIndexForm(testCase)
            A = PhasorArray.random(2, 2, 2);
            whole = setHarmonic(A, 1, ones(2));
            entry = setHarmonic(A, [1 2 1], 5);
            braced = setHarmonic(A, {1, 2, 1}, 7);

            testCase.verifyEqual(whole{:,:,1}, ones(2), 'AbsTol', testCase.tol);
            testCase.verifyEqual(entry{1,2,1}, 5, 'AbsTol', testCase.tol);
            testCase.verifyEqual(braced{1,2,1}, 7, 'AbsTol', testCase.tol);
            % Untouched harmonics must stay put.
            testCase.verifyEqual(entry{:,:,2}, A{:,:,2}, 'AbsTol', testCase.tol);
        end

        function testSetHarmonicMirrorsWhenRealnessIsRequested(testCase)
            % isReal defaults to isreal(A): writing +k must mirror onto -k so
            % the array stays a real signal.
            A = PhasorArray(zeros(1,1,5));
            B = setHarmonic(A, 1, 2 + 3i);
            testCase.verifyEqual(B{1,1,-1}, conj(B{1,1,1}), 'AbsTol', testCase.tol, ...
                'the conjugate mirror was not written');
            testCase.verifyTrue(isrealp(B), 'the result should still be a real signal');
        end

        function testSinCosFormEvaluatesToTheSameSignal(testCase)
            % SinCosForm re-encodes the same real signal on a [sin; cos] basis,
            % so evaluating it in that form must return the original samples.
            A = PhasorArray.random(2, 3, 3);
            th = [0.21, 1.4, 2.9, 5.05];
            direct = evalp(A, th);
            viaSinCos = PhasorArray2time(A.SinCosForm(), 2*pi, th, ...
                'providedPhasorForm', "SinCos");
            testCase.verifyEqual(viaSinCos, direct, 'AbsTol', testCase.tol, ...
                'the SinCos form does not describe the same signal');
        end

        function testSquareOfSineMatchesHalfAngle(testCase)
            % sin^2(th) = (1 - cos(2 th))/2
            S = testCase.build(@(th) sin(th));
            P = S*S;

            testCase.verifyEqual(testCase.coef(P,1,1, 0),  0.5,  'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(P,1,1,+2), -0.25, 'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(P,1,1,-2), -0.25, 'AbsTol', testCase.tol);
        end

        function testSymmetryClosureReportsTheRequestedProperty(testCase)
            held = symmetryClosure("real");
            testCase.verifyTrue(ismember("real", held), ...
                'the closure of {real} must contain real');
        end

    end

    methods (Test, TestTags = {'Install'})
        % Smoke set for a fresh install: no optional toolbox, a few seconds,
        % one check per layer. Run with run_all_tests("install").
        function testScalarCoefficientsAreExact(testCase)
            f = @(th) 2 + 3*cos(th) - 5*sin(2*th);
            A = testCase.build(f);

            testCase.verifyEqual(A.h, 2, 'reduce should prune to the true bandwidth');
            testCase.verifyEqual(testCase.coef(A,1,1, 0),  2,        'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(A,1,1,+1),  1.5,      'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(A,1,1,-1),  1.5,      'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(A,1,1,+2), +2.5i,     'AbsTol', testCase.tol);
            testCase.verifyEqual(testCase.coef(A,1,1,-2), -2.5i,     'AbsTol', testCase.tol);
        end

        function testEvalp(testCase)
            c = PhasorArray.cos();
            % cos(0) = 1
            testCase.verifyTrue(abs(evalp(c, 0) - 1) < testCase.tol, 'cos(0)');
            % cos(pi) = -1
            testCase.verifyTrue(abs(evalp(c, pi) - (-1)) < testCase.tol, 'cos(pi)');
            % cos(pi/2) = 0
            testCase.verifyTrue(abs(evalp(c, pi/2)) < testCase.tol, 'cos(pi/2)');
        end
    end
end
