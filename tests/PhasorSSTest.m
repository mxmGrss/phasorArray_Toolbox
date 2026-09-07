classdef PhasorSSTest < matlab.unittest.TestCase
    %PHASORSSTEST  PhasorSS: dimensions, static/real predicates, angle evaluation.
    %
    %   PhasorSS had no dedicated test class before this one -- it was only
    %   exercised through Exemples/ and templates/, which count as coverage
    %   for the constructor and most operators but not for these accessors.

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

        function testNxNuNyMatchTheMatrices(testCase)
            % nx/nu/ny read off size(A), size(B,2), size(C,1) -- a 2-state,
            % 1-input, 1-output system pins all three independently, so a
            % transposed or swapped dimension would be caught.
            A = PhasorArray(cat(3,[0 1;-2 -3],[1/2 0;0 0]),"isreal",true);
            B = [0;1]; C = [1 0]; D = 0;
            P = PhasorSS(A,B,C,D,0.1,'isReal',true);
            testCase.verifyEqual(P.nx, 2, 'nx must be the number of states');
            testCase.verifyEqual(P.nu, 1, 'nu must be the number of inputs');
            testCase.verifyEqual(P.ny, 1, 'ny must be the number of outputs');
        end

        function testNxNuNyOnAMimoSystem(testCase)
            % Same, with nx ~= nu ~= ny ~= each other so no accessor could
            % accidentally return another one's value and still pass.
            A = PhasorArray(-eye(3));
            B = PhasorArray(randn(3,2));
            C = PhasorArray(randn(4,3));
            D = PhasorArray(randn(4,2));
            P = PhasorSS(A,B,C,D,1);
            testCase.verifyEqual(P.nx, 3, 'nx');
            testCase.verifyEqual(P.nu, 2, 'nu');
            testCase.verifyEqual(P.ny, 4, 'ny');
        end

        function testIsStaticTrueOnlyWhenAIsEmpty(testCase)
            % isStatic checks isempty(obj.A) verbatim -- pin both branches so
            % a future refactor cannot silently invert the condition.
            dynamic = PhasorSS(PhasorArray(-1), 1, 1, 0, 1);
            testCase.verifyFalse(dynamic.isStatic, 'a system with a state matrix is not static');

            static = PhasorSS([], [], [], 2, 1);
            testCase.verifyTrue(static.isStatic, 'a system with empty A must be reported static');
        end

        function testRealifyProducesTheRealPartOfEachMatrix(testCase)
            % realify must equal mreal() applied to A, B, C, D individually --
            % not, say, drop the imaginary part globally or touch only A.
            A = PhasorArray((1+0.3i) * [0 1; -2 -3]);
            B = PhasorArray((0.5-0.2i) * [0; 1]);
            C = PhasorArray((1+1i) * [1 0]);
            D = PhasorArray(0.1i);
            P  = PhasorSS(A,B,C,D,1);
            testCase.verifyFalse(P.isReal, 'the fixture must be complex for this test to mean anything');

            Pr = P.realify();
            testCase.verifyTrue(Pr.isReal, 'realify must flag the result as real');
            testCase.verifyEqual(value(Pr.A), value(mreal(A)), 'AbsTol', testCase.tol, 'A');
            testCase.verifyEqual(value(Pr.B), value(mreal(B)), 'AbsTol', testCase.tol, 'B');
            testCase.verifyEqual(value(Pr.C), value(mreal(C)), 'AbsTol', testCase.tol, 'C');
            testCase.verifyEqual(value(Pr.D), value(mreal(D)), 'AbsTol', testCase.tol, 'D');
        end

        function testEvalAngleMatchesEvalpAtEveryMatrix(testCase)
            % evalAngle's argument IS the PhasorArray time/angle variable --
            % it is evalp under another name, applied to all four matrices at
            % once. A, B, C, D are given the same h on purpose: evalAngle
            % combines them with one shared harmonic-combination vector sized
            % off A. Each is also given a nonzero cosine coefficient so none
            % gets silently reduced to h=0 by the constructor -- when that
            % happens, evalAngle's harmonicCombine degenerates into a scalar
            % broadcast and returns the wrong shape (see the flagged defect
            % in TODO-architecture.md); testing that path is out of scope
            % here, this test only pins the case where h truly matches.
            A = PhasorArray(cat(3,[0 1;-2 -3],[1/2 0;0 0]),   "isreal", true);
            B = PhasorArray(cat(3,[0;1],[0.01;0]),            "isreal", true);
            C = PhasorArray(cat(3,[1 0],[0.01 0]),            "isreal", true);
            D = PhasorArray(cat(3,0.1,0.05),                  "isreal", true);
            P = PhasorSS(A,B,C,D,0.1,'isReal',true);

            for theta = [0, 0.7, pi/3, 2.1]
                [Ae,Be,Ce,De] = P.evalAngle(theta);
                testCase.verifyEqual(Ae, evalp(A,theta), 'AbsTol', testCase.tol, sprintf('A at theta=%g', theta));
                testCase.verifyEqual(Be, evalp(B,theta), 'AbsTol', testCase.tol, sprintf('B at theta=%g', theta));
                testCase.verifyEqual(Ce, evalp(C,theta), 'AbsTol', testCase.tol, sprintf('C at theta=%g', theta));
                testCase.verifyEqual(De, evalp(D,theta), 'AbsTol', testCase.tol, sprintf('D at theta=%g', theta));
            end
        end

        function testAddOutputAppendsRowsWithoutTouchingTheExisting(testCase)
            % addOutput is a concatenation: C gains rows, D gains rows, A and
            % B do not move. The old rows must survive byte for byte -- that
            % is the part a naive reimplementation (e.g. rebuilding C from
            % scratch) could get wrong while still passing a size check.
            A = PhasorArray(-eye(2));
            B = PhasorArray([1 0; 0 1]);
            C = PhasorArray([1 0]);
            D = PhasorArray([0 0]);
            P  = PhasorSS(A,B,C,D,1);
            P2 = P.addOutput([0 1], [0 0]);

            testCase.verifyEqual(P2.ny, P.ny + 1, 'ny must grow by the number of new rows');
            testCase.verifyEqual(value(P2.A), value(P.A), 'A must be untouched');
            testCase.verifyEqual(value(P2.B), value(P.B), 'B must be untouched');
            testCase.verifyEqual(value(P2.C(1:P.ny,:)), value(P.C), 'the old rows of C must be preserved exactly');
            testCase.verifyEqual(value(P2.D(1:P.ny,:)), value(P.D), 'the old rows of D must be preserved exactly');
            testCase.verifyEqual(value(P2.C(P.ny+1:end,:)), [0 1], 'the new row of C');
            testCase.verifyEqual(value(P2.D(P.ny+1:end,:)), [0 0], 'the new row of D');
        end

        function testAddInputAppendsColumnsWithoutTouchingTheExisting(testCase)
            % Symmetric to addOutput: B gains columns, D gains columns, A and
            % C are untouched, and the old columns keep their values.
            A = PhasorArray(-eye(2));
            B = PhasorArray([1; 0]);
            C = PhasorArray([1 0; 0 1]);
            D = PhasorArray([0; 0]);
            P  = PhasorSS(A,B,C,D,1);
            P2 = P.addInput([0; 1], [0; 0]);

            testCase.verifyEqual(P2.nu, P.nu + 1, 'nu must grow by the number of new columns');
            testCase.verifyEqual(value(P2.A), value(P.A), 'A must be untouched');
            testCase.verifyEqual(value(P2.C), value(P.C), 'C must be untouched');
            testCase.verifyEqual(value(P2.B(:,1:P.nu)), value(P.B), 'the old columns of B must be preserved exactly');
            testCase.verifyEqual(value(P2.D(:,1:P.nu)), value(P.D), 'the old columns of D must be preserved exactly');
            testCase.verifyEqual(value(P2.B(:,P.nu+1:end)), [0;1], 'the new column of B');
            testCase.verifyEqual(value(P2.D(:,P.nu+1:end)), [0;0], 'the new column of D');
        end

        function testHmqDcGainMatchesTheClassicFormula(testCase)
            % For an LTI system (h=0 throughout), the harmonic DC gain must
            % reduce to the ordinary state-space DC gain -C*A^-1*B + D.
            % Compared against a reference built independently with plain
            % numeric linear algebra, not against the toolbox's own
            % dcgain/toeplitzSS machinery.
            %
            % h_compute=0 exercises a real defect fixed the same session:
            % extractBlocksAndLabels used to reject h=0 outright
            % (mustBePositive instead of mustBeNonnegative), even though it
            % is the simplest valid case -- a single-harmonic block.
            % formatInputRange had a second, independent defect: it silently
            % dropped every selector that was not a channel name (':', a
            % numeric index, and the paired harmonic list itself), so
            % hmqDcGain returned an empty [0x0] matrix regardless of range.
            A = [-2 1; 0 -3]; B = [0;1]; C = [1 0]; D = 0.5;
            P = PhasorSS(PhasorArray(A), PhasorArray(B), PhasorArray(C), PhasorArray(D), 2*pi);
            expected = -C*(A\B) + D;
            [dcGain, ~] = P.hmqDcGain(0, 2*pi);
            testCase.verifyEqual(dcGain, expected, 'AbsTol', 1e-8, ...
                'hmqDcGain must reduce to -C*A^-1*B + D for an LTI system');
        end

        function testHmqBodeMatchesFreqrespOfTheKnownLti(testCase)
            % Same idea for the frequency response: at h=0 the harmonic Bode
            % is exactly the ordinary transfer function's freqresp.
            openBefore = findobj('Type','figure');
            closer = onCleanup(@() close(setdiff(findobj('Type','figure'), openBefore))); %#ok<NASGU>

            A = [-2 1; 0 -3]; B = [0;1]; C = [1 0]; D = 0.5;
            P = PhasorSS(PhasorArray(A), PhasorArray(B), PhasorArray(C), PhasorArray(D), 2*pi);
            sysRef = ss(A,B,C,D);
            freqs = [0.1 1 5];
            [Hsel, ~, ~] = P.HmqBode(0, 2*pi, 0, 'freqRange', freqs);
            testCase.verifyEqual(size(findobj('Type','figure')), size(openBefore), ...
                'HmqBode must not open a figure when its outputs are requested (nargout>0)');
            Href = squeeze(freqresp(sysRef, 2*pi*freqs));
            testCase.verifyEqual(squeeze(Hsel), Href, 'AbsTol', 1e-6, ...
                'HmqBode at h=0 must match the plain freqresp of the same LTI system');
        end

        function testLftMatchesTheClassicFeedbackFormula(testCase)
            % Reference independent of the code under test: for sys2 a pure
            % static gain K (no states), the lower LFT reduces to the
            % textbook feedback-interconnection formulas
            %   Acl = A + B2*K*(I-D22*K)^-1*C2,  Bcl = B1 + B2*K*(I-D22*K)^-1*D21
            %   Ccl = C1 + D12*K*(I-D22*K)^-1*C2, Dcl = D11 + D12*K*(I-D22*K)^-1*D21
            % sys1: 2 inputs [w1 u1], 2 outputs [z1 y1], each scalar.
            A = [-1 0.5; 0 -2]; B = [1 0; 0 1]; C = [1 0; 0 1]; D = zeros(2,2);
            sys1 = PhasorSS(PhasorArray(A), PhasorArray(B), PhasorArray(C), PhasorArray(D), 1);
            K = 0.3;
            sys2 = PhasorSS([], [], [], PhasorArray(K), 1);   % static gain, no states

            closedLoop = lft(sys1, sys2, 1, 1);

            B1 = B(:,1); B2 = B(:,2);
            C1 = C(1,:); C2 = C(2,:);
            D11 = D(1,1); D12 = D(1,2); D21 = D(2,1); D22 = D(2,2);
            invTerm = K / (1 - D22*K);
            Acl = A + B2*invTerm*C2;
            Bcl = B1 + B2*invTerm*D21;
            Ccl = C1 + D12*invTerm*C2;
            Dcl = D11 + D12*invTerm*D21;

            testCase.verifyEqual(value(closedLoop.A), Acl, 'AbsTol', 1e-8, 'Acl');
            testCase.verifyEqual(value(closedLoop.B), Bcl, 'AbsTol', 1e-8, 'Bcl');
            testCase.verifyEqual(value(closedLoop.C), Ccl, 'AbsTol', 1e-8, 'Ccl');
            testCase.verifyEqual(value(closedLoop.D), Dcl, 'AbsTol', 1e-8, 'Dcl');
        end

        function testEvalAngleHandlesAReducedMatrix(testCase)
            % Regression test for a real defect found and fixed in
            % harmonicCombine (see its docstring): when a feedthrough matrix
            % is constant -- or otherwise gets reduced to h=0 by the
            % constructor -- while A has h>0, evalAngle used to hand a 1x1 D
            % to harmonicCombine alongside a [2h+1 x 1] combination vector.
            % MATLAB's scalar-broadcast rule for '*' made that call succeed
            % silently with the wrong shape ([1,1,2h+1]) instead of erroring.
            % harmonicCombine now center-truncates the mismatched side
            % instead, so this must come back as a plain scalar equal to
            % evalp(D, angle) -- D being constant, the same value at every
            % angle.
            A = PhasorArray(cat(3,[0 1;-2 -3],[1/2 0;0 0]), "isreal", true);
            B = PhasorArray(cat(3,[0;1],[0;0]),             "isreal", true);
            C = PhasorArray(cat(3,[1 0],[0 0]),             "isreal", true);
            D = PhasorArray(0.3);   % constructed directly at h=0, unlike A/B/C above
            P = PhasorSS(A,B,C,D,0.1,'isReal',true);
            testCase.verifyEqual(P.D.h, 0, 'fixture assumption: D must be reduced to h=0');

            [~,~,~,De] = P.evalAngle(0.7);
            testCase.verifySize(De, [1 1], 'De must come back as a scalar, not a per-harmonic broadcast');
            testCase.verifyEqual(De, evalp(D, 0.7), 'AbsTol', testCase.tol);
        end

    end
end
