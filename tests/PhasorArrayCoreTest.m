classdef PhasorArrayCoreTest < matlab.unittest.TestCase
    %PHASORARRAYCORETEST  Construction, indexing, arithmetic and payload type.
    %
    %   Constructors and properties, paren and brace indexing, the arithmetic and
%   concatenation operators, and the rule that a sym or sdpvar payload survives
%   every one of them. A cast back to double leaves the wrapper intact and only
%   the payload class reveals it.

    properties
        tol = 1e-10;
    end

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            srcFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcFolder, 'IncludingSubfolders', true));
        end
    end

    methods (Access = private)
        function needSymbolic(testCase)
            testCase.assumeTrue(logical(exist('sym', 'file')), 'Symbolic Math Toolbox required');
        end

        function needYalmip(testCase)
            testCase.assumeTrue(exist('sdpvar', 'file') == 2, 'YALMIP required');
        end

        function verifySym(testCase, x, label)
            % Taking the operand as an argument sidesteps MATLAB refusing dot
            % indexing straight after a parenthesised expression.
            testCase.verifyClass(x.value, 'sym', sprintf('%s lost its sym payload', label));
        end

        function verifySdp(testCase, x, label)
            % ndsdpvar and sdpvar are the same payload family at different ranks.
            testCase.verifyTrue(isa(x.value, 'ndsdpvar') || isa(x.value, 'sdpvar'), ...
                sprintf('%s lost its sdpvar payload, got %s', label, class(x.value)));
        end

    end

    methods (Test)

        function testAddPad(testCase)
            A = PhasorArray.random(2, 2, 3);  % h=3
            B = PhasorArray.random(2, 2, 5);  % h=5
            C = A + B;
            testCase.verifyTrue(C.h == 5, 'Result should have h=max(3,5)=5');
            % A's h=4,5 harmonics are zero, so C{:,:,4} = B{:,:,4}
            testCase.verifyTrue(max(abs(C{:, :,4} - B{:,:,4}), [], 'all') < testCase.tol, 'Padding mismatch');
        end

        function testBlkdiag(testCase)
            A = PhasorArray.random(2, 2, 3);
            B = PhasorArray.random(3, 3, 5);
            C = blkdiag(A, B);
            testCase.verifyTrue(size(C, 1) == 5, 'Rows = 2+3');
            testCase.verifyTrue(size(C, 2) == 5, 'Cols = 2+3');
        end

        function testConcatenationMixedH(testCase)
            A = PhasorArray.random(2, 2, 1);
            B = PhasorArray.random(2, 2, 3);
            
            % Horzcat
            H = [A, B];
            testCase.verifyTrue(size(H, 1) == 2 && size(H, 2) == 4, 'horzcat size wrong');
            testCase.verifyTrue(H.h == 3, 'horzcat should pad to h=3');
            
            % Vertcat
            V = [A; B];
            testCase.verifyTrue(size(V, 1) == 4 && size(V, 2) == 2, 'vertcat size wrong');
            testCase.verifyTrue(V.h == 3, 'vertcat should pad to h=3');
            
            % Blkdiag
            D = blkdiag(A, B);
            testCase.verifyTrue(size(D, 1) == 4 && size(D, 2) == 4, 'blkdiag size wrong');
            testCase.verifyTrue(D.h == 3, 'blkdiag should pad to h=3');
            
            % Verify preservation of off-diagonal zeros in blkdiag
            valD = squeeze(D.value);
            off_diag = valD(1:2, 3:4, :);
            testCase.verifyTrue(max(abs(off_diag(:))) == 0, 'blkdiag off-diagonal should be zero across all harmonics');
        end


        function testConstructorEye(testCase)
            A = PhasorArray.eye(3);
            testCase.verifyTrue(isequal(size(A, 1), 3), 'Size mismatch');
            val = A{1:3, 1:3, 0};  % DC component
            err = max(abs(val - eye(3)), [], 'all');
            testCase.verifyTrue(err < testCase.tol, 'DC should be identity');
        end

        function testConstructorIsreal(testCase)
            V = randn(2, 2, 4);  % h will be interpreted as positive-only phasors
            A = PhasorArray(V, 'isreal', true);
            testCase.verifyTrue(isreal(A), 'Should be flagged as real');
            % Evaluate: should produce real values
            At = evalp(A, pi/4);
            testCase.verifyTrue(max(abs(imag(At(:)))) < testCase.tol, 'Evaluation should be real');
        end

        function testConstructorRandom(testCase)
            A = PhasorArray.random(4, 3, 7);
            testCase.verifyTrue(isequal(size(A, 1), 4), 'Row mismatch');
            testCase.verifyTrue(isequal(size(A, 2), 3), 'Col mismatch');
            testCase.verifyTrue(A.h == 7, 'Harmonic order mismatch');
        end

        function testConstructorScalar(testCase)
            v = [0.5, 1, 0.5];   % DC=1, h1=0.5
            A = ScalarPhasorArray(v);
            testCase.verifyTrue(A.h == 1, 'h should be 1');
            % DC coefficient should be 1
            dc = A{1,1,0};
            testCase.verifyTrue(abs(dc - 1) < testCase.tol, 'DC mismatch');
        end

        function testConstructorTrig(testCase)
            c = PhasorArray.cos();
            s = PhasorArray.sin();
            testCase.verifyTrue(c.h == 1, 'cos should have h=1');
            testCase.verifyTrue(s.h == 1, 'sin should have h=1');
            % Evaluate at t=0: cos(0)=1, sin(0)=0
            c0 = evalp(c, 0);
            s0 = evalp(s, 0);
            testCase.verifyTrue(abs(c0 - 1) < testCase.tol, 'cos(0) should be 1');
            testCase.verifyTrue(abs(s0) < testCase.tol, 'sin(0) should be 0');
        end

        function testConstructorZero(testCase)
            A = PhasorArray.zeros(3, 4);
            testCase.verifyTrue(isequal(size(A, 1), 3), 'Row count');
            testCase.verifyTrue(isequal(size(A, 2), 4), 'Col count');
            testCase.verifyTrue(A.h == 0, 'Zero matrix should have h=0');
        end

        function testFloatFastPathSymmetry(testCase)
            A = PhasorArray.random(2, 2, 2);
            % A + A' should automatically be recognized as symmetric and real (if A is real)
            S = A + A';
            % isreal() native override
            testCase.verifyTrue(isreal(S), 'isreal() method should return true');
        end

        function testHermitian(testCase)
            A = PhasorArray.random(3, 2, 5);
            Ah = A';
            testCase.verifyTrue(size(Ah, 1) == 2 && size(Ah, 2) == 3, 'Hermitian size');
            % For Hermitian: (A')_k = conj(A_{-k})'
            % At k=0: conj(A_0)' 
            testCase.verifyTrue(max(abs(Ah{:, :,0} - conj(A{:,:,0}).'), [], 'all') < testCase.tol, 'Hermitian DC');
        end

        function testHorzcat(testCase)
            A = PhasorArray.random(2, 2, 3);
            B = PhasorArray.random(2, 3, 5);
            C = [A, B];
            testCase.verifyTrue(size(C, 1) == 2, 'Rows');
            testCase.verifyTrue(size(C, 2) == 5, 'Cols = 2+3');
            testCase.verifyTrue(C.h == 5, 'h = max(3,5)');
        end

        function testIndexingBracesAbsent(testCase)
            A = PhasorArray.random(2, 2, 3);  % h=3
            % Request harmonic far beyond stored range
            val = A{1, 1, 100};
            testCase.verifyTrue(abs(val) < testCase.tol, 'Absent harmonic should return 0');
        end

        function testIndexingBracesRead(testCase)
            V = randn(3, 3, 11);
            A = PhasorArray(V);
            % Brace access uses harmonic order: k=0 -> DC -> slice index h+1=6
            dc = A{1, 1, 0};
            testCase.verifyTrue(abs(dc - V(1, 1, 6)) < testCase.tol, 'Brace DC mismatch');
            % k=1 -> slice index h+2=7
            h1 = A{1, 1, 1};
            testCase.verifyTrue(abs(h1 - V(1, 1, 7)) < testCase.tol, 'Brace h=1 mismatch');
        end

        function testIndexingBracesWrite(testCase)
            A = PhasorArray.zeros(2, 2);   % h=0
            A{1, 1, 2} = 5;               % assign to harmonic k=2 (auto-expands)
            testCase.verifyTrue(A.h >= 2, 'h should have grown');
            testCase.verifyTrue(abs(A{1, 1, 2} - 5) < testCase.tol, 'Assignment mismatch');
        end



        function testInvOfConstantMatrix(testCase)
            A = PhasorArray([2 1; 1 3]); % Constant LTI matrix
            B = inv(A);
            testCase.verifyTrue(isa(B, 'PhasorArray'), 'inv should return PhasorArray');
            
            I_check = A * B;
            valI = squeeze(I_check.value);
            err = max(abs(valI - eye(2)), [], 'all');
            testCase.verifyTrue(err < 1e-6, sprintf('A * inv(A) != I (err = %e)', err));
        end

        function testIsreal(testCase)
            v = randn(2, 2, 4);
            A = PhasorArray(v, 'isreal', true);
            testCase.verifyTrue(isreal(A), 'Should be real');
        end

        function testLeftDivideVsInvMulRect(testCase)
            [A, BmatLeft, BvecLeft] = make_rect_case_4x3_h15();

            % A rectangle (n x m), b vector (n x 1)
            check_left_divide_equivalence(A, BvecLeft, 1e-2, 1e-2, 'rectangular: b vector (n x 1)', false);
            % A rectangle (n x m), b rectangle (n x l)
            check_left_divide_equivalence(A, BmatLeft, 1e-2, 1e-2, 'rectangular: b rectangle (n x l)', false);
        end

        function testLeftDivideVsInvMulSquare(testCase)
            [A, Bsquare, BvecLeft] = make_pd_case_3x3_h15();
            Brect = PhasorArray.random(3, 5, 15); % n x m

            % A square (n x n), b vector (n x 1)
            check_left_divide_equivalence(A, BvecLeft, 3e-3, 3e-3, 'square: b vector (n x 1)', true);
            % A square (n x n), b square (n x n)
            check_left_divide_equivalence(A, Bsquare, 3e-3, 3e-3, 'square: b square (n x n)', true);
            % A square (n x n), b rectangle (n x m)
            check_left_divide_equivalence(A, Brect, 3e-3, 3e-3, 'square: b rectangle (n x m)', true);
        end

        function testMldivide(testCase)
            A = PhasorArray([4 1; 1 3]);
            B = PhasorArray.random(2, 1, 1);
            % A \ B
            X = A \ B;
            % Verify A * X == B
            err = max(abs(value((A * X) - B)), [], 'all');
            testCase.verifyTrue(err < 1e-6, sprintf('A \\ B failed verification (err = %e)', err));
        end


        function testNdsdpvarAdd(testCase)
            P = PhasorArray.ndsdpvar(2, 2, 1, "symmetry", "real");
            Q = PhasorArray.ndsdpvar(2, 2, 1, "symmetry", "real");
            R = P + Q;
            testCase.verifyTrue(isa(R, 'PhasorArray'), 'Sum should be PhasorArray');
            testCase.verifyTrue(isa(R.value, 'ndsdpvar') || isa(R.value, 'sdpvar'), 'Sum value should be symbolic');
        end

        function testNdsdpvarComplexCombined(testCase)
            % A deeply combined test challenging padding, alignment, and Cauchy product
            % simultaneously on mixed sdpvar/double arrays.
            has_sdpvar = exist('sdpvar', 'file') == 2;
            testCase.assumeTrue(has_sdpvar, 'YALMIP is required for sdpvar tests');
    
            % Dimensions and varying harmonics
            n = 2; m = 1;
            P = PhasorArray.ndsdpvar(n, n, 2); % h=2
            A = PhasorArray.random(n, n, 1); % h=1
            B = PhasorArray.random(n, m, 3); % h=3
            K = PhasorArray.ndsdpvar(m, n, 1, "symmetry", "real"); % h=1
            Q = PhasorArray.random(n, n, 0); % h=0
    
            % Operation: LMI_like = A'*P + P*A + B*K + K'*B' + Q
            % 1. A'*P -> h = 1 + 2 = 3
            ATP = A' * P;
            % 2. P*A -> h = 2 + 1 = 3
            PA = P * A;
            % 3. B*K -> h = 3 + 1 = 4
            BK = B * K;
            % 4. K'*B' -> h = 1 + 3 = 4
            KTBT = K' * B';
    
            % 5. Sum them all + Q(h=0). The maximum h in the sum is 4.
            % This forces heavy zero-padding on ATP (3->4), PA (3->4), and Q (0->4).
            LMI_term = ATP + PA + BK + KTBT + Q;
    
            % Verifications
            testCase.verifyTrue(isa(LMI_term.value, 'ndsdpvar') || isa(LMI_term.value, 'sdpvar'), ...
                'Complex combined LMI term must strictly preserve the sdpvar graph');
            testCase.verifyTrue(LMI_term.h == 4, 'Combined term should have padded up to max h=4');
            testCase.verifyTrue(size(LMI_term, 1) == n && size(LMI_term, 2) == n, 'Combined term should be nxn');
            testCase.verifyTrue(isreal(LMI_term), 'The symmetric construction should remain real in time domain');
        end

        function testNdsdpvarConstruct(testCase)
            P = PhasorArray.ndsdpvar(3, 3, 2);
            testCase.verifyTrue(isa(P, 'PhasorArray'), 'Should be PhasorArray');
            testCase.verifyTrue(P.h == 2, 'Harmonic order should be 2');
            testCase.verifyTrue(size(P, 1) == 3 && size(P, 2) == 3, 'Size should be 3x3');
        end

        function testNdsdpvarDiag(testCase)
            P = PhasorArray.ndsdpvar(3, 3, 1);
            d = diag(P);
            testCase.verifyTrue(size(d, 1) == 3 && size(d, 2) == 1, 'diag should extract 3x1 vector');
            % blkdiag with ndsdpvar
            Q = PhasorArray.ndsdpvar(2, 2, 0, "symmetry", "real");
            R = blkdiag(P, Q);
            testCase.verifyTrue(size(R, 1) == 5 && size(R, 2) == 5, 'blkdiag should give 5x5');
        end

        function testNdsdpvarMixedType(testCase)
            has_sdpvar = exist('sdpvar', 'file') == 2;
            testCase.assumeTrue(has_sdpvar, 'YALMIP is required for sdpvar tests');
    
            % Test varying harmonic horizons to trigger padding/alignment logic
            % (Crucial: padding with zeros can silently cast sdpvar to double if not careful)
            for h_pair = [1, 3; 3, 1]'
                hP = h_pair(1);
                hD = h_pair(2);
        
                P = PhasorArray.ndsdpvar(2, 2, hP, "symmetry", "real");
                D = PhasorArray.random(2, 2, hD);
        
                % Addition (triggers padding to max(hP, hD))
                S1 = P + D;
                testCase.verifyTrue(isa(S1.value, 'ndsdpvar') || isa(S1.value, 'sdpvar'), sprintf('P(%d)+D(%d) should preserve sdpvar', hP, hD));
                testCase.verifyTrue(S1.h == max(hP, hD), 'Addition should pad to max h');
        
                S2 = D + P;
                testCase.verifyTrue(isa(S2.value, 'ndsdpvar') || isa(S2.value, 'sdpvar'), sprintf('D(%d)+P(%d) should preserve sdpvar', hD, hP));
        
                % Multiplication (triggers Cauchy product, h = hP + hD)
                M1 = P * D;
                testCase.verifyTrue(isa(M1.value, 'ndsdpvar') || isa(M1.value, 'sdpvar'), sprintf('P(%d)*D(%d) should preserve sdpvar', hP, hD));
                testCase.verifyTrue(M1.h == hP + hD, 'Multiplication should sum h');
        
                M2 = D * P;
                testCase.verifyTrue(isa(M2.value, 'ndsdpvar') || isa(M2.value, 'sdpvar'), sprintf('D(%d)*P(%d) should preserve sdpvar', hD, hP));
        
                % Concatenation
                C = [P, D];
                testCase.verifyTrue(isa(C.value, 'ndsdpvar') || isa(C.value, 'sdpvar'), sprintf('[P(%d), D(%d)] should preserve sdpvar', hP, hD));
                testCase.verifyTrue(C.h == max(hP, hD), 'Concatenation should pad to max h');
        
                % Kronecker (h = hP + hD)
                K = kron(P, D);
                testCase.verifyTrue(isa(K.value, 'ndsdpvar') || isa(K.value, 'sdpvar'), sprintf('kron(P(%d), D(%d)) should preserve sdpvar', hP, hD));
                testCase.verifyTrue(K.h == hP + hD, 'Kronecker should sum h');
            end
        end

        function testProjectionsKeepTheSdpvarPayload(testCase)
            testCase.assumeTrue(exist('sdpvar', 'file') == 2, 'YALMIP required');
            % A hermitian variable is complex in general, so all three keep a
            % live expression. Under "real" the imaginary part is structurally
            % zero and YALMIP collapses it to a numeric zero, which is the right
            % answer rather than a lost payload.
            P = PhasorArray.ndsdpvar(2, 2, 1, "symmetry", "hermitian");
            for op = {@mreal, @mimag, @mconj}
                out = op{1}(P);
                testCase.verifyTrue(isa(out.value, 'ndsdpvar') || isa(out.value, 'sdpvar'), ...
                    sprintf('%s dropped the sdpvar payload', func2str(op{1})));
            end

            R = PhasorArray.ndsdpvar(2, 2, 1, "symmetry", "real");
            testCase.verifyEqual(max(abs(mimag(R).value), [], 'all'), 0, ...
                'the imaginary part of a real signal must vanish');
            testCase.verifyTrue(isa(mreal(R).value, 'ndsdpvar') || isa(mreal(R).value, 'sdpvar'), ...
                'mreal dropped the sdpvar payload');
        end

        function testProperties(testCase)
            A = PhasorArray.random(3, 4, 5);
            testCase.verifyTrue(size(A, 1) == 3, 'size(A,1)');
            testCase.verifyTrue(size(A, 2) == 4, 'size(A,2)');
            testCase.verifyTrue(A.h == 5, 'h');
            testCase.verifyTrue(dim(A) == 12, 'dim = n1*n2');
        end

        function testRealnessPredicates(testCase)
            R = PhasorArray.random(2, 2, 2);                       % real in time
            I = 1i * R;                                            % purely imaginary
            C = PhasorArray(complex(randn(2,2,5), randn(2,2,5)));  % neither

            testCase.verifyTrue(isrealp(R),  'a real signal must be real');
            testCase.verifyFalse(isimag(R),  'a non-zero real signal is not imaginary');
            testCase.verifyTrue(isimag(I),   'i*real must be imaginary');
            testCase.verifyFalse(isrealp(I), 'i*real is not real');
            testCase.verifyTrue(iscomplex(C), 'a generic array is complex');
        end

        function testRightDivideVsMulInvRect(testCase)
            [A, ~, ~, BmatRight, BvecRight] = make_rect_case_4x3_h15();

            % xA = b with A rectangle (n x m), b vector (1 x m)
            check_right_divide_equivalence(BvecRight, A, 1e-2, 1e-2, 'rectangular: b vector (1 x m)', false);
            % xA = b with A rectangle (n x m), b rectangle (l x m)
            check_right_divide_equivalence(BmatRight, A, 1e-2, 1e-2, 'rectangular: b rectangle (l x m)', false);
        end

        function testRightDivideVsMulInvSquare(testCase)
            [A, Bsquare, ~, BvecRight] = make_pd_case_3x3_h15();
            Brect = PhasorArray.random(5, 3, 15); % l x n

            % xA = b with A square (n x n), b vector (1 x n)
            check_right_divide_equivalence(BvecRight, A, 3e-3, 3e-3, 'square: b vector (1 x n)', true);
            % xA = b with A square (n x n), b square (n x n)
            check_right_divide_equivalence(Bsquare, A, 3e-3, 3e-3, 'square: b square (n x n)', true);
            % xA = b with A square (n x n), b rectangle (l x n)
            check_right_divide_equivalence(Brect, A, 3e-3, 3e-3, 'square: b rectangle (l x n)', true);
        end

        function testScalarMul(testCase)
            A = PhasorArray.random(2, 2, 3);
            C = 3 * A;
            testCase.verifyTrue(max(abs(C{:, :,0} - 3 * A{:,:,0}), [], 'all') < testCase.tol, 'Scalar multiply DC');
        end

        function testSdpvarPayloadSurvivesArithmetic(testCase)
            testCase.needYalmip();
            P = PhasorArray.ndsdpvar(2, 2, 1, "symmetry", "real");
            D = PhasorArray.random(2, 2, 1);

            testCase.verifySdp(P - D, 'P-D');
            testCase.verifySdp(D - P, 'D-P');
            testCase.verifySdp(-P,    '-P');
            testCase.verifySdp(2 * P, '2*P');
        end

        function testSdpvarPayloadSurvivesRestructuring(testCase)
            testCase.needYalmip();
            P = PhasorArray.ndsdpvar(2, 2, 1, "symmetry", "real");
            D = PhasorArray.random(2, 2, 3);

            testCase.verifySdp(P.',      'transpose');
            testCase.verifySdp(P',       'ctranspose');
            testCase.verifySdp([P; D],   'vertcat with a double');
            testCase.verifySdp(P{1,1},   'brace indexing');
            testCase.verifySdp(trace(P), 'trace');
        end

        function testSdpvarRoundTripThroughSdpval(testCase)
            testCase.needYalmip();
            % After solving, sdpval must hand back a numeric PhasorArray of the
            % same shape, otherwise a gain extracted from P is unusable.
            P = PhasorArray.ndsdpvar(2, 2, 2, "symmetry", "real");
            assign(P.value, ones(2, 2, 5));
            Pv = sdpval(P);
            testCase.verifyClass(Pv, 'PhasorArray');
            testCase.verifyTrue(isnumeric(Pv.value), 'sdpval should return a numeric payload');
            testCase.verifyEqual(Pv.h, P.h, 'sdpval changed the harmonic order');
        end

        function testSdpvarSurvivesToeplitzOperators(testCase)
            testCase.needYalmip();
            % T_tb and F_tb feed every LMI in the toolbox: they must return
            % optimisation variables, not their current numeric value.
            P = PhasorArray.ndsdpvar(2, 2, 2, "symmetry", "real");
            testCase.verifyTrue(isa(P.T_tb(4), 'sdpvar'), 'T_tb collapsed to a constant');
            testCase.verifyTrue(isa(P.F_tb(4), 'sdpvar'), 'F_tb collapsed to a constant');
        end

        function testSubsasgnMismatchH(testCase)
            A = PhasorArray.random(2, 2, 1); % h=1
            B = PhasorArray.random(1, 2, 3); % h=3
            % Assign B to first row of A using braces
            A{1, :} = B;
            testCase.verifyTrue(A.h == 3, 'Assignment should pad to max harmonic h=3');
            testCase.verifyTrue(size(A, 1) == 2 && size(A, 2) == 2, 'Size should remain 2x2');
            
            % Verify that row 2 (which had h=1) was padded with zeros correctly.
            % h=3 gives 7 slices, k=-3..3, so the DC sits at index 4. The
            % harmonics outside the original h=1 are k=+-2,+-3: indices 1,2,6,7.
            valA = squeeze(A.value);
            row2_high_harmonics = valA(2, :, [1 2 6 7]);
            testCase.verifyTrue(max(abs(row2_high_harmonics(:))) < testCase.tol, 'Padding for assigned rows failed');
        end

        function testSubsrefDC(testCase)
            A = PhasorArray.random(2, 2, 2);
            % Extract DC component (h=0) by evaluating at t=0 of the DC part?
            % Actually, let's just test value extraction
            val = A.value;
            testCase.verifyTrue(isnumeric(val), 'value extraction should return numeric');
        end

        function testSubsrefParentheses(testCase)
            A = PhasorArray.random(3, 3, 2);
            % Extract a 2x2 submatrix using braces (which returns a PhasorArray)
            subA = A{1:2, 2:3};
            testCase.verifyTrue(isa(subA, 'PhasorArray'), 'Extraction should return PhasorArray');
            testCase.verifyTrue(size(subA, 1) == 2 && size(subA, 2) == 2, 'Extracted size is wrong');
            testCase.verifyTrue(subA.h == 2, 'Harmonic order should be preserved');
            
            % Check value preservation
            valRef = squeeze(A.value);
            valSub = squeeze(subA.value);
            testCase.verifyTrue(max(abs(valSub(1, 1, :) - valRef(1, 2, :))) < testCase.tol, 'Values corrupted during subsref');
        end

        function testSymArithmetic(testCase)
            A = PhasorArray.sym(2, 2, 1, "A","isreal",true);
            B = PhasorArray.sym(2, 2, 1, "B","isreal",true);
            C = A + B;   % should work symbolically
            D = A * B;   % Cauchy product, h should grow
            testCase.verifyTrue(isa(C, 'PhasorArray'), 'Sum should be PhasorArray');
            testCase.verifyTrue(isa(D, 'PhasorArray'), 'Product should be PhasorArray');
            testCase.verifyTrue(D.h == 2, 'Product h should be 1+1=2');
        end

        function testSymConstruct(testCase)
            A = PhasorArray.sym(2, 2, 3, "M","isreal",true);
            testCase.verifyTrue(isa(A, 'PhasorArray'), 'Should be PhasorArray');
            testCase.verifyTrue(A.h == 3, 'Harmonic order should be 3');
            testCase.verifyTrue(size(A, 1) == 2 && size(A, 2) == 2, 'Size should be 2x2');
        end

        function testSymDiag(testCase)
            v = PhasorArray.sym(3, 1, 2, "v","isreal",true);
            D = diag(v);
            testCase.verifyTrue(size(D, 1) == 3 && size(D, 2) == 3, 'diag(3x1) should be 3x3');
        end

        function testSymPayloadSurvivesArithmetic(testCase)
            testCase.needSymbolic();
            A = PhasorArray.sym(2, 2, 1, "A", "isreal", true);
            B = PhasorArray.sym(2, 2, 1, "B", "isreal", true);

            testCase.verifySym(A + B, 'A+B');
            testCase.verifySym(A - B, 'A-B');
            testCase.verifySym(A * B, 'A*B');
            testCase.verifySym(-A,    '-A');
            testCase.verifySym(3 * A, '3*A');
        end

        function testSymPayloadSurvivesMixedWithDouble(testCase)
            testCase.needSymbolic();
            % Zero padding to a common h is where a cast to double hides.
            A = PhasorArray.sym(2, 2, 1, "A", "isreal", true);
            D = PhasorArray.random(2, 2, 3);

            testCase.verifySym(A + D, 'sym(h=1) + double(h=3)');
            testCase.verifySym(D + A, 'double(h=3) + sym(h=1)');
            testCase.verifySym(A * D, 'sym * double');
            testCase.verifySym(D * A, 'double * sym');

            padded = A + D;
            product = A * D;
            testCase.verifyEqual(padded.h, 3, 'padding to max h failed');
            testCase.verifyEqual(product.h, 4, 'Cauchy product must add the orders');
        end

        function testSymPayloadSurvivesRestructuring(testCase)
            testCase.needSymbolic();
            A = PhasorArray.sym(2, 2, 1, "A", "isreal", true);
            B = PhasorArray.sym(2, 2, 1, "B", "isreal", true);

            testCase.verifySym(A',            'ctranspose');
            testCase.verifySym(A.',           'transpose');
            testCase.verifySym([A, B],        'horzcat');
            testCase.verifySym([A; B],        'vertcat');
            testCase.verifySym(kron(A, B),    'kron');
            testCase.verifySym(blkdiag(A, B), 'blkdiag');
            testCase.verifySym(A{1,1},        'brace indexing');
        end

        function testTranspose(testCase)
            A = PhasorArray.random(3, 2, 5);
            At = A.';
            testCase.verifyTrue(size(At, 1) == 2 && size(At, 2) == 3, 'Transposed size');
            % Check DC transpose
            testCase.verifyTrue(max(abs(At{:, :,0} - A{:,:,0}.'), [], 'all') < testCase.tol, 'DC transpose');
        end

        function testTransposes(testCase)
            A = PhasorArray.random(2, 3, 2);
            % Complex conjugate transpose
            AH = A';
            testCase.verifyTrue(size(AH, 1) == 3 && size(AH, 2) == 2, 'ctranspose size wrong');
            
            % Non-conjugate transpose
            AT = A.';
            testCase.verifyTrue(size(AT, 1) == 3 && size(AT, 2) == 2, 'transpose size wrong');
        end

        function testVertcat(testCase)
            A = PhasorArray.random(2, 3, 3);
            B = PhasorArray.random(4, 3, 5);
            C = [A; B];
            testCase.verifyTrue(size(C, 1) == 6, 'Rows = 2+4');
            testCase.verifyTrue(size(C, 2) == 3, 'Cols');
            testCase.verifyTrue(C.h == 5, 'h = max(3,5)');
        end


        function testReduceRefusesASymbolicPayload(testCase)
            % Reducing compares magnitudes, which an optimisation variable has
            % none of. It used to surface a YALMIP strict-inequality error or a
            % MATLAB dimension error, neither naming the real problem.
            testCase.needYalmip();
            P = PhasorArray.ndsdpvar(2, 2, 2, "symmetry", "real");
            testCase.verifyError(@() reduce(P), 'PhasorArray:reduce:symbolicPayload');
            testCase.verifyError(@() reduce(P, 1), 'PhasorArray:reduce:symbolicPayload');
        end
    end

    methods (Test, TestTags = {'Install'})
        % Smoke set for a fresh install: no optional toolbox, a few seconds,
        % one check per layer. Run with run_all_tests("install").
        function testConstructor3d(testCase)
            V = randn(3, 3, 11);   % h = 5
            A = PhasorArray(V);
            testCase.verifySize(A, [3, 3, 11], 'Size mismatch');
            testCase.verifyTrue(A.h == 5, 'Harmonic order mismatch');
            Aval = A.value;
            testCase.verifyTrue(max(abs(Aval(:) - V(:))) < testCase.tol, 'Value mismatch');
        end

        function testAdd(testCase)
            A = PhasorArray.random(2, 2, 3);
            B = PhasorArray.random(2, 2, 3);
            C = A + B;
            % Check DC: C{:,:,0} = A{:,:,0} + B{:,:,0}
            testCase.verifyTrue(max(abs(C{:, :,0} - A{:,:,0} - B{:,:,0}), [], 'all') < testCase.tol, 'DC addition');
        end

        function testMul(testCase)
            A = PhasorArray.random(2, 2, 2);  % h=2
            B = PhasorArray.random(2, 2, 3);  % h=3
            C = A * B;
            testCase.verifyTrue(C.h == 5, 'h_C should be h_A + h_B = 5');
            % Verify via time-domain sampling
            t = linspace(0, 2*pi, 200);
            for k = 1:5
                tk = t(k*30);
                At = evalp(A, tk);
                Bt = evalp(B, tk);
                Ct = evalp(C, tk);
                err = max(abs(At * Bt - Ct), [], 'all');
                testCase.verifyTrue(err < 1e-6, sprintf('Time-domain verification failed at t=%g', tk));
            end
        end

        function testIndexingParens(testCase)
            V = randn(3, 3, 11);
            A = PhasorArray(V);
            % Parentheses access the raw 3D array
            slice = A(1, 1, 6);   % (1,1, h+1) = DC for h=5
            testCase.verifyTrue(abs(slice - V(1, 1, 6)) < testCase.tol, 'Parentheses indexing mismatch');
        end

        function testInv(testCase)
            A = PhasorArray.eye(2) + 0.1 * PhasorArray.random(2, 2, 2);
            Ainv = inv(A);
            I_approx = A * Ainv;  % should be ≈ identity
            % Check DC of result is close to I
            I_dc = I_approx{:,:,0};
            testCase.verifyTrue(max(abs(I_dc - eye(2)), [], 'all') < 0.05, 'A*inv(A) DC should be ~I');
        end
    end

end

function [A, Bmat, BvecLeft, BvecRight] = make_pd_case_3x3_h15()
% Build a 3x3 periodic SPD-like matrix with h=15 by rejection on det(A(t)).
    n = 3;
    h = 15;
    tGrid = linspace(0, 2*pi, 61);
    maxTries = 30;

    for k = 1:maxTries
        R = PhasorArray.random(n, n, h);
        S = 0.2 * (0.5 / n) * (R + R');          % Hermitian periodic perturbation
        A = PhasorArray.eye(n) + S;
        Bmat = PhasorArray.random(n, n, h);
        BvecLeft = PhasorArray.random(n, 1, h);
        BvecRight = PhasorArray.random(1, n, h);

        detA = det(A);
        detVals = zeros(size(tGrid));
        minSv = inf;
        for i = 1:numel(tGrid)
            detVals(i) = real(evalp(detA, tGrid(i)));
            s = svd(evalp(A, tGrid(i)));
            minSv = min(minSv, min(s));
        end

        if all(isfinite(detVals)) && all(detVals > 1e-2) && minSv > 0.2
            return
        end
    end

    error('Could not generate a valid 3x3 periodic matrix A with strictly positive det(A(t)).');
end

function [A, BmatLeft, BvecLeft, BmatRight, BvecRight] = make_rect_case_4x3_h15()
% Build a well-conditioned periodic rectangular matrix A (4x3), h=15.
    nrow = 4;
    ncol = 3;
    h = 15;
    tGrid = linspace(0, 2*pi, 61);
    maxTries = 120;

    A0 = [eye(ncol); zeros(nrow - ncol, ncol)];  % full column-rank baseline

    for k = 1:maxTries
        R = PhasorArray.random(nrow, ncol, h);
        A = PhasorArray(A0) + 0.02 * R;

        minSv = inf;
        for i = 1:numel(tGrid)
            s = svd(evalp(A, tGrid(i)));
            minSv = min(minSv, min(s));
        end

        if isfinite(minSv) && minSv > 0.2
            BmatLeft = PhasorArray.random(nrow, 5, h);     % for A\B (n x l)
            BvecLeft = PhasorArray.random(nrow, 1, h);     % for A\B
            BmatRight = PhasorArray.random(5, ncol, h);    % for B/A (l x m)
            BvecRight = PhasorArray.random(1, ncol, h);    % for B/A
            return
        end
    end

    error('Could not generate a well-conditioned rectangular periodic matrix A (4x3).');
end

function check_left_divide_equivalence(A, B, tolCons, tolRes, label, enforceSmallResidual)

    cleanupWarnings = suppress_expected_divide_warnings(~enforceSmallResidual);
    Xdiv = mldivide(A, B, "autoUpdateh", true, "verbose", 0);

    clear cleanupWarnings;
    Xinv = inv(A) * B;

    assert(size(Xdiv,1) == size(Xinv,1) && size(Xdiv,2) == size(Xinv,2), ...
        'Size mismatch in left-division comparison (%s).', label);

    eDiv = energy(A * Xdiv - B);
    eInv = energy(A * Xinv - B);
    relEnergy = abs(eDiv - eInv) / max([eDiv, eInv, 1]);
    assert(relEnergy < 1e-6, ...
        'Residual energy mismatch (%s): rel=%e, Eslash=%e, Einv=%e', ...
        label, relEnergy, eDiv, eInv);
end

function check_right_divide_equivalence(B, A, tolCons, tolRes, label, enforceSmallResidual)
    cleanupWarnings = suppress_expected_divide_warnings(~enforceSmallResidual);
    Xdiv = mrdivide(B, A, "autoUpdateh", true, "verbose", 0);
    clear cleanupWarnings;
    Xinv = B * inv(A);

    assert(size(Xdiv,1) == size(Xinv,1) && size(Xdiv,2) == size(Xinv,2), ...
        'Size mismatch in right-division comparison (%s).', label);

    eDiv = energy(Xdiv * A - B);
    eInv = energy(Xinv * A - B);
    relEnergy = abs(eDiv - eInv) / max([eDiv, eInv, 1]);
    assert(relEnergy < 1e-6, ...
        'Residual energy mismatch (%s): rel=%e, Eslash=%e, Einv=%e', ...
        label, relEnergy, eDiv, eInv);
end

function cleanupObj = suppress_expected_divide_warnings(shouldSuppress)
    if ~shouldSuppress
        cleanupObj = onCleanup(@()[]);
        return;
    end

    warningIds = {
        'phasorArray:mlHmcDivide:stagnation'
        'phasorArray:mlHmcDivide:residual'
        'phasorArray:mlHmcDivide:maxhReached'
        'phasorArray:mlHmcDivide:slowConvergence'
        'phasorArray:mlHmcDivide:unreachable'
    };

    prevStates = cell(size(warningIds));
    for k = 1:numel(warningIds)
        prevStates{k} = warning('query', warningIds{k});
        warning('off', warningIds{k});
    end

    cleanupObj = onCleanup(@() restore_warning_states(prevStates));
end




function restore_warning_states(prevStates)
    for k = 1:numel(prevStates)
        warning(prevStates{k}.state, prevStates{k}.identifier);
    end
end

function e = time_domain_consistency(X1, X2)
% Compare two PhasorArray objects in time domain when harmonic bases differ.
    tGrid = linspace(0, 2*pi, 41);
    e = 0;
    for i = 1:numel(tGrid)
        D = evalp(X1, tGrid(i)) - evalp(X2, tGrid(i));
        e = max(e, norm(D, 'fro'));
    end
end

function e = time_domain_residual_left(A, X, B)
    tGrid = linspace(0, 2*pi, 41);
    e = 0;
    for i = 1:numel(tGrid)
        D = evalp(A, tGrid(i)) * evalp(X, tGrid(i)) - evalp(B, tGrid(i));
        e = max(e, norm(D, 'fro'));
    end
end

function e = time_domain_residual_right(X, A, B)
    tGrid = linspace(0, 2*pi, 41);
    e = 0;
    for i = 1:numel(tGrid)
        D = evalp(X, tGrid(i)) * evalp(A, tGrid(i)) - evalp(B, tGrid(i));
        e = max(e, norm(D, 'fro'));
    end

end
