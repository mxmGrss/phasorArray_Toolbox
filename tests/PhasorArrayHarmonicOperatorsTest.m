classdef PhasorArrayHarmonicOperatorsTest < matlab.unittest.TestCase
    %PHASORARRAYHARMONICOPERATORSTEST  Toeplitz and block-Toeplitz operators, and their inverses.
    %
    %   T_tb, T_bt, BT, TB, the Hankel builders and the N and S operators, each against
%   its extractor and against the product it is meant to implement. Shapes are
%   deliberately non-square: a row/column swap is invisible on an n x n operand.

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
        function v = payload(~, x)
            if isa(x, 'PhasorArray'), v = x.value; else, v = x; end
        end

    end

    methods (Test)
        function testBlockToeplitzActionMatchesTheProduct(testCase)
            A = PhasorArray.random(2, 3, 2);
            x = PhasorArray.random(3, 1, 2);
            h = 8;
            testCase.verifyEqual(A.BT(h) * x.F_bt(h), F_bt(A * x, h), 'AbsTol', 1e-9, ...
                'BT does not implement the product on a rectangular operand');
        end

        function testBt2Tb(testCase)
            % BT and TB are the same operator in two block layouts: BT_2_TB is the
            % perfect-shuffle change of layout and must map one onto the other.
            A = PhasorArray.random(2, 2, 2);
            h = 3;
            err = norm(BT_2_TB(A.BT(h), 2*h+1) - A.TB(h), inf);
            testCase.verifyTrue(err < testCase.tol, 'BT_2_TB(BT(A)) ~= TB(A)');
        end

        function testBt2arrayRoundtrip(testCase)
            % Extracting the harmonic slices back from BT(A) must recover A.
            A = PhasorArray.random(2, 2, 2);
            h = A.h;
            pA = BT2array(A.BT(h), 2*h+1);
            % BT(A,h) holds harmonics up to 2h; original A occupies the center
            hBig = 2*h;
            center = pA(:, :, hBig+1-h : hBig+1+h);
            ref = A.value;
            testCase.verifyTrue(norm(center(:) - ref(:), inf) < testCase.tol, 'BT2array(BT(A)) ~= A');
        end

        function testBtExtractionGivesTheArray(testCase)
            A = PhasorArray.random(3, 3, 2);
            h = A.h;
            pA = BT2array(A.BT(h), 2*h + 1);
            hBig = (size(pA, 3) - 1) / 2;
            centre = pA(:, :, hBig+1-h : hBig+1+h);
            testCase.verifyEqual(centre, A.value, 'AbsTol', testCase.tol, ...
                'BT2array(BT(A)) does not recover A');
        end

        function testBtFbtMul(testCase)
            % Documented key property of F_bt: F_bt(A*x) = BT(A,h) * F_bt(x,h)
            % (truncation-exact when h >= hA + hx).
            A = PhasorArray.random(2, 2, 2);
            x = PhasorArray.random(2, 1, 2);
            h = 4;
            lhs = A.BT(h) * x.F_bt(h);
            rhs = F_bt(A * x, h);
            testCase.verifyTrue(norm(lhs - rhs, inf) < testCase.tol, 'BT * F_bt ~= F_bt(A*x)');
        end

        function testBtNbtSpectrum(testCase)
            % Harmonic dynamics in BT layout: eig(BT(A) - N_bt) must contain the
            % prescribed Floquet exponents (omega = 1, T = 2*pi).
            mus = [-0.7; -0.3];
            A = PhasorArray.randomWithNPole(mus, 6, seed=1);
            H = 12;
            M = A.BT(H) - N_bt(2, H, 2*pi);
            e = eig(full(M));
            err = max(arrayfun(@(m) min(abs(e - m)), mus));
            testCase.verifyTrue(err < 1e-8, 'BT - N_bt spectrum misses Floquet exponents');
        end

        function testBtScalarEqTb(testCase)
            % On a scalar PhasorArray the block layouts are trivial: BT and TB must
            % return the exact same Toeplitz matrix.
            a = PhasorArray.random(1, 1, 3);
            h = 5;
            testCase.verifyTrue(norm(a.BT(h) - a.TB(h), inf) < testCase.tol, 'scalar BT ~= TB');
        end

        function testBthankelScalar(testCase)
            % Scalar case: BT and TB layouts coincide, Hankel matrices too.
            a = PhasorArray.random(1, 1, 3);
            [HpJ1, JHm1, Hp1, Hm1] = a.TBHankel(4);
            [HpJ2, JHm2, Hp2, Hm2] = a.BTHankel(4);
            err = max([norm(HpJ1-HpJ2,inf), norm(JHm1-JHm2,inf), ...
                       norm(Hp1-Hp2,inf),   norm(Hm1-Hm2,inf)]);
            testCase.verifyTrue(err < testCase.tol, 'scalar BTHankel ~= TBHankel');
        end

        function testFtbReadBackGivesTheArray(testCase)
            for shape = {[2 3], [3 1], [1 4], [3 3]}
                n = shape{1}(1); m = shape{1}(2);
                A = PhasorArray.random(n, m, 2);
                back = testCase.payload(F_tb_2_PhasorArray(A.F_tb(A.h), n, m));
                testCase.verifyEqual(back, A.value, 'AbsTol', testCase.tol, ...
                    sprintf('F_tb round trip failed on %dx%d', n, m));
            end
        end

        function testLayoutConversionsAreInverse(testCase)
            % BT and TB are the same operator under a perfect shuffle, so the
            % two converters must undo each other in both directions.
            A = PhasorArray.random(2, 2, 2);
            h = 3;
            nb = 2*h + 1;
            BTm = A.BT(h);
            TBm = A.TB(h);

            testCase.verifyEqual(TB_2_BT(BT_2_TB(BTm, nb), nb), BTm, 'AbsTol', testCase.tol, ...
                'TB_2_BT o BT_2_TB is not the identity');
            testCase.verifyEqual(BT_2_TB(TB_2_BT(TBm, nb), nb), TBm, 'AbsTol', testCase.tol, ...
                'BT_2_TB o TB_2_BT is not the identity');
            testCase.verifyEqual(BT_2_TB(BTm, nb), TBm, 'AbsTol', testCase.tol, ...
                'the two layouts do not describe the same operator');
        end

        function testNTb(testCase)
            n = 3;
            h = 5;
            T = 2*pi;
            N = N_tb(n, h, T);
            expectedSize = n * (2*h + 1);
            testCase.verifySize(N, [expectedSize, expectedSize], 'N_tb size');
            % In TB formalism: N = I_n ⊗ diag(jkω)
            % The k=0 block is at row/col indices (h*n+1):(h*n+n) in each n×n sub-block
            % But in TB formalism, the diagonal is organized as n blocks of (2h+1)×(2h+1)
            % The (h+1)-th diagonal entry in each block corresponds to k=0 and should be 0
            % Simpler check: the diagonal of N at positions corresponding to k=0 should be 0
            Ndiag = diag(N);
            % In TB: I_n ⊗ diag(jkω), so diagonal = repmat(jkω vector, n, 1)
            kvals = -h:h;
            omega = 2*pi/T;
            expected_diag = repmat(1i * kvals(:) * omega, n, 1);
            testCase.verifyTrue(max(abs(Ndiag - expected_diag)) < testCase.tol, 'N_tb diagonal mismatch');
        end

        function testNbtHardcoded(testCase)
            % N_bt is the same operator in the Block-Toeplitz layout, so the
            % Kronecker order is reversed: kron(diag(j*k*omega), eye(n)).
            T = 2*pi;

            % For a scalar, N_tb and N_bt are identical
            expected_N_bt_scalar = diag(1i*(-1:1));
            N_bt_scalar = N_bt(1, 1, T);
            testCase.verifyEqual(N_bt_scalar, expected_N_bt_scalar, 'AbsTol', testCase.tol, ...
                'N_bt(Nc=1) does not match hardcoded expectation for scalar.');

            expected_N_bt_matrix = kron(diag(1i*(-1:1)), eye(2));
            N_bt_matrix = N_bt(2, 1, T);
            testCase.verifyEqual(N_bt_matrix, expected_N_bt_matrix, 'AbsTol', testCase.tol, ...
                'N_bt(Nc=1) does not match hardcoded expectation for 2x2 matrix.');
        end


        function testSTb(testCase)
            % S_tb(h, phi, n) builds kron(I_n, diag(exp(jk*phi)))
            h = 5;
            phi = pi/4;
            n = 2;
            S = S_tb(h, phi, n);
            expectedSize = n * (2*h + 1);
            testCase.verifySize(S, [expectedSize, expectedSize], 'S_tb size');
            % Diagonal should be exp(jk*phi) repeated n times
            Sdiag = diag(S);
            kvals = -h:h;
            expected_diag = repmat(exp(1i * kvals(:) * phi), n, 1);
            testCase.verifyTrue(max(abs(Sdiag - expected_diag)) < 1e-12, 'S_tb diagonal entries');
        end

        function testSTbPhaseshift(testCase)
            % Verify: S_tb(h, phi, n) * F_tb(A, h) == F_tb(PhaseShift(A, phi), h)
            % This tests that the Toeplitz dephasing operator S_tb is consistent
            % with the PhasorArray PhaseShift method.
            A = PhasorArray.random(2, 1, 3);  % 2x1 column vector, h=3
            phi = pi/6;
            hTrunc = 10;
    
            % Method 1: PhaseShift then Fourier vector
            A_shifted = A.PhaseShift(phi);
            F_shifted = A_shifted.F_tb(hTrunc);
    
            % Method 2: Fourier vector then multiply by S_tb
            F_A = A.F_tb(hTrunc);
            n = size(A, 1);
            S = S_tb(hTrunc, phi, n);
            F_via_S = S * F_A;
    
            err = max(abs(F_shifted - F_via_S));
            testCase.verifyTrue(err < testCase.tol, sprintf('S_tb * F_tb vs PhaseShift(A,phi).F_tb error: %e', err));
        end

        function testSpbthankel(testCase)
            A = PhasorArray.random(2, 2, 2);
            [HpJ1, JHm1, Hp1, Hm1] = A.BTHankel(4);
            [HpJ2, JHm2, Hp2, Hm2] = A.spBTHankel(4);
            testCase.verifyTrue(issparse(HpJ2), 'spBTHankel not sparse');
            err = max([norm(full(HpJ2)-HpJ1,inf), norm(full(JHm2)-JHm1,inf), ...
                       norm(full(Hp2)-Hp1,inf),   norm(full(Hm2)-Hm1,inf)]);
            testCase.verifyTrue(err < testCase.tol, 'spBTHankel ~= BTHankel');
        end

        function testSptbhankel(testCase)
            % Sparse and dense Hankel builders must agree entrywise.
            A = PhasorArray.random(2, 2, 2);
            [HpJ1, JHm1, Hp1, Hm1] = A.TBHankel(4);
            [HpJ2, JHm2, Hp2, Hm2] = A.spTBHankel(4);
            testCase.verifyTrue(issparse(HpJ2), 'spTBHankel not sparse');
            err = max([norm(full(HpJ2)-HpJ1,inf), norm(full(JHm2)-JHm1,inf), ...
                       norm(full(Hp2)-Hp1,inf),   norm(full(Hm2)-Hm1,inf)]);
            testCase.verifyTrue(err < testCase.tol, 'spTBHankel ~= TBHankel');
        end

        function testTBt(testCase)
            A = PhasorArray.random(2, 2, 3);
            hTrunc = 8;
            T = A.T_bt(hTrunc);
            n = size(A, 1);
            Nt = 2*hTrunc + 1;
            expectedSize = n * Nt;
            testCase.verifySize(T, [expectedSize, expectedSize], 'T_bt size');
        end


        function testTTbRoundtrip(testCase)
            A = PhasorArray.random(2, 2, 3);
            hTrunc = 10;
            T = A.T_tb(hTrunc);
            % Reconstruct: the central column blocks should contain A's phasors
            n = size(A, 1);
            Nt = 2*hTrunc + 1;
            expectedSize = n * Nt;
            testCase.verifySize(T, [expectedSize, expectedSize], 'T_tb size mismatch');
        end

        function testTbExtractionGivesTheArray(testCase)
            % TB(A,h) carries harmonics up to 2h; the original sits in the centre.
            A = PhasorArray.random(2, 2, 2);
            h = A.h;
            pA = TB2array(A.TB(h), size(A, 1));
            hBig = (size(pA, 3) - 1) / 2;
            centre = pA(:, :, hBig+1-h : hBig+1+h);
            testCase.verifyEqual(centre, A.value, 'AbsTol', testCase.tol, ...
                'TB2array(TB(A)) does not recover A');
        end

        function testTbmtimes(testCase)
            % Ground truth: exact product in the harmonic domain, then Toeplitz.
            A = PhasorArray.random(2, 2, 2);
            B = PhasorArray.random(2, 2, 2);
            h = 4;
            Ttrue = TB(A*B, h);
            testCase.verifyTrue(norm(A.TBmtimes(B, h) - Ttrue, inf) < 1e-10, 'TBmtimes ~= TB(A*B)');
            testCase.verifyTrue(norm(full(A.spTBmtimes(B, h)) - Ttrue, inf) < 1e-10, 'spTBmtimes ~= TB(A*B)');
        end

        function testTbtHardcodedMatrix(testCase)
            % Hardcoded test for T_bt (Block Toeplitz) with a 2x2 PhasorArray
            M0 = [1, 2; 3, 4];
            M1 = [5+1i, 6; 7, 8-2i];
            Mm1 = [5-1i, 6; 7, 8+2i];
            
            hmc = cat(3, Mm1, M0, M1);
            A = PhasorArray(hmc);
            
            % Truncation Nc = 1 (size should be 3 blocks of 2x2 = 6x6)
            Z = zeros(2,2);
            expected_T_bt_1 = [
                M0,  Mm1, Z;
                M1,  M0,  Mm1;
                Z,   M1,  M0
            ];
            
            T_bt_1 = T_bt(A, 1);
            testCase.verifyEqual(T_bt_1, expected_T_bt_1, 'AbsTol', testCase.tol, ...
                'T_bt(Nc=1) does not match hardcoded expectation for 2x2 matrix.');
        end

        function testTftbReadBackGivesTheArray(testCase)
            for shape = {[2 3], [4 1], [3 3]}
                n = shape{1}(1); m = shape{1}(2);
                A = PhasorArray.random(n, m, 2);
                back = testCase.payload(TFTB_2_array(array2TFTB(A, A.h), n, m));
                testCase.verifyEqual(back, A.value, 'AbsTol', testCase.tol, ...
                    sprintf('array2TFTB round trip failed on %dx%d', n, m));
            end
        end

        function testToeplitzActionMatchesTheProduct(testCase)
            % T_tb(A) * F_tb(x) = F_tb(A x), on a non-square A.
            A = PhasorArray.random(2, 3, 2);
            x = PhasorArray.random(3, 1, 2);
            h = 8;
            testCase.verifyEqual(A.T_tb(h) * x.F_tb(h), F_tb(A * x, h), 'AbsTol', 1e-9, ...
                'T_tb does not implement the product on a rectangular operand');
        end

        function testTtbHardcodedMatrix(testCase)
            % Hardcoded test for T_tb (Toeplitz Block) with a 2x2 PhasorArray
            M0 = [1, 2; 3, 4];
            M1 = [5+1i, 6; 7, 8-2i];
            Mm1 = [5-1i, 6; 7, 8+2i];
            
            hmc = cat(3, Mm1, M0, M1);
            A = PhasorArray(hmc);
            
            % For T_tb, we build a 2x2 block matrix where each block is 3x3 (for Nc=1)
            % Element (1,1): c0=1, c1=5+1i, cm1=5-1i
            T11 = [1, 5-1i, 0; 5+1i, 1, 5-1i; 0, 5+1i, 1];
            % Element (1,2): c0=2, c1=6, cm1=6
            T12 = [2, 6, 0; 6, 2, 6; 0, 6, 2];
            % Element (2,1): c0=3, c1=7, cm1=7
            T21 = [3, 7, 0; 7, 3, 7; 0, 7, 3];
            % Element (2,2): c0=4, c1=8-2i, cm1=8+2i
            T22 = [4, 8+2i, 0; 8-2i, 4, 8+2i; 0, 8-2i, 4];
            
            expected_T_tb_1 = [
                T11, T12;
                T21, T22
            ];
            
            T_tb_1 = T_tb(A, 1);
            testCase.verifyEqual(T_tb_1, expected_T_tb_1, 'AbsTol', testCase.tol, ...
                'T_tb(Nc=1) does not match hardcoded expectation for 2x2 matrix.');
        end


    end

    methods (Test, TestTags = {'Install'})
        % Smoke set for a fresh install: no optional toolbox, a few seconds,
        % one check per layer. Run with run_all_tests("install").
        function testTtbHardcodedScalar(testCase)
            % Hardcoded test for T_tb (Toeplitz Block) with a scalar PhasorArray
            % c_0 = 2, c_1 = 3 + 1i, c_-1 = 3 - 1i
            hmc = cat(3, 3-1i, 2, 3+1i);
            A = PhasorArray(hmc);
            
            % Truncation Nc = 1
            expected_T_tb_1 = [
                2,      3-1i,   0;
                3+1i,   2,      3-1i;
                0,      3+1i,   2
            ];
            
            T_tb_1 = T_tb(A, 1);
            testCase.verifyEqual(T_tb_1, expected_T_tb_1, 'AbsTol', testCase.tol, ...
                'T_tb(Nc=1) does not match hardcoded expectation for scalar.');
            
            % Truncation Nc = 2
            expected_T_tb_2 = [
                2,      3-1i,   0,      0,      0;
                3+1i,   2,      3-1i,   0,      0;
                0,      3+1i,   2,      3-1i,   0;
                0,      0,      3+1i,   2,      3-1i;
                0,      0,      0,      3+1i,   2
            ];
            
            T_tb_2 = T_tb(A, 2);
            testCase.verifyEqual(T_tb_2, expected_T_tb_2, 'AbsTol', testCase.tol, ...
                'T_tb(Nc=2) does not match hardcoded expectation for scalar.');
        end

        function testTTbMul(testCase)
            % A*B in PhasorArray should be consistent with T_tb(A) * F_tb(B)
            A = PhasorArray.random(2, 2, 2);
            B = PhasorArray.random(2, 1, 2);
            C = A * B;
    
            hTrunc = 10;
            T_A = A.T_tb(hTrunc);
            F_B = B.F_tb(hTrunc);
            F_C_expected = T_A * F_B;
            F_C_actual = C.F_tb(hTrunc);
    
            err = max(abs(F_C_expected - F_C_actual));
            testCase.verifyTrue(err < 1e-8, 'T_tb * F_tb inconsistency');
        end

        function testNtbHardcoded(testCase)
            % N_tb is the harmonic derivative, diag(j*k*omega), laid out
            % Toeplitz-Block: kron(eye(n), diag(j*k*omega)).
            T = 2*pi;            % omega = 1, so the entries are j*k

            % Nc = 2, scalar
            expected_N_tb_scalar = diag(1i*(-2:2));
            N_tb_scalar = N_tb(1, 2, T);
            testCase.verifyEqual(N_tb_scalar, expected_N_tb_scalar, 'AbsTol', testCase.tol, ...
                'N_tb(Nc=2) does not match hardcoded expectation for scalar.');

            % Nc = 1 for 2x2 matrix: harmonics run fastest, blocks outermost
            expected_N_tb_matrix = kron(eye(2), diag(1i*(-1:1)));
            N_tb_matrix = N_tb(2, 1, T);
            testCase.verifyEqual(N_tb_matrix, expected_N_tb_matrix, 'AbsTol', testCase.tol, ...
                'N_tb(Nc=1) does not match hardcoded expectation for matrix.');
        end
    end
end
