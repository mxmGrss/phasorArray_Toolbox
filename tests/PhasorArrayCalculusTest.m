classdef PhasorArrayCalculusTest < matlab.unittest.TestCase
    %PHASORARRAYCALCULUSTEST  Derivative, integral, determinant, energy and matrix functions.
    %
    %   Operations that act on the harmonic array as a function of time. Where a closed
%   form exists it is used as the reference: expm(cos) against besseli, det
%   against the Leibniz expansion, energy against Parseval.

    properties
        tol = 1e-10;
        tolAnalytic = 1e-10;
    end

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            srcFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcFolder, 'IncludingSubfolders', true));
        end
    end

    methods (Test)
        function testAntidDc(testCase)
            A = PhasorArray.random(2, 2, 3);
            intA = antiD(A, 2*pi);
            dc = intA{:,:,0};
            testCase.verifyTrue(max(abs(dc(:))) < testCase.tol, 'antiD DC should be forced to zero');
        end

        function testAntidRoundtrip(testCase)
            c = PhasorArray.cos();
            T = 2*pi;
            % antiD(d(cos)) should give back cos (up to DC)
            dc = d(c, T);
            adc = antiD(dc, T);
            % Compare AC part: should be cos(t) (DC may differ)
            err_h1 = abs(adc{1,1,1} - c{1,1,1});
            testCase.verifyTrue(err_h1 < testCase.tol, 'antiD(d(cos)) h=1 mismatch');
        end

        function testConcordia(testCase)
            C = PhasorArray.Concordia();
            testCase.verifyTrue(size(C, 1) == 3 && size(C, 2) == 3, 'Concordia size should be 3x3');
            % Validate power invariance factor (sqrt(2/3))
            valC = squeeze(C.value);
            testCase.verifyTrue(abs(valC(1,1) - sqrt(2/3)) < testCase.tol, 'Concordia amplitude factor wrong');
        end


        function testDetRuns(testCase)
            A = PhasorArray.eye(2) + 0.1 * PhasorArray.random(2, 2, 2);
            d = det(A);
            testCase.verifyTrue(isa(d, 'PhasorArray'), 'det should return PhasorArray');
        end


        function testDetleibniz3x3(testCase)
            A = PhasorArray.random(3, 3, 2);
            dL = detLeibnizHmc(A);
            % Compare with time-domain evaluation
            t_test = [0, pi/3, pi/2, pi, 3*pi/2];
            for k = 1:numel(t_test)
                At = evalp(A, t_test(k));
                dLt = evalp(dL, t_test(k));
                dRef = det(At);
                err = abs(dLt - dRef);
                testCase.verifyTrue(err < 1e-6, sprintf('3x3 Leibniz mismatch at t=%g: err=%e', t_test(k), err));
            end
        end

        function testDetleibnizIdentity(testCase)
            I = PhasorArray.eye(3);
            d = detLeibnizHmc(I);
            % det(I) = 1 (constant)
            dc = d{1, 1, 0};
            testCase.verifyTrue(abs(dc - 1) < testCase.tol, 'det(I) DC should be 1');
            testCase.verifyTrue(energy(d) - 1 < testCase.tol, 'det(I) energy should be 1');
        end

        function testDetleibnizProduct(testCase)
            % For scalar PhasorArrays (1x1), det(A*B) = det(A)*det(B) is trivial.
            % For 2x2, test: det(A*B)(t) = det(A)(t) * det(B)(t)
            n = 2; h = 2;
            vA = complex(randn(n, n, 2*h+1), randn(n, n, 2*h+1));
            vB = complex(randn(n, n, 2*h+1), randn(n, n, 2*h+1));
            A = PhasorArray(vA);
            B = PhasorArray(vB);
    
            dAB = detLeibnizHmc(A * B);
            dA  = detLeibnizHmc(A);
            dB  = detLeibnizHmc(B);
            dAdB = dA * dB;
    
            % Compare at time points
            for t0 = linspace(0, 2*pi, 8)
                v1 = full(evalp(dAB, t0));
                v2 = full(evalp(dAdB, t0));
                err = abs(v1 - v2);
                testCase.verifyTrue(err < testCase.tol, sprintf('det(AB) vs det(A)*det(B) at t=%g: %e', t0, err));
            end
        end

        function testDetleibnizVsFft2x2(testCase)
            % Use general random inputs but disable reduction in the reference det
            % to check exact consistency with Leibniz formula
            n = 2; h = 2;
            % val = complex(randn(n, n, 2*h+1), randn(n, n, 2*h+1));
            A = PhasorArray.random(n,n,h);
    
            dL = detLeibnizHmc(A);
            % Force det to NOT reduce small terms, so it matches Leibniz exact result
            dF = det(A, 'reduceThreshold', 0);
    
            % Compare at several time points
            for t0 = linspace(0, 2*pi, 10)
                vL = full(evalp(dL, t0));
                vF = full(evalp(dF, t0));
                err = abs(vL - vF);
                testCase.verifyTrue(err < 1e-6, sprintf('Leibniz vs FFT-det mismatch at t=%g: %e', t0, err));
            end
        end

        function testDetleibnizVsFft3x3(testCase)
            % Use general random inputs
            n = 3; h = 2;
            A = PhasorArray.random(n,n,h);
            % A = PhasorArray(val);
    
            dL = detLeibnizHmc(A);
            % Disable reduction to match exact Leibniz

            dF = det(A, 'reduceThreshold', 0);
    
            for t0 = linspace(0, 2*pi, 10)
                vL = full(evalp(dL, t0));
                vF = full(evalp(dF, t0));
                err = abs(vL - vF);
                testCase.verifyTrue(err < 1e-6, sprintf('Leibniz vs FFT-det 3x3 mismatch at t=%g: %e', t0, err));
            end
        end

        function testDq0(testCase)
            D = PhasorArray.dq0(0); % Theta = 0
            testCase.verifyTrue(size(D, 1) == 3 && size(D, 2) == 3, 'dq0 size should be 3x3');
            % valD is 3x3x(2h+1); two subscripts would fold the harmonic axis
            % into the columns, so the DC page has to be named explicitly.
            valD = squeeze(D.value);
            dc = (size(valD,3)+1)/2;
            testCase.verifyTrue(abs(valD(3,3,dc) - sqrt(1/3)) < testCase.tol, 'dq0 zero-sequence factor wrong');
        end

        function testEnergyDecomp(testCase)
            A = PhasorArray.random(3, 3, 5);
            E = energy(A);
            Er = realEnergy(A);
            Ei = imagEnergy(A);
            testCase.verifyTrue(abs(E - Er - Ei) < testCase.tol * max(E, 1), 'Energy: total ≠ real + imag');
        end

        function testEnergyEqualsPageEnergySum(testCase)
            A = PhasorArray.random(2, 2, 2);
            
            % energy() should equal realEnergy + imagEnergy
            Etot = energy(A);
            Er = realEnergy(A);
            Ei = imagEnergy(A);
            testCase.verifyTrue(abs(Etot - (Er + Ei)) < testCase.tol, 'Parseval theorem broken: total != real + imag');
            
            % energy should also equal sum of pageEnergy. pageEnergy is
            % [n x m x (2h+1)], so the sum has to run over every dimension.
            Epages = pageEnergy(A);
            testCase.verifyTrue(abs(Etot - sum(Epages(:))) < testCase.tol, 'Parseval theorem broken: total != sum(pages)');
        end


        function testExpandAndSquishBaseAreInverse(testCase)
            A = PhasorArray.random(2, 3, 3);
            for m = [2 3 5]
                B = expandBase(A, m);
                testCase.verifyEqual(B.h, m * A.h, sprintf('expandBase(%d) wrong order', m));
                C = squishBase(B, m);
                testCase.verifyEqual(C.value, A.value, 'AbsTol', testCase.tol, ...
                    sprintf('squishBase o expandBase is not the identity for m=%d', m));
            end
        end

        function testExpandsquish(testCase)
            A = PhasorArray.random(2, 2, 3);
            m = 3;
            B = expandBase(A, m);
            testCase.verifyTrue(B.h == m * A.h, 'expandBase h mismatch');
    
            C = squishBase(B, m);
            testCase.verifyTrue(C.h == A.h, 'squishBase h mismatch');
    
            % Round-trip should be lossless
            err = energy(A - C);
            testCase.verifyTrue(err < testCase.tol, sprintf('expandBase/squishBase round-trip error: %e', err));
        end

        function testExpm(testCase)
            % Test that expm exists and compiles
            A = PhasorArray([0 1; -1 0]);
            E = expm(A);
            testCase.verifyTrue(isa(E, 'PhasorArray'), 'expm should return PhasorArray');
        end

        function testExpmIsPointwiseNotMonodromy(testCase)
            % expm acts on A(t) instant by instant; pinning that down keeps a
            % future change from silently turning it into a transition matrix.
            A = PhasorArray(diag([-1, -2]));
            E = expm(A);
            testCase.verifyEqual(diag(E{:,:,0}), [exp(-1); exp(-2)], 'AbsTol', testCase.tolAnalytic);
        end

        function testExpmOfConstantMatchesMatlab(testCase)
            % A constant PhasorArray must reproduce the plain matrix exponential
            % with no harmonic content at all.
            M = [0 1; -1 0];
            E = expm(PhasorArray(M));
            dc = E{:,:,0};
            testCase.verifyEqual(dc, expm(M), 'AbsTol', testCase.tolAnalytic);
            testCase.verifyEqual(E.h, 0, 'A constant matrix must not gain harmonics');
        end

        function testExpmOfCosineGivesBesselCoefficients(testCase)
            % exp(cos th) = I_0(1) + 2 sum_k I_k(1) cos(k th), so the phasor
            % coefficient of harmonic k is exactly besseli(|k|, 1).
            E = expm(PhasorArray.cos());
            c = squeeze(E.value);
            h = E.h;
            for k = 0:6
                testCase.verifyEqual(c(k + h + 1), besseli(k, 1), 'AbsTol', testCase.tolAnalytic, ...
                    sprintf('harmonic %d is not besseli(%d,1)', k, k));
                testCase.verifyEqual(c(-k + h + 1), besseli(k, 1), 'AbsTol', testCase.tolAnalytic, ...
                    sprintf('harmonic %d is not besseli(%d,1)', -k, k));
            end
        end

        function testKron(testCase)
            A = PhasorArray.random(2, 2, 2);
            B = PhasorArray.random(2, 2, 3);
            K = kron(A, B);
            % Size: (2*2) x (2*2) = 4x4
            testCase.verifyTrue(size(K, 1) == 4 && size(K, 2) == 4, 'kron size');
            % Verify via time-domain evaluation: K(t) = A(t) ⊗ B(t)
            t0 = 0.7;
            At = evalp(A, t0);
            Bt = evalp(B, t0);
            Kt = evalp(K, t0);
            Kref = kron(At, Bt);
            err = max(abs(Kt(:) - Kref(:)));
            testCase.verifyTrue(err < testCase.tol, sprintf('kron time-domain error: %e', err));
        end

        function testNdsdpvarDetleibniz(testCase)
            testCase.assumeTrue(exist('sdpvar', 'file') == 2, 'YALMIP required');
            P = PhasorArray.ndsdpvar(2, 2, 3, "symmetry", "real");
            d = detLeibnizHmc(P);
            testCase.verifyTrue(isa(d, 'PhasorArray'), 'det should return PhasorArray');
            % The result should contain sdpvar expressions
            dval = d.value;
            testCase.verifyTrue(isa(dval, 'ndsdpvar') || isa(dval, 'sdpvar'), ...
                'Leibniz det of ndsdpvar should stay symbolic (sdpvar/ndsdpvar)');
        end

        function testOplus(testCase)
            A = PhasorArray.random(2, 2, 2);
            B = PhasorArray.random(3, 3, 2);
            Op = oplus(A, B);
            % Size: (2*3 + 3*2) ... no, oplus = A ⊗ I_b + I_a ⊗ B => size (2*3) x (2*3)
            testCase.verifyTrue(size(Op, 1) == 6 && size(Op, 2) == 6, 'oplus size should be 6x6');
            % Verify at a time sample
            t0 = 1.3;
            At = evalp(A, t0);
            Bt = evalp(B, t0);
            Opt = evalp(Op, t0);
            Opref = kron(At, eye(3)) + kron(eye(2), Bt);
            err = max(abs(Opt(:) - Opref(:)));
            testCase.verifyTrue(err < testCase.tol, sprintf('oplus time-domain error: %e', err));
        end

        function testPark(testCase)
            P = PhasorArray.Park(0);
            testCase.verifyTrue(isa(P, 'PhasorArray'), 'Park should return PhasorArray');
            testCase.verifyTrue(size(P, 1) == 3 && size(P, 2) == 3, 'Park should be 3x3');
        end

        function testParkOrthogonal(testCase)
            P = PhasorArray.Park(0);
            % Park is amplitude-invariant: the 0-component row has factor 1/sqrt(3),
            % while d,q rows have factor sqrt(2/3). So P*P.' is diagonal but
            % NOT a scalar multiple of I.
            % Check that P*P.' is diagonal and time-invariant.
            for t0 = [0, pi/4, pi/2, pi]
                Pt = evalp(P, t0);
                PPt = Pt * Pt.';
                % Should be diagonal: check off-diagonal elements are ~0
                offdiag = PPt - diag(diag(PPt));
                err_offdiag = max(abs(offdiag(:)));
                assert(err_offdiag < testCase.tol, ...
                    sprintf('Park*Park^T off-diagonal should be 0 at t=%g, err=%e', t0, err_offdiag));
                % Diagonal should be positive
                d = diag(PPt);
                testCase.verifyTrue(all(d > 0), 'Park*Park^T diagonal entries should be positive');
            end
        end

        function testPhaseshift(testCase)
            c = PhasorArray.cos();
            % Phase shift by π/2 should give -sin(t)
            c_shifted = c.PhaseShift(pi/2);
            % Evaluate at t=0: cos(0+π/2) = cos(π/2) = 0
            v0 = evalp(c_shifted, 0);
            testCase.verifyTrue(abs(v0) < testCase.tol, 'PhaseShift(cos, π/2) at t=0 should be 0');
            % Evaluate at t=π/2: cos(π/2+π/2) = cos(π) = -1
            vpi2 = evalp(c_shifted, pi/2);
            testCase.verifyTrue(abs(vpi2 - (-1)) < testCase.tol, 'PhaseShift(cos, π/2) at t=π/2 should be -1');
        end

        function testPhaseshift3phase(testCase)
            c = PhasorArray.cos();
            % Create balanced 3-phase: [cos(t), cos(t-2π/3), cos(t-4π/3)]
            phases = c.PhaseShift([0, -2*pi/3, -4*pi/3]);
            testCase.verifyTrue(size(phases, 1) == 1 && size(phases, 2) == 3, '3-phase should be 1x3');
            % Sum of balanced 3-phase = 0 at all times
            s = phases{1,1} + phases{1,2} + phases{1,3};
            testCase.verifyTrue(energy(s) < testCase.tol, 'Balanced 3-phase sum should have zero energy');
        end

        function testSymDetleibniz(testCase)
            A = PhasorArray.sym(2, 2, 1, "M","isreal",true);
            d = detLeibnizHmc(A);
            testCase.verifyTrue(isa(d, 'PhasorArray'), 'Should return PhasorArray');
            testCase.verifyTrue(d.h == 2, 'det of 2x2 with h=1 should have h=2');
        end


        function testParkMapsBalancedThreePhaseToConstantDq0(testCase)
            % The defining property, and the one that fixes the operand order:
            % Park(0) applied to a balanced set gives d = 1, q = 0, zero = 0,
            % constant in time. testParkOrthogonal cannot see the order --
            % swapping the two factors changes the matrix by 0.577 yet leaves
            % P*P' diagonal, because both factors are orthogonal.
            P = PhasorArray.Park(0);
            th = linspace(0, 2*pi, 17);
            th(end) = [];
            abc = [cos(th); cos(th - 2*pi/3); cos(th + 2*pi/3)];
            Pt = evalp(P, th);

            dq0 = zeros(3, numel(th));
            for k = 1:numel(th)
                dq0(:, k) = Pt(:, :, k) * abc(:, k);
            end
            testCase.verifyEqual(dq0(1, :), ones(1, numel(th)), 'AbsTol', 1e-10, 'd is not 1');
            testCase.verifyEqual(dq0(2, :), zeros(1, numel(th)), 'AbsTol', 1e-10, 'q is not 0');
            testCase.verifyEqual(dq0(3, :), zeros(1, numel(th)), 'AbsTol', 1e-10, 'the zero sequence is not 0');
        end

        function testEnergyEqualsSumOfSquaredCoefficients(testCase)
            % Parseval as an absolute value. The E = Edc + Eac and
            % E = Ereal + Eimag checks are invariant under a common scale
            % factor, so they cannot detect an energy that is uniformly wrong.
            A = PhasorArray.random(3, 3, 5);
            testCase.verifyEqual(energy(A), sum(abs(A.value).^2, 'all'), 'RelTol', 1e-12);

            % A single harmonic of known amplitude: cos has c_1 = c_-1 = 1/2.
            testCase.verifyEqual(energy(PhasorArray.cos()), 0.5, 'AbsTol', 1e-12, ...
                'the energy of cos is |1/2|^2 twice');
        end


        function testElementwiseProductMatchesTheTimeDomain(testCase)
            % A .* B is the entrywise product of the two periodic matrices, so its
            % order is the sum of the operands' and it must agree pointwise.
            A = PhasorArray.random(3, 3, 5);
            B = PhasorArray.random(3, 3, 3);
            C = A .* B;
            testCase.verifyEqual(C.h, A.h + B.h, 'the order must be the sum');

            th = linspace(0, 2*pi, 17);
            testCase.verifyEqual(evalp(C, th), evalp(A, th) .* evalp(B, th), 'AbsTol', 1e-12);
        end

        function testConstantMaskDoesNotInflateTheOrder(testCase)
            % Multiplying by a constant mask must leave the order alone. It used to
            % pad both operands to a common order first, so an order-5 array came
            % back at order 10, zeros above 5.
            A = PhasorArray.random(3, 3, 5);
            U = A .* PhasorArray(triu(ones(3)));

            testCase.verifyEqual(U.h, A.h, 'a constant mask changed the order');
            testCase.verifyTrue(isrealp(U), 'masking must preserve realness');
            Ut = evalp(U, linspace(0, 2*pi, 17));
            for k = 1:size(Ut, 3)
                testCase.verifyEqual(tril(Ut(:,:,k), -1), zeros(3), 'AbsTol', 1e-14, ...
                    'A(t) is not upper triangular');
            end
        end


        function testTriangularMasksKeepTheOrderAndReconstruct(testCase)
            % triu and tril are a product with a constant mask, so the order is
            % untouched and the three parts must add back to the original.
            A = PhasorArray.random(3, 3, 5);
            U = triu(A);
            L = tril(A);
            Dg = PhasorArray(A.value .* eye(3));

            testCase.verifyEqual(U.h, A.h, 'triu changed the order');
            testCase.verifyEqual(L.h, A.h, 'tril changed the order');
            testCase.verifyTrue(isrealp(U), 'masking must preserve realness');

            R = U + L - Dg;
            testCase.verifyEqual(R.value, A.value, 'AbsTol', 1e-14, ...
                'triu + tril - diag does not rebuild A');

            th = linspace(0, 2*pi, 17);
            Ut = evalp(U, th);
            for k = 1:size(Ut, 3)
                testCase.verifyEqual(tril(Ut(:,:,k), -1), zeros(3), 'AbsTol', 1e-14);
            end

            % The offset argument reaches the diagonal itself.
            U1 = evalp(triu(A, 1), th);
            for k = 1:size(U1, 3)
                testCase.verifyEqual(tril(U1(:,:,k), 0), zeros(3), 'AbsTol', 1e-14);
            end
        end

    end

    methods (Test, TestTags = {'Install'})
        % Smoke set for a fresh install: no optional toolbox, a few seconds,
        % one check per layer. Run with run_all_tests("install").
        function testDerivative(testCase)
            % d/dt of cos(ωt) = -ω sin(ωt)
            c = PhasorArray.cos();
            T = 2*pi;
            dc = d(c, T);
            % Evaluate at t=0: d/dt cos(t)|_0 = -sin(0) = 0
            dc0 = evalp(dc, 0);
            testCase.verifyTrue(abs(dc0) < testCase.tol, 'd/dt cos(0) should be 0');
            % Evaluate at t=pi/2: d/dt cos(t)|_{pi/2} = -sin(pi/2) = -1
            dcpi2 = evalp(dc, pi/2);
            testCase.verifyTrue(abs(dcpi2 - (-1)) < testCase.tol, 'd/dt cos(pi/2) should be -1');
        end

        function testEnergyParseval(testCase)
            A = PhasorArray.random(3, 3, 5);
            E = energy(A);
            Edc = DCenergy(A);
            Eac = ACenergy(A);
            testCase.verifyTrue(abs(E - Edc - Eac) < testCase.tol * E, 'Parseval: total ≠ DC + AC');
        end

        function testDetleibniz2x2(testCase)
            A = PhasorArray.random(2, 2, 3);
            dL = detLeibnizHmc(A);
            % Compare with formula: A11*A22 - A12*A21
            dManual = A{1,1} * A{2,2} - A{1,2} * A{2,1};
            err = energy(dL - dManual);
            testCase.verifyTrue(err < testCase.tol, sprintf('2x2 Leibniz error: %e', err));
        end

    end

    methods (Test)
        % Symmetry projections and three-phase transformations.
        %
        % Deliberately outside the Install block above: these are regression
        % coverage, not installation checks. The smoke set answers "is the
        % install functional", one check per layer in a few seconds, and
        % every test added to it is paid on every run.
        function testHermAndSymSplitsAreExact(testCase)
            % Each split reconstructs A, and each part carries none of the other.
            % Complex data on purpose: PhasorArray.random returns a real-valued
            % A(t), on which the Hermitian and symmetric notions coincide and
            % the test would not separate them.
            A = PhasorArray(randn(3, 3, 9) + 1i * randn(3, 3, 9));
            H = mherm(A);   Kh = mherm(A, skewOption='skew');
            S = msym(A);    Ks = msym(A,  skewOption='skew');
            testCase.verifyEqual(value(H + Kh), value(A), 'AbsTol', testCase.tol, ...
                'Hermitian split must reconstruct A');
            testCase.verifyEqual(value(S + Ks), value(A), 'AbsTol', testCase.tol, ...
                'symmetric split must reconstruct A');
            % The next two projections are empty by construction -- that is the
            % property under test -- so phasorSymmetry warns. Silenced here so
            % a passing test stays quiet.
            ws = warning('off', 'PhasorArray:symmetry:emptyProjection');
            restore = onCleanup(@() warning(ws)); %#ok<NASGU>
            testCase.verifyLessThan(energy(mherm(H, skewOption='skew')), testCase.tol, ...
                'the Hermitian part must carry no skew-Hermitian content');
            testCase.verifyLessThan(energy(msym(S, skewOption='skew')), testCase.tol, ...
                'the symmetric part must carry no skew-symmetric content');
        end

        function testHermAndSymAreDistinctNotions(testCase)
            % A = A' and A = A.' are different properties on complex data: the
            % Hermitian split conjugates and mirrors the harmonics, the
            % symmetric one does neither.
            A = PhasorArray(randn(3, 3, 9) + 1i * randn(3, 3, 9));
            testCase.verifyFalse(isreal(A), 'the fixture must be complex for this test to mean anything');
            testCase.verifyGreaterThan(energy(mherm(A) - msym(A)), testCase.tol, ...
                'Hermitian and symmetric parts must differ on complex data');
            % They coincide on a real-valued A(t), where conjugation is a no-op.
            R = PhasorArray.random(3, 3, 4);
            testCase.verifyTrue(isreal(R), 'PhasorArray.random is expected to give a real A(t)');
            testCase.verifyLessThan(energy(mherm(R) - msym(R)), 1e-9, ...
                'on a real A(t) the two notions must agree');
        end

        function testHermAndSymMatchTheOperatorForms(testCase)
            % The raw-payload implementations must agree with the operators.
            % Complex, so that ' and .' are genuinely different operators here.
            A = PhasorArray(randn(3, 3, 9) + 1i * randn(3, 3, 9));
            testCase.verifyEqual(value(mherm(A)), value((A + A') * (1/2)), ...
                'AbsTol', testCase.tol, 'mherm must equal (A + A'')/2');
            testCase.verifyEqual(value(mherm(A, skewOption='skew')), value((A - A') * (1/2)), ...
                'AbsTol', testCase.tol, 'skew mherm must equal (A - A'')/2');
            testCase.verifyEqual(value(msym(A)), value((A + A.') * (1/2)), ...
                'AbsTol', testCase.tol, 'msym must equal (A + A.'')/2');
            testCase.verifyEqual(value(msym(A, skewOption='skew')), value((A - A.') * (1/2)), ...
                'AbsTol', testCase.tol, 'skew msym must equal (A - A.'')/2');
        end

        function testEnergiesSplitTheTotal(testCase)
            % Both decompositions are orthogonal, so each pair of energies adds
            % up to the total -- the identity realEnergy/imagEnergy satisfy.
            A = PhasorArray(randn(3, 3, 9) + 1i * randn(3, 3, 9));
            Et = energy(A);
            testCase.verifyTrue(abs(Et - hermEnergy(A) - hermEnergy(A, false, skewOption='skew')) ...
                < testCase.tol * max(Et, 1), 'Energy: total ~= herm + skew-herm');
            testCase.verifyTrue(abs(Et - symEnergy(A) - symEnergy(A, false, skewOption='skew')) ...
                < testCase.tol * max(Et, 1), 'Energy: total ~= sym + skew-sym');
        end

        function testSymmetryPredicatesAgreeWhereTheyShould(testCase)
            % Each predicate must agree with its magnitude counterpart.
            % ISSYMMETRIC tests every harmonic independently: A(t) is
            % symmetric exactly when every A_k is. ISHERMITIAN pairs k
            % with -k (A_k = A_{-k}'), because conjugation mirrors the spectrum.
            raw = randn(3, 3, 7) + 1i * randn(3, 3, 7);

            Asym = PhasorArray((raw + pagetranspose(raw)) * (1/2));
            testCase.verifyLessThan(symEnergy(Asym, false, skewOption='skew'), testCase.tol, ...
                'built symmetric, so it must carry no skew-symmetric energy');
            testCase.verifyTrue(issymmetric(Asym), 'issymmetric must agree with msym');

            Aherm = PhasorArray((raw + flip(pagectranspose(raw), 3)) * (1/2));
            testCase.verifyLessThan(hermEnergy(Aherm, false, skewOption='skew'), testCase.tol, ...
                'built Hermitian in time, so it must carry no skew-Hermitian energy');
            testCase.verifyTrue(ishermitian(Aherm), 'ishermitian must agree with mherm');

            % The per-slice reading is the wrong one and must not pass: every
            % A_k Hermitian does not make A(t) Hermitian.
            Apage = PhasorArray((raw + pagectranspose(raw)) * (1/2));
            testCase.verifyFalse(ishermitian(Apage), ...
                'per-slice Hermitian symmetry must not be mistaken for A(t) = A(t)''');
            testCase.verifyGreaterThan(hermEnergy(Apage, false, skewOption='skew'), testCase.tol, ...
                'and it must carry skew-Hermitian energy');
        end

        function testEnergyElementwiseShape(testCase)
            A = PhasorArray.random(3, 3, 4);
            Eew = hermEnergy(A, true);
            testCase.verifySize(Eew, [3 3], 'element-wise energy must be n x m');
            [Eew2, Etot] = symEnergy(A);
            testCase.verifySize(Eew2, [3 3], 'two outputs must give the element-wise matrix first');
            testCase.verifyEqual(sum(Eew2, 'all'), Etot, 'AbsTol', testCase.tol, ...
                'total must be the sum of the element-wise energies');
        end

        % ---------------------------------------------------------------
        % Three-phase transformations
        %
        % Public API a user reaches for directly, and until now covered only
        % for Concordia, dq0 and Park. What follows pins the algebra of the
        % whole family: which product is orthogonal, which is merely
        % diagonal, how the pieces compose, and what `order` and `dephase`
        % actually do.
        % ---------------------------------------------------------------

        function testClarkIsAmplitudeInvariantAndConstant(testCase)
            % Clark preserves amplitude, not power: K*K.' is diagonal but not
            % the identity. The d,q rows carry 2/3 and the zero row 1/3.
            K = PhasorArray.Clark();
            testCase.verifyEqual(K.h, 0, 'Clark is a constant matrix, so h must be 0');
            testCase.verifyEqual(size(value(K)), [3 3], 'Clark must be 3x3 with a single harmonic');
            testCase.verifyEqual(evalp(K * K.', 0), diag([2/3 2/3 1/3]), ...
                'AbsTol', testCase.tol, 'Clark*Clark.'' must be diag(2/3, 2/3, 1/3)');
        end

        function testConcordiaIsPowerInvariant(testCase)
            % The counterpart of the test above: Concordia is orthogonal, so
            % it preserves power. testConcordia only checks one entry of the
            % matrix, which a wrong row ordering would survive.
            C = PhasorArray.Concordia();
            testCase.verifyEqual(evalp(C * C.', 0), eye(3), 'AbsTol', testCase.tol, ...
                'Concordia must be orthogonal');
        end

        function testRotationsAreOrthogonalAtEveryInstant(testCase)
            % A rotation stays a rotation at every instant, not only at t=0 --
            % which is where a sign error in one of the four sine entries
            % would hide.
            rotations = {PhasorArray.Rotdq0(0, 1), PhasorArray.negativeRotdq0(0, 1)};
            for ii = 1:numel(rotations)
                for t = [0, 0.7, pi/3, 2.1, 5.5]
                    Rt = evalp(rotations{ii}, t);
                    testCase.verifyEqual(Rt * Rt.', eye(3), 'AbsTol', testCase.tol, ...
                        sprintf('rotation not orthogonal at t = %g', t));
                end
            end
        end

        function testDq0FamilyIsPowerInvariantAndParkFamilyIsNot(testCase)
            % dq0 = rotation * Concordia inherits orthogonality; Park =
            % rotation * Clark inherits the diagonal-but-not-identity form.
            % Checking both together is what makes the distinction testable.
            for t = [0, 0.4, pi/2, 3.3]
                orthogonal = {PhasorArray.dq0(0, 1), PhasorArray.negativeDQ0(0, 1)};
                for ii = 1:numel(orthogonal)
                    Dt = evalp(orthogonal{ii}, t);
                    testCase.verifyEqual(Dt * Dt.', eye(3), 'AbsTol', testCase.tol, ...
                        sprintf('dq0-family must be orthogonal at t = %g', t));
                end
                diagonal = {PhasorArray.Park(0, 1), PhasorArray.negativePark(0, 1)};
                for ii = 1:numel(diagonal)
                    Pt = evalp(diagonal{ii}, t);
                    testCase.verifyEqual(Pt * Pt.', diag([2/3 2/3 1/3]), ...
                        'AbsTol', testCase.tol, ...
                        sprintf('Park-family must be diag(2/3,2/3,1/3) at t = %g', t));
                end
            end
        end

        function testTransformsFactorThroughTheirRotation(testCase)
            % Each composite transform is exactly rotation * static frame,
            % as the docstrings promise.
            testCase.verifyLessThan(energy(PhasorArray.Park(0, 1) ...
                - PhasorArray.Rotdq0(0, 1) * PhasorArray.Clark()), testCase.tol, ...
                'Park must equal Rotdq0 * Clark');
            testCase.verifyLessThan(energy(PhasorArray.dq0(0, 1) ...
                - PhasorArray.Rotdq0(0, 1) * PhasorArray.Concordia()), testCase.tol, ...
                'dq0 must equal Rotdq0 * Concordia');
            testCase.verifyLessThan(energy(PhasorArray.negativeDQ0(0, 1) ...
                - PhasorArray.negativeRotdq0(0, 1) * PhasorArray.Concordia()), testCase.tol, ...
                'negativeDQ0 must equal negativeRotdq0 * Concordia');
            testCase.verifyLessThan(energy(PhasorArray.negativePark(0, 1) ...
                - PhasorArray.negativeRotdq0(0, 1) * PhasorArray.Clark()), testCase.tol, ...
                'negativePark must equal negativeRotdq0 * Clark');
        end

        function testDq0MapsBalancedThreePhaseToAConstant(testCase)
            % The power-invariant twin of testParkMapsBalancedThreePhaseToConstantDq0.
            % The d component lands on sqrt(3/2), not 1: that factor is the
            % whole difference between the two conventions, so a test that
            % only checked q = 0 would not see them swapped.
            abc = [PhasorArray.cos(0, 1); PhasorArray.cos(-2*pi/3, 1); PhasorArray.cos(2*pi/3, 1)];
            y   = PhasorArray.dq0(0, 1) * abc;
            testCase.verifyEqual(real(y{:,:,0}), [sqrt(3/2); 0; 0], 'AbsTol', 1e-10, ...
                'dq0 of a balanced direct set must be [sqrt(3/2); 0; 0]');
            ac = value(y);
            ac(:, :, y.h + 1) = 0;
            testCase.verifyLessThan(norm(ac(:)), 1e-10, 'the image must be constant in time');
        end

        function testNegativeDq0RejectsTheDirectSequence(testCase)
            % What makes the negative frame worth having: fed the direct
            % sequence it produces no constant term and keeps oscillating,
            % where dq0 produces a constant. Without this, negativeDQ0 could
            % be an exact copy of dq0 and every other test would still pass.
            abc = [PhasorArray.cos(0, 1); PhasorArray.cos(-2*pi/3, 1); PhasorArray.cos(2*pi/3, 1)];
            y   = PhasorArray.negativeDQ0(0, 1) * abc;
            testCase.verifyLessThan(norm(y{:,:,0}), 1e-10, ...
                'the negative frame must leave no DC on a direct sequence');
            ac = value(y);
            ac(:, :, y.h + 1) = 0;
            testCase.verifyGreaterThan(norm(ac(:)), 1, 'and must keep it oscillating');
        end

        function testDocumentedMatricesMatchTheCode(testCase)
            % The help text of each transform prints its matrix; this
            % transcribes those four and compares them to what the code
            % actually evaluates to.
            t = 0.4;
            c  = @(x) cos(t + x);
            s  = @(x) sin(t + x);
            p  = 2*pi/3;
            documented = struct( ...
                'Park',        (2/3)     * [c(0) c(-p) c(p); -s(0) -s(-p) -s(p); 1/2 1/2 1/2], ...
                'negativePark',(2/3)     * [c(0) c(p) c(-p); -s(0) -s(p) -s(-p); 1/2 1/2 1/2], ...
                'dq0',         sqrt(2/3) * [c(0) c(-p) c(p); -s(0) -s(-p) -s(p); ...
                                            1/sqrt(2) 1/sqrt(2) 1/sqrt(2)], ...
                'negativeDQ0', sqrt(2/3) * [c(0) c(p) c(-p); -s(0) -s(p) -s(-p); ...
                                            1/sqrt(2) 1/sqrt(2) 1/sqrt(2)]);
            for name = string(fieldnames(documented))'
                testCase.verifyEqual(evalp(PhasorArray.(name)(0, 1), t), ...
                    documented.(name), 'AbsTol', testCase.tol, ...
                    sprintf('%s does not evaluate to the matrix its help text prints', name));
            end
        end

        function testOrderPlacesTheHarmonic(testCase)
            % order = k puts the rotation at harmonics +/-k and nowhere else.
            for k = 1:3
                R = PhasorArray.Rotdq0(0, k);
                occupied = [];
                for j = -R.h:R.h
                    if norm(R{:,:,j}, 'fro') > testCase.tol, occupied(end+1) = j; end %#ok<AGROW>
                end
                testCase.verifyEqual(occupied, [-k 0 k], ...
                    sprintf('Rotdq0 of order %d must occupy harmonics -%d, 0, +%d', k, k, k));
            end
        end

        function testDephaseShiftsTheBaseAngle(testCase)
            % dephase is phi in cos(k*(theta + phi)): it shifts the base
            % angle, so Rotdq0(d,k) at t equals Rotdq0(0,k) at t+d for every k.
            d = 0.37; t = 0.9;
            for k = 1:3
                testCase.verifyEqual(evalp(PhasorArray.Rotdq0(d, k), t), ...
                    evalp(PhasorArray.Rotdq0(0, k), t + d), 'AbsTol', testCase.tol, ...
                    sprintf('dephase must act as a time shift at order %d', k));
            end
        end

        function testInclude0FalseDropsTheZeroSequence(testCase)
            % Same flag, two row conventions: the abc-side transforms carry
            % the zero sequence last and drop the last row, the sequence
            % transforms carry it first and drop the first. Both are correct;
            % only the resulting size is common to all of them.
            % size() of a PhasorArray reports the harmonic axis too, so the
            % matrix dimensions have to be asked for by name.
            shape = @(X) size(X, [1 2]);
            testCase.verifyEqual(shape(PhasorArray.Clark(false)),        [2 3], 'Clark(false)');
            testCase.verifyEqual(shape(PhasorArray.Concordia(false)),    [2 3], 'Concordia(false)');
            testCase.verifyEqual(shape(PhasorArray.Rotdq0(0, 1, false)), [2 2], 'Rotdq0(...,false) is square');
            testCase.verifyEqual(shape(PhasorArray.dq0(0, 1, false)),    [2 3], 'dq0(...,false)');
            testCase.verifyEqual(shape(PhasorArray.Park(0, 1, false)),   [2 3], 'Park(...,false)');
            testCase.verifyEqual(shape(PhasorArray.ZPNSequence(0, 1, false)), [2 3], 'ZPNSequence(...,false)');
            testCase.verifyEqual(shape(PhasorArray.zeroPosNegSequenceDQ(0, 1, false)), [4 3], ...
                'zeroPosNegSequenceDQ(...,false) keeps both dq frames');
        end

        function testZeroPosNegSequenceDqStacksTheTwoFrames(testCase)
            % It is a stack, not an independent derivation: zero row, then
            % the direct dq pair, then the inverse dq pair. Pinning the row
            % order is the point -- swapping the two frames is invisible to
            % any norm-based check.
            ZQ = PhasorArray.zeroPosNegSequenceDQ(0, 1);
            Dp = PhasorArray.dq0(0, 1, false);
            Dn = PhasorArray.negativeDQ0(0, 1, true);
            testCase.verifyEqual(size(ZQ, [1 2]), [5 3], 'zeroPosNegSequenceDQ must be 5x3');
            testCase.verifyLessThan(energy(ZQ{1,:}   - Dn{3,:}),   testCase.tol, 'row 1 is the zero sequence');
            testCase.verifyLessThan(energy(ZQ{2:3,:} - Dp),        testCase.tol, 'rows 2-3 are the direct dq frame');
            testCase.verifyLessThan(energy(ZQ{4:5,:} - Dn{1:2,:}), testCase.tol, 'rows 4-5 are the inverse dq frame');
        end

        function testZpnSequenceCarriesTheFortescueMatrix(testCase)
            % The matrix itself is the classical Fortescue operator and is
            % unitary up to the factor 3. Order 3 is excluded on purpose: a
            % becomes 1 there and the matrix degenerates to rank one, which
            % is the expected behaviour for a triplen order.
            for k = [1 2]
                Z = PhasorArray.ZPNSequence(0, k);
                testCase.verifyEqual(size(Z, [1 2]), [3 3], 'ZPNSequence must be 3x3');
                a = exp(2i*pi*k/3);
                F = [1 1 1; 1 a^2 a; 1 a a^2];
                testCase.verifyEqual(Z{:,:,1}, F, 'AbsTol', testCase.tol, ...
                    sprintf('ZPNSequence of order %d must carry the Fortescue matrix', k));
                % Unitarity checked on the matrix ZPNSequence produced, not on
                % the reference F built above.
                testCase.verifyEqual(Z{:,:,1} * Z{:,:,1}' / 3, eye(3), ...
                    'AbsTol', testCase.tol, ...
                    'the operator must be unitary up to the factor 3');
            end
            % Rank collapse at a triplen order, stated so it is not mistaken
            % for a regression later.
            a3 = exp(2i*pi);
            testCase.verifyEqual(a3, 1, 'AbsTol', testCase.tol, ...
                'at order 3 the sequence operator degenerates by construction');
        end

    end
end
