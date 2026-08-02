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
end
