classdef PhasorArraySolversTest < matlab.unittest.TestCase
    %PHASORARRAYSOLVERSTEST  Lyapunov, Sylvester, Riccati, pole placement and their LMI counterparts.
    %
    %   Each solver is held to the residual of the equation it claims to solve,
%   recomputed from the returned solution rather than read back from the solver.
%   Positive definiteness alone hides a wrong equation, which is how the ACDC
%   template ran for months with a period where a pulsation was due.

    properties
        tol = 1e-10;
        T   = 2*pi;
        h   = 10;
        A
        Q
    end

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            srcFolder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcFolder, 'IncludingSubfolders', true));
        end
    end


    methods (TestMethodSetup)
        function buildSystem(testCase)
            A0 = [-3 0.5; 0 -2];
            A1 = [0.1 0; 0 0.05];
            A2 = [0 0.3; 0.1 0];
            testCase.A = PhasorArray(cat(3, conj(A2), conj(A1), A0, A1, A2));
            testCase.Q = PhasorArray(eye(2));
        end
    end

    methods (Access = private)
        function r = residual(testCase, X)
            % Toeplitz residual of dP + A'P + PA + Q = 0, dP/dt = [N, P].
            hh = testCase.h;
            N  = N_tb(size(testCase.A, 1), hh, testCase.T);
            XT = X.T_tb(hh);
            r  = norm(T_tb(testCase.A' * X + X * testCase.A, hh) ...
                      + testCase.Q.T_tb(hh) + (N*XT - XT*N), 'fro');
        end

    end

    methods (Test)
        function testAdaptivehInfofields(testCase)
            % All three solvers now share the driver, so all three must publish the
            % same diagnostics contract (mlHmcDivide gained time_history in the process).
            need = {'status','statusMsg','resrelnorm','resnorm','h', ...
                    'h_history','resrel_history','res_history','time_history', ...
                    'regime_history','s_alg_history','s_exp_history'};

            A = PhasorArray.random(2,2,2) - 2*PhasorArray.eye(2);
            [~, iL] = lyap(A, PhasorArray.eye(2), autoUpdateh=true, maxh=8);

            E = 0.1*PhasorArray.random(2,2,1) + PhasorArray.eye(2);
            [~, iG] = lyapG(A, PhasorArray.eye(2), E, autoUpdateh=true, maxh=8);

            Ad = 0.2*PhasorArray.random(2,2,1) + PhasorArray.eye(2);
            [~, iD] = mlHmcDivide(Ad, PhasorArray.random(2,1,1), autoUpdateh=true, maxh=8);

            for s = {iL, iG, iD}
                missing = need(~isfield(s{1}, need));
                testCase.verifyTrue(isempty(missing), sprintf('info is missing: %s', strjoin(missing, ', ')));
                assert(numel(s{1}.h_history) == numel(s{1}.time_history), ...
                    'h_history and time_history lengths disagree.');
                testCase.verifyTrue(ismember(s{1}.status, [0 1 2 4]), ...
                    sprintf('Unexpected status %d from an adaptive run.', s{1}.status));
            end
        end

        function testAdaptivehProgress(testCase)
            % The refinement order must strictly increase at every step, whatever the
            % step rule decides. This is the guard that makes a non-terminating loop
            % structurally impossible; it is exercised here through a deliberately
            % pathological callback whose residual never improves, so the step rule
            % lands in the 'stagnated' regime and proposes deltah = 1 forever.
            stuck = @(h) deal(PhasorArray.eye(2), 1, 1, PhasorArray.eye(2));  % residual never moves
            cfg = struct('thresholdResidual', 1e-12, 'maxh', 8, ...
                'stagnationWindow', 100, 'stagnationRatio', 0, ...
                'updateMethod', 'adaptive', 'verbose', false, 'hOp', 1, ...
                'hOutFcn', @(h) h, 'preamble', '', 'label', 'h');

            [best, trace] = adaptiveHSolve(stuck, 0, cfg);

            assert(all(diff(trace.h_history) >= 1), ...
                sprintf('h_history is not strictly increasing: %s', mat2str(trace.h_history)));
            testCase.verifyTrue(trace.h_history(end) <= 8, 'maxh was exceeded.');
            testCase.verifyTrue(trace.status == 2, sprintf('Expected status 2 (maxh), got %d.', trace.status));
            testCase.verifyTrue(best.h == trace.h_history(1), 'Best iterate should be the first (no improvement occurred).');
        end

        function testAdaptivehStatus(testCase)
            % Convergence must be reported as status 0 and must publish the CONVERGED
            % iterate, not merely the smallest-residual one seen so far.
            % Residual model: e(h) = 10^-h, so the threshold is crossed at a known h.
            solve = @(h) deal(PhasorArray.eye(2)*(h+1), 10^-h, 10^-h, PhasorArray.eye(2));
            cfg = struct('thresholdResidual', 1e-3, 'maxh', 20, ...
                'stagnationWindow', 100, 'stagnationRatio', 0, ...
                'updateMethod', 'incremental', 'verbose', false, 'hOp', 0, ...
                'hOutFcn', @(h) h, 'preamble', '', 'label', 'h');

            [best, trace] = adaptiveHSolve(solve, 0, cfg);

            testCase.verifyTrue(trace.status == 0, sprintf('Expected status 0 (converged), got %d.', trace.status));
            testCase.verifyTrue(best.resrelnorm <= 1e-3, 'Returned iterate does not meet the threshold.');
            testCase.verifyTrue(best.h == 3, sprintf('Expected convergence at h=3, got h=%d.', best.h));
            % best.sol must be the solution computed AT best.h
            testCase.verifyTrue(isequal(value(best.sol), value(PhasorArray.eye(2)*(best.h+1))), ...
                'best.sol does not correspond to best.h.');
            testCase.verifyTrue(numel(trace.h_history) == numel(trace.resrel_history), 'History lengths disagree.');
            testCase.verifyTrue(numel(trace.h_history) == numel(trace.time_history), 'time_history length disagrees.');
        end

        function testCheckSolvers(testCase)
            % Checks if YALMIP and MOSEK are detected
            has_sdpvar = exist('sdpvar', 'file') == 2;
            if has_sdpvar
                fprintf('    [INFO] YALMIP detected.\n');
            else
                fprintf('    [INFO] YALMIP NOT detected.\n');
            end
    
            has_mosek = exist('mosekopt', 'file') == 3;
            if has_mosek
                fprintf('    [INFO] MOSEK detected.\n');
            else
                fprintf('    [INFO] MOSEK NOT detected.\n');
            end
        end

        function testHareStabilisesTheClosedLoop(testCase)
            A0 = [-3 0.5; 0 -2];
            A1 = [0.1 0; 0 0.05];
            A  = PhasorArray(cat(3, conj(A1), A0, A1));
            B  = PhasorArray([1; 0.5]);
            Q  = PhasorArray(diag([10, 1]));
            R  = PhasorArray(1);

            [K, X] = hare(A, B, Q, R);
            testCase.verifyClass(K, 'PhasorArray');

            h = 12;
            Acl = A - B * K;
            e = eig(Acl.T_tb(h) - N_tb(2, h, 2*pi));
            testCase.verifyLessThan(max(real(e)), 0, 'the closed loop is not stable');

            % K must be the gain the returned X induces: K = R^-1 B' X.
            Kfromx = R \ (B.' * X);
            relErr = norm(K.T_tb(h) - Kfromx.T_tb(h), 'fro') / (norm(K.T_tb(h), 'fro') + eps);
            % Kleinman stops on its own residual threshold, so the gain matches
            % the returned X to the solver tolerance, not to eps.
            testCase.verifyLessThan(relErr, 1e-4, 'K does not match R^-1 B'' X');
        end

        function testLmiRouteReachesTheSameResidual(testCase)
            testCase.assumeTrue(exist('sdpvar', 'file') == 2, 'YALMIP required');
            hh = testCase.h;
            nx = size(testCase.A, 1);

            P  = PhasorArray.ndsdpvar(nx, nx, hh);
            PT = P.T_tb(hh);
            N  = N_tb(nx, hh, testCase.T);
            Op = T_tb(testCase.A' * P + P * testCase.A, hh) + testCase.Q.T_tb(hh) + (N*PT - PT*N);

            sol = optimize([Op <= 0, PT >= 0], trace(PT), sdpsettings('verbose', 0));
            testCase.assumeTrue(sol.problem == 0, sprintf('solver failed: %s', sol.info));

            X_lmi  = sdpval(P);
            X_lyap = lyap(testCase.A, testCase.Q, 'h', hh, 'autoUpdateh', true, 'T', testCase.T);

            testCase.verifyLessThan(testCase.residual(X_lmi), 1e-4, ...
                'the LMI optimum does not satisfy the Lyapunov equation');
            testCase.verifyLessThan( ...
                norm(X_lmi.T_tb(hh) - X_lyap.T_tb(hh), 'fro') / norm(X_lyap.T_tb(hh), 'fro'), 1e-5, ...
                'the two routes disagree on the solution');
        end


        function testLyapLtp(testCase)
            % LTP Lyapunov: dP + A'(t)P + PA(t) + Q = 0 with periodic A(t)
            A0 = [-3 0.5; 0 -2];
            A1 = [0.1 0; 0 0.05];
            A2 = [0 0.3; 0.1 0];
            A = PhasorArray(cat(3, conj(A2), conj(A1), A0, A1, A2));  % h=1
            Q = PhasorArray(eye(2));
    
            [X, residual] = lyap(A, Q, 'h', 2,'autoUpdateh',true);
    
            % Residual should be moderate (truncation at h=10)
            assert(residual.resnorm < 1, ...
                sprintf('LTP Lyapunov residual: %e', residual.resnorm));
    
            E = energy(X);
            testCase.verifyTrue(isfinite(E) && E > 0, 'LTP Lyapunov solution should have finite positive energy');
    
            % Symmetry residual should be small (P should be symmetric)
            assert(residual.resPsym < 1e-6, ...
                sprintf('LTP Lyapunov symmetry residual: %e', residual.resPsym));
        end

        function testLyapSolutionIsSymmetricAndPositive(testCase)
            X = lyap(testCase.A, testCase.Q, 'h', testCase.h, 'autoUpdateh', true, 'T', testCase.T);
            XT = X.T_tb(testCase.h);
            testCase.verifyLessThan(norm(XT - XT', 'fro') / norm(XT, 'fro'), 1e-8, ...
                'P is not Hermitian');
            testCase.verifyGreaterThan(min(real(eig(XT))), 0, ...
                'P is not positive definite for a stable A and Q = I');
        end


        function testPlaceAssignsTheRequestedFloquetExponents(testCase)
            A0 = [0 1; -2 -3];
            A1 = [0.05 0; 0 0.05];
            A  = PhasorArray(cat(3, conj(A1), A0, A1));
            B  = PhasorArray([0; 1]);
            mus = [-4; -5];   % away from the open-loop poles {-1,-2}: sharing them makes the Sylvester equation singular

            % Default settings. hG used to default to A.h, which scattered the
            % placed exponents between 1e-2 and 7e-1 depending on the random draw
            % of G, silently. A constant G lands within 1e-9 on every draw.
            [K, ~, res] = place(A, B, mus);
            testCase.verifyLessThan(sqrt(energy(res)), 1e-6, ...
                'the Sylvester residual place reports is large');
            Acl = A - B * K;
            h = 20;
            e = eig(Acl.T_tb(h) - N_tb(2, h, 2*pi));
            err = max(arrayfun(@(m) min(abs(e - m)), mus));
            testCase.verifyLessThan(err, 1e-6, ...
                sprintf('placed exponents are off by %.2e', err));

            % Stable across draws: the old default was not.
            for seed = 1:4
                rng(seed);
                Kk = place(A, B, mus);
                Ak = A - B * Kk;
                ek = eig(Ak.T_tb(h) - N_tb(2, h, 2*pi));
                testCase.verifyLessThan(max(arrayfun(@(m) min(abs(ek - m)), mus)), 1e-6, ...
                    sprintf('placement depends on the random G (seed %d)', seed));
            end
        end

        function testRandomSpdRespectsTheLowerBound(testCase)
            alpha = 0.5;
            P = PhasorArray.randomSPD(3, 4, 'alphaMin', alpha, 'seed', 1);
            th = linspace(0, 2*pi, 201);
            Pt = evalp(P, th);
            lam = arrayfun(@(k) min(real(eig(Pt(:,:,k)))), 1:numel(th));
            testCase.verifyGreaterThanOrEqual(min(lam), alpha - 1e-8, ...
                sprintf('P(t) - alpha I is not positive semidefinite, min eig %.4f', min(lam)));
            testCase.verifyTrue(isrealp(P), 'randomSPD must return a real signal');
        end

        function testRandomwithnpoleScalar(testCase)
            % Scalar case: exponent = DC coefficient, output must be a genuinely
            % periodic real signal for every Method.
            mu = -0.7;
            h  = 8;
            tol = 1e-8;
            for method = ["structured", "generic", "truncated"]
                a = PhasorArray.randomWithNPole(mu, h, Method=method, seed=1);
                v = squeeze(a.value);
                % The qualification diagnostic is a single string and comes before any
                % name-value pair; it is not a printf-style format plus arguments.
                testCase.verifyTrue(~any(isnan(v(:))), sprintf('NaN output (%s)', method));
                testCase.verifyEqual(abs(v(h+1) - mu), 0, sprintf('DC ~= mu (%s)', method), 'AbsTol', tol);
                testCase.verifyTrue(a.h == h, sprintf('wrong bandwidth (%s)', method));
                testCase.verifyTrue(norm(v(h+2:end)) > 0, sprintf('constant output, no AC part (%s)', method));
                testCase.verifyTrue(isreal(a), sprintf('a(t) not real (%s)', method));
            end
        end

        function testResidualRejectsAWrongPeriod(testCase)
            % The check that the ACDC template needed: solving with the period
            % where the pulsation belongs still returns a positive definite P,
            % so only the residual tells the two apart.
            X_ok    = lyap(testCase.A, testCase.Q, 'h', testCase.h, 'T', 2*pi);
            X_wrong = lyap(testCase.A, testCase.Q, 'h', testCase.h, 'T', 1);

            testCase.verifyLessThan(testCase.residual(X_ok), 1e-5);
            testCase.verifyGreaterThan(testCase.residual(X_wrong), 1e-3, ...
                'a wrong period must show up in the residual');
        end

        function testRiccatiLmiVsKleinman(testCase)
        % Verifies that RicHarmonicKlein converges to the same LQR gain as the
        % direct Schur-complement Riccati LMI (YALMIP).
        %
        % HARE:  (A-N)^* P + P(A-N) - P B R^{-1} B^* P + Q = 0
        % Schur LMI:
        %   [ He((A-N)^*P) + Q,  PB ] >= 0,   max tr(PT)
        %   [ B^*P,               R  ]
        %
        % System: stable open-loop (A0 has eigs -3, -2) with periodic perturbation.
        % K0=0 is valid; Lyapunov solver converges cleanly.

            T  = 2*pi;
            h  = 8;
            A0 = [-3 0.5; 0 -2];
            A1 = [0.1 0; 0 0.05];
            A  = PhasorArray(cat(3, conj(A1), A0, A1));
            B  = PhasorArray([1; 0.5]);
            Q  = PhasorArray(diag([10, 1]));
            R  = PhasorArray(1);

            % --- 1. RicHarmonicKlein (K0=0 valid since A is open-loop stable) ---
            [K_klein] = RicHarmonicKlein(A, B, Q, R, ...
                PhasorArray(zeros(1, 2)), T, ...
                'maxIter', 200, 'thresholdResidual', 1e-8, 'autoUpdateh', true);

            % --- 2. Riccati Schur LMI (wiki template §4) ---
            P  = PhasorArray.ndsdpvar(2, 2, h);
            PT = P.T_tb(h);
            N  = N_tb(2, h, T);

            % He((A-N)^* P) + Q  in Toeplitz block
            % (A-N)^* = A^* + N  since N^* = -N (skew-Hermitian)
            % => T_tb(A'P + PA) + N*PT - PT*N + Q_tb
            F11 = T_tb(A'*P + P*A, h) + (N*PT - PT*N) + Q.T_tb(h);
            F12 = T_tb(P*B, h);
            F22 = R.T_tb(h);

            F = [PT >= 0, [F11, F12; F12', F22] >= 0];
            sol = optimize(F, -trace(PT), sdpsettings('solver', 'mosek', 'verbose', 0));
            testCase.verifyTrue(sol.problem == 0, sprintf('YALMIP Riccati LMI failed: %s', sol.info));

            P_opt = sdpval(P);
            K_lmi = R \ (B.' * P_opt);

            % --- 3. Compare at truncation h ---
            rel_err = norm(K_lmi.T_tb(h) - K_klein.T_tb(h), 'fro') / ...
                      (norm(K_lmi.T_tb(h), 'fro') + eps);
            fprintf('    [INFO] Riccati LMI vs Kleinman: rel_err = %.2e\n', rel_err);
            assert(rel_err < 5e-2, ...
                sprintf('RicHarmonicKlein vs Schur LMI mismatch: %.2e', rel_err));
        end

        function testSylvester(testCase)
            % Sylvester: dM + A(t)M + MB(t) + C = 0
            A = PhasorArray([-2 0; 0 -3]);
            B = PhasorArray([-1 0; 0 -4]);
            C = PhasorArray(randn(2, 2));
    
            [~, residual] = lyap(A, B, C);
    
            % Residual should be small for LTI case
            assert(residual.resnorm < 1e-6, ...
                sprintf('Sylvester residual: %e', residual.resnorm));
        end

        function testYalmipLmi(testCase)
            % TEST_YALMIP_LMI Checks LMI stability condition: P > 0, \dot{P} + A'P + PA < 0
            % For a simple constant stable system A = -I.
    
            nx = 2;
            h = 1;
            T = 2*pi; % Period
    
            % Stable system A = -eye(2)
            A = PhasorArray(-eye(nx));
    
            % Lyapunov Matrix P(t) as ndsdpvar
            P = PhasorArray.ndsdpvar(nx, nx, h);
    
            PT = P.T_tb(h);
            PA  = P*A;
            ATP = A'*P;
            Q   = PhasorArray(eye(nx));
            QT  = Q.T_tb(h);
    
            % Derivative operator matrix for harmonic domain
            % Corresponds to d/dt in frequency domain: N = diag(j*k*omega)
            N = N_tb(nx, h, T);
    
            % LMI: \dot{P} + A'P + PA <= -epsilon * I
            % \dot{P} term corresponds to commutator with N: N*PT - PT*N
            % Using the formula from Exemples/Exemple_Toolbox_LMI.m line 47:
            Sum = ATP + PA;
            LMI_lhs = Sum.T_tb(h) + QT + (N*PT - PT*N);
    
            % Constraints
            % P >= epsilon
            % LMI <= -epsilon
            epsilon = 1e-3;
    
            F = [PT >= epsilon*eye(size(PT)), ...
                 LMI_lhs <= -epsilon*eye(size(LMI_lhs))];
         
            % Solve with minimal output
            opts = sdpsettings('verbose', 0);
            sol = optimize(F, [], opts);
    
            if sol.problem ~= 0
               % If solver issue, warn
               warning('YALMIP optimize returned problem code %d: %s', sol.problem, sol.info);
            else
                % Verify P is positive definite
                P_val = value(PT);
                min_eig_P = min(real(eig(P_val)));

                % Verify LMI is negative definite
                LMI_val = value(LMI_lhs);
                max_eig_LMI = max(real(eig(LMI_val)));

                % Assertions
                testCase.verifyTrue(min_eig_P >= epsilon - 1e-6, 'P should be positive definite');
                testCase.verifyTrue(max_eig_LMI <= -epsilon + 1e-6, 'LMI should be negative definite');
            end
        end

        function testYalmipLmiVsLyap(testCase)
            % TEST_YALMIP_LMI_VS_LYAP Compares X from lyap(A,Q) with X from LMI.
            % Solves A'X + XA + dX/dt + Q = 0.
    
            % System definition (from test_lyap_LTP)
            A0 = [-3 0.5; 0 -2];
            A1 = [0.1 0; 0 0.05];
            A2 = [0 0.3; 0.1 0];
            % A(t) has frequencies -2, -1, 0, 1, 2. So h=2.
            A = PhasorArray(cat(3, conj(A2), conj(A1), A0, A1, A2));
            Q = PhasorArray(eye(2));
    
            h_sol = 10;
            nx = 2;
            T = 2*pi; % Assume default period T=2*pi (omega=1)
    
            % 1. Solve using lyap (Arithmetic)
            % Note: lyap(A, Q, ...) uses defaults if T not specified.
            [X_lyap, ~] = lyap(A, Q, 'h', h_sol,'autoUpdateh',true,'T',T);
    
            % 2. Solve using YALMIP LMI (Optimization)
            % Equation: A'X + XA + dX/dt + Q = 0
            % We solve this as a feasibility problem with equality constraints.
    
            X_yalmip = PhasorArray.ndsdpvar(nx, nx, h_sol);
    
            XT = X_yalmip.T_tb(h_sol);
            ATX = A'*X_yalmip;
            XA = X_yalmip*A;
            QT = Q.T_tb(h_sol);
    
            % Derivative operators
            N = N_tb(nx, h_sol, T);
    
            % dX/dt term -> N*XT - XT*N
            % Equation (A'X + XA + Q + dX/dt).T_tb = 0
    
            Sum = ATX + XA;
            Lyap_Op = Sum.T_tb(h_sol) + QT + (N*XT - XT*N);
    
            % Constraint: Lyap_Op <= 0 (Relaxed equality to inequality)
            % We minimize trace(XT) to find the 'smallest' solution that satisfies 
            % the inequality (which converges to the equality solution for stable A).
            F = [Lyap_Op <= 0; XT>=0];
    
            % Solve: minimize trace(XT)
            opts = sdpsettings('verbose', 0);
            sol = optimize(F, trace(XT), opts);
    
            if sol.problem ~= 0
                warning('YALMIP optimize returned problem code %d: %s', sol.problem, sol.info);
            else
                % Compare solutions
                X_val = value(XT);
                X_ref = X_lyap.T_tb(h_sol);
        
                diff = max(abs(X_val - X_ref), [], 'all');
        
                % Normalize error by norm of solution
                rel_err = diff / max(abs(X_ref), [], 'all');
        
                % Assertion
                % Consistency should be high if T matches.
                testCase.verifyTrue(rel_err < 1e-5, sprintf('Mismatch between lyap() and LMI solution. Rel Err: %e', rel_err));
            end
        end


        function testLyapSolvesBothIntegrationDirections(testCase)
            % backward integrates P(t) = int_t^inf Phi' Q Phi, forward
            % P(t) = int_-inf^t Phi Q Phi'. Both are exposed and they are not the
            % same solution: each satisfies its own equation and fails the other.
            hh = testCase.h;
            N  = N_tb(size(testCase.A, 1), hh, testCase.T);
            fwdSign = [1 -1];
            names = {'backward', 'forward'};

            P = cell(1, 2);
            for d = 1:2
                P{d} = lyap(testCase.A, testCase.Q, 'h', hh, 'T', testCase.T, ...
                    'autoUpdateh', true, 'direction', names{d});
            end
            testCase.verifyGreaterThan(max(abs(P{1}.value - P{2}.value), [], 'all'), 1e-6, ...
                'the two directions returned the same solution');

            for d = 1:2
                XT = P{d}.T_tb(hh);
                own   = norm(T_tb(testCase.A'*P{d} + P{d}*testCase.A, hh) ...
                             + testCase.Q.T_tb(hh) + fwdSign(d)*(N*XT - XT*N), 'fro');
                other = norm(T_tb(testCase.A'*P{d} + P{d}*testCase.A, hh) ...
                             + testCase.Q.T_tb(hh) - fwdSign(d)*(N*XT - XT*N), 'fro');
                % Truncation at hh, not the solver's own accuracy: the two
                % residuals separate by five orders, which is the point.
                testCase.verifyLessThan(own, 1e-5, ...
                    sprintf('%s does not satisfy its own equation', names{d}));
                testCase.verifyGreaterThan(other, 1e-3, ...
                    sprintf('%s also satisfies the opposite equation', names{d}));
            end
        end


        function testPlaceRefinesTheOrderOnAStrongPerturbation(testCase)
            % At the fixed hLyap = 4*A.h a strongly periodic system places its
            % exponents to 5e-05 only, and says so through the Sylvester
            % residual. Refinement is on by default so a caller who never
            % thought about the order still gets the poles asked for.
            A0 = [0 1; -2 -3];
            A1 = [0.6 0.3; 0.2 0.5];
            A  = PhasorArray(cat(3, conj(A1), A0, A1));
            B  = PhasorArray([0; 1]);
            mus = [-4; -5];
            hv = 80;

            placedError = @(K) max(arrayfun(@(m) ...
                min(abs(eig(subsref(A - B*K, substruct('.','T_tb','()',{hv})) ...
                             - N_tb(2, hv, 2*pi)) - m)), mus));

            [Kfix, ~, ~, iFix] = place(A, B, mus, 'autoUpdateh', false);
            [Kad,  ~, ~, iAd]  = place(A, B, mus);

            testCase.verifyEqual(iFix.status, 3, 'fixed h should report status 3');
            testCase.verifyEqual(iAd.status, 0, ...
                sprintf('refinement did not converge: %s', iAd.statusMsg));
            testCase.verifyGreaterThan(iAd.h, iFix.h, 'refinement should raise the order');

            errFix = placedError(Kfix);
            errAd  = placedError(Kad);
            testCase.verifyLessThan(errAd, 1e-10, ...
                sprintf('refined placement is off by %.2e', errAd));
            testCase.verifyLessThan(errAd, errFix / 100, ...
                'refinement should beat the fixed order by orders of magnitude');
            testCase.verifyLessThan(iAd.resrelnorm, 1e-10, 'residual not converged');
        end


        function testEverySolverRefinesByDefault(testCase)
            % The five entry points share one driver and must share one default,
            % otherwise the same problem converges or not depending on which one
            % the caller reached for. mrHmcDivide used to default to a fixed
            % order while its mirror mlHmcDivide refined.
            Ad = 0.2*PhasorArray.random(2, 2, 1) + PhasorArray.eye(2);
            Bd = PhasorArray.random(2, 1, 1);
            E  = 0.1*PhasorArray.random(2, 2, 1) + PhasorArray.eye(2);

            [~, iL]  = lyap(testCase.A, testCase.Q);
            [~, iG]  = lyapG(testCase.A, testCase.Q, E);
            [~, iMl] = mlHmcDivide(Ad, Bd);
            [~, iMr] = mrHmcDivide(PhasorArray.random(1, 2, 1), Ad);
            [~, ~, ~, iP] = place(testCase.A, PhasorArray([1; 0.5]), [-6; -7]);

            names = {'lyap', 'lyapG', 'mlHmcDivide', 'mrHmcDivide', 'place'};
            infos = {iL, iG, iMl, iMr, iP};
            for k = 1:numel(infos)
                testCase.verifyNotEqual(infos{k}.status, 3, ...
                    sprintf('%s did not refine with its default settings', names{k}));
                testCase.verifyEqual(infos{k}.status, 0, ...
                    sprintf('%s did not converge: %s', names{k}, infos{k}.statusMsg));
            end
        end

        function testOperatorFormPointsAtTheMethodWhenItFails(testCase)
            % \ and / take no name-value pairs, so a caller who hits a
            % non-converged solve needs to be told where the knobs live.
            % Independent of the ambient warning state: a suite run with
            % warnings off would otherwise see nothing to catch.
            prev = warning('on', 'PhasorArray:divide:notConverged');
            restore = onCleanup(@() warning(prev));

            Ad = 0.2*PhasorArray.random(2, 2, 1) + PhasorArray.eye(2);
            Bd = PhasorArray.random(2, 1, 1);
            testCase.verifyWarning( ...
                @() mldivide(Ad, Bd, 'maxh', 1, 'thresholdResidual', 1e-14), ...
                'PhasorArray:divide:notConverged');
            % A converged solve stays silent.
            testCase.verifyWarningFree(@() mldivide(Ad, Bd));
        end


        function testRelativeResidualIsScaleInvariant(testCase)
            % Every solver divides the residual by the right-hand side, never by
            % the solution. That makes resrelnorm invariant when the whole
            % problem is scaled: multiplying Q by 1000 multiplies P, and the
            % residual, by the same 1000. A solver normalising by anything else
            % -- the solution, or nothing at all -- would move by three orders.
            k = 1000;

            [~, i1] = lyap(testCase.A, testCase.Q);
            [~, i2] = lyap(testCase.A, k * testCase.Q);
            testCase.verifyEqual(i2.resrelnorm, i1.resrelnorm, 'RelTol', 1e-6, ...
                'lyap: relative residual moved with the scale of Q');

            Ad = 0.2*PhasorArray.random(2, 2, 1) + PhasorArray.eye(2);
            Bd = PhasorArray.random(2, 1, 1);
            [~, i3] = mlHmcDivide(Ad, Bd);
            [~, i4] = mlHmcDivide(Ad, k * Bd);
            testCase.verifyEqual(i4.resrelnorm, i3.resrelnorm, 'RelTol', 1e-6, ...
                'mlHmcDivide: relative residual moved with the scale of B');

            % And the absolute one must move, otherwise nothing was scaled.
            testCase.verifyEqual(i2.resnorm / i1.resnorm, k, 'RelTol', 1e-3, ...
                'lyap: the absolute residual did not follow the scaling');
        end

        function testEnergyAndFrobeniusAgreeOnThePayload(testCase)
            % The solvers write the same quantity two ways -- norm(X.value,'fro')
            % in lyap and mlHmcDivide, sqrt(energy(X)) in place. They must be one
            % metric, which by Parseval is the L2 norm of X(t) over a period.
            A = PhasorArray.random(3, 3, 4);
            v = A.value;
            testCase.verifyEqual(sqrt(energy(A)), norm(v(:)), 'RelTol', 1e-12);
            testCase.verifyEqual(sqrt(energy(A)), norm(v, 'fro'), 'RelTol', 1e-12);
        end

    end

    methods (Test, TestTags = {'Install'})
        % Smoke set for a fresh install: no optional toolbox, a few seconds,
        % one check per layer. Run with run_all_tests("install").
        function testLyapLti(testCase)
            % LTI Lyapunov: A'P + PA + Q = 0 with constant A, Q
            A = PhasorArray([-2 1; 0 -3]);
            Q = PhasorArray(eye(2));
    
            [X, residual] = lyap(A, Q);
    
            % Residual should be small
            assert(residual.resnorm < 1e-6, ...
                sprintf('LTI Lyapunov residual: %e', residual.resnorm));
    
            % The DC part should match MATLAB's built-in lyap
            n = size(A, 1);
            A_const = [-2 1; 0 -3];
            X_dc = X{1:n, 1:n, 0};
            X_ref = lyap(A_const', eye(n));
            err = max(abs(X_dc - X_ref), [], 'all');
            testCase.verifyTrue(err < 1e-6, sprintf('LTI Lyapunov DC vs MATLAB lyap: %e', err));
        end

        function testLyapSolvesFastWithSmallResidual(testCase)
            tic;
            [X, info] = lyap(testCase.A, testCase.Q, ...
                'h', testCase.h, 'autoUpdateh', true, 'T', testCase.T);
            elapsed = toc;

            testCase.verifyLessThan(info.resnorm, 1e-5, ...
                sprintf('lyap reports a large residual: %.2e', info.resnorm));
            testCase.verifyLessThan(testCase.residual(X), 1e-5, ...
                'the recomputed residual disagrees with the reported one');
            testCase.verifyEqual(info.status, 0, ...
                sprintf('lyap did not converge, status %d', info.status));
            % Measured at 0.05 s; the bound only catches a collapse in the
            % adaptive-order rule, not ordinary machine-to-machine spread.
            testCase.verifyLessThan(elapsed, 5, ...
                sprintf('lyap took %.2f s on a 2x2 h=10 problem', elapsed));
        end
    end
end
