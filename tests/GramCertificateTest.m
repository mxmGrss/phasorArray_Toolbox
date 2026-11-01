classdef GramCertificateTest < matlab.unittest.TestCase
    % Exact Gram certificates of periodic matrix inequalities (gramConstraint)
    % and their a posteriori validation (checkGramCertificate).
    methods (TestClassSetup)
        function sourcePath(testCase)
            source = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(source,'IncludingSubfolders',true));
            testCase.assumeTrue(exist('sdpvar','file') == 2, 'YALMIP required');
        end
    end
    methods (Test)
        function exactBoundOfAScalarPolynomial(testCase)
            % m(t) = 2 + cos(t) has minimum 1: margin 0.9 is feasible, 1.1 is not.
            M = PhasorArray(reshape([0.5 2 0.5],1,1,3));
            for basis = ["hermitian","real"]
                [F,G,c] = M.gramConstraint(0.9,basis);
                testCase.verifyEqual(GramCertificateTest.solve(F),0,basis);
                r = M.checkGramCertificate(0.9,value(G),c);
                testCase.verifyTrue(r.accepted,basis);
                testCase.verifyGreaterThan(r.lower,0.9-1e-6,basis);
                F = M.gramConstraint(1.1,basis);
                testCase.verifyNotEqual(GramCertificateTest.solve(F),0,basis);
            end
        end
        function maximalMarginIsTheMinimumOfTheSymbol(testCase)
            % 2x2 real symmetric, odd degree 3: half-integer real basis.
            rng(7); C = randn(2,2,7)/6; C = (C + conj(flip(C,3)))/2;
            C = (C + permute(C,[2 1 3]))/2; C(:,:,4) = C(:,:,4) + 2*eye(2);
            M = PhasorArray(C);
            samples = evalp(M,linspace(0,2*pi,4001));
            trueMin = min(arrayfun(@(k) min(eig(samples(:,:,k))),1:size(samples,3)));
            for basis = ["hermitian","real"]
                mu = sdpvar(1);
                F = M.gramConstraint(mu,basis);
                testCase.verifyEqual(GramCertificateTest.solve(F,-mu),0,basis);
                testCase.verifyEqual(value(mu),trueMin,'AbsTol',1e-5,basis);
            end
        end
        function sectionAcceptsWhatGramRejects(testCase)
            % m(t) = 1 + 1.2 cos(t) is negative near t = pi, yet its 1x1 section
            % (the mean) is positive: sections are necessary conditions only.
            M = PhasorArray(reshape([0.6 1 0.6],1,1,3));
            testCase.verifyGreaterThan(min(eig(M.T_tb(0))),0);
            testCase.verifyNotEqual(GramCertificateTest.solve(M.gramConstraint(0)),0);
        end
        function lyapunovLmiReachesTheLyapunovSolution(testCase)
            % a = -1, q(t) = 1 + cos(t)/2: the solution of dP/dt - 2P + q = 0
            % is the degree-1 polynomial p_k = q_k/(2 - jk); min p_0 over the
            % LMI dP/dt + A'P + PA + Q <= 0 is attained there.
            A = PhasorArray(-1); Q = PhasorArray(reshape([0.25 1 0.25],1,1,3));
            P = PhasorArray.ndsdpvar(1,1,1);
            L = -(P.d(2*pi) + A.'*P + P*A + Q);
            [F1,G1,c1] = L.gramConstraint(0,"real");
            F2 = P.gramConstraint(1e-3,"real");
            testCase.verifyEqual(GramCertificateTest.solve([F1,F2],phas(P,0)),0);
            testCase.verifyEqual(value(phas(P,0)),0.5,'AbsTol',1e-6);
            r = sdpval(L).checkGramCertificate(0,value(G1),c1);
            testCase.verifyGreaterThan(r.lower,-1e-6);
        end
        function nonHermitianDataRejected(testCase)
            M = PhasorArray(reshape([1 2 0],1,1,3));
            testCase.verifyError(@() M.gramConstraint(0),'PhasorArray:gramConstraint:notHermitian');
            testCase.verifyError(@() PhasorArray(ones(2,3)).gramConstraint(0), ...
                'PhasorArray:gramConstraint:notSquare');
            % A scalar Hermitian symbol is always real; a complex one needs n > 1.
            testCase.verifyError(@() PhasorArray([2 1i; -1i 2]).gramConstraint(0,"real"), ...
                'PhasorArray:gramConstraint:notReal');
        end
    end
    methods (Static, Access = private)
        function problem = solve(F,objective)
            if nargin < 2, objective = []; end
            sol = optimize(F,objective,sdpsettings('verbose',0));
            problem = sol.problem;
        end
    end
end
