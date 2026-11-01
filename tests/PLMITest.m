classdef PLMITest < matlab.unittest.TestCase
    % Periodic matrix inequality collections (pLMI): construction, Gram and
    % section compilation, validation reports and the Lyapunov/Riccati factories.
    methods (TestClassSetup)
        function sourcePath(testCase)
            source = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(source,'IncludingSubfolders',true));
            testCase.assumeTrue(exist('sdpvar','file') == 2, 'YALMIP required');
        end
    end
    methods (Test)
        function scalarConvention(t)
            X = sdpvar(2,2,'symmetric'); P = PhasorArray(X);
            C = P >= 1; D = P >= eye(2); E = 1 <= P;
            assign(X,eye(2));
            t.verifyEqual(double(value(C.Symbols{1})),eye(2)-ones(2),'AbsTol',1e-12);
            t.verifyEqual(double(value(D.Symbols{1})),zeros(2),'AbsTol',1e-12);
            t.verifyEqual(double(value(E.Symbols{1})),eye(2)-ones(2),'AbsTol',1e-12);
        end
        function collection(t)
            P = PhasorArray.ndsdpvar(2,2,1);
            C = [P >= eye(2), P <= 3*eye(2)];
            t.verifyEqual(C.count(),2);
            D = [C;P >= 0]; t.verifyEqual(D.count(),3);
        end
        function sharedVariablesAcrossVelocities(t)
            P = PhasorArray.ndsdpvar(2,2,1);
            delta = value(P.d()-P.d(2*pi));
            t.verifyLessThan(norm(getbase(reshape(delta,[],1)),'fro'),1e-12);
            A = PhasorArray(-eye(2));
            C = [P.d()+A'*P+P*A <= -eye(2), 2*P.d()+A'*P+P*A <= -eye(2)];
            ids = getvariables(reshape(value(P),[],1));
            for k = 1:2
                t.verifyEqual(getvariables(reshape(value(C.Symbols{k}),[],1)),ids);
            end
            omega = sdpvar(1);
            t.verifyError(@() pLMI(PhasorArray(omega)*P.d()+A'*P+P*A),'pLMI:notAffine');
        end
        function rejectMalformed(t)
            X = sdpvar(2,2,'full');
            t.verifyError(@() pLMI(PhasorArray(X)),'pLMI:notHermitian');
            t.verifyError(@() pLMI(PhasorArray(ones(2,3))),'pLMI:notSquare');
            x = sdpvar(1);
            t.verifyError(@() pLMI(PhasorArray(x*x)),'pLMI:notAffine');
            P = PhasorArray.ndsdpvar(2,2,0);
            t.verifyError(@() (P > 0),'pLMI:strictInequality');
            t.verifyError(@() (P < 0),'pLMI:strictInequality');
            t.verifyError(@() (P >= ones(3)),'pLMI:dimensions');
        end
        function numericComparisonsAreUnchanged(t)
            P = PhasorArray([1 2;3 4]);
            [actual,fraction] = ge(P,2);
            t.verifyClass(actual,'logical');
            t.verifyEqual(actual,[false true;true true]);
            t.verifySize(fraction,[2 2]);
        end
        function lyapunovAtTwoVelocities(t)
            % One P certifies the same plant at angular velocities 1 and 2.
            A = PhasorArray([-1 1;-1 -2])+0.2*PhasorArray.cos()*[0 1;1 0];
            P = PhasorArray.ndsdpvar(2,2,1); mu = 0.1;
            C = [P >= mu*eye(2), P.d()+A'*P+P*A <= -mu*eye(2), 2*P.d()+A'*P+P*A <= -mu*eye(2)];
            [F,g] = C.gram("real");
            sol = optimize(F,[],sdpsettings('verbose',0));
            r = g.check(sol);
            t.verifyEqual(sol.problem,0);
            t.verifyTrue(all([r.strictlyCertified]));
        end
        function complexHermitian(t)
            x = sdpvar(1); C = pLMI(PhasorArray(x*eye(2)+[0 1i;-1i 0]));
            t.verifyError(@() C.gram("real"),'pLMI:notReal');
            [F,g] = C.gram();
            sol = optimize([F,x == 3],[],sdpsettings('verbose',0));
            r = g.check(sol);
            t.verifyEqual(sol.problem,0);
            t.verifyTrue(r.strictlyCertified);
            t.verifyEqual(r.lowerBound,2,'AbsTol',1e-5);
        end
        function sectionIsOnlyNecessary(t)
            % x + 2cos(theta) with x = 1: feasible order-0 section, negative at pi.
            x = sdpvar(1); M = PhasorArray(x)+2*PhasorArray.cos(); C = pLMI(M);
            [F,s] = C.section(0);
            sol = optimize([F,x == 1],[],sdpsettings('verbose',0));
            r = s.check(sol);
            t.verifyEqual(sol.problem,0); t.verifyTrue(r.sectionFeasible);
            t.verifyFalse(r.strictlyCertified);
            t.verifyEqual(r.status,"section_only");
            t.verifyEqual(sdpval(M).evalp(pi),-1,'AbsTol',1e-7);
            [F,g] = C.gram("real");
            sol = optimize([F,x == 1],[],sdpsettings('verbose',0));
            r = g.check(sol); t.verifyNotEqual(sol.problem,0);
            t.verifyEqual(r.status,"solver_failed");
        end
        function sectionOfHigherOrder(t)
            x = sdpvar(1); C = pLMI(PhasorArray(x)+2*PhasorArray.cos());
            [F,s] = C.section(2);
            t.verifySize(s.Matrices{1},[5 5]);
            sol = optimize([F,x == 3],[],sdpsettings('verbose',0));
            r = s.check(sol); t.verifyEqual(sol.problem,0);
            t.verifyTrue(r.sectionFeasible); t.verifyFalse(r.strictlyCertified);
            t.verifyEqual(r.sectionMinimum,3-sqrt(3),'AbsTol',1e-6);
        end
        function requestedMargin(t)
            x = sdpvar(1); P = PhasorArray(x*eye(2)); C = P >= eye(2);
            [~,g] = C.gram();
            assign(x,0.5); assign(g.Matrices{1},zeros(2));
            r = g.check(struct('problem',0));
            t.verifyFalse(r.strictlyCertified);
            t.verifyLessThanOrEqual(r.lowerBound,-0.49);
            assign(x,1); r = g.check(struct('problem',0));
            t.verifyFalse(r.strictlyCertified);
            assign(x,2); assign(g.Matrices{1},eye(2));
            r = g.check(struct('problem',0)); t.verifyTrue(r.strictlyCertified);
        end
        function diagnostics(t)
            P = PhasorArray.ndsdpvar(2,2,0); C = P >= eye(2);
            [~,g] = C.gram();
            r = g.check(struct('problem',1)); t.verifyEqual(r.status,"solver_failed");
            r = g.check(struct('problem',0)); t.verifyEqual(r.status,"nonfinite");
            t.verifyError(@() g.check(struct()),'pLMI:diagnostics');
        end
        function constantSection(t)
            C = pLMI(eye(2)); [F,s] = C.section(1);
            t.verifyEmpty(F); r = s.check(struct('problem',0));
            t.verifyTrue(r.sectionFeasible); t.verifyFalse(r.strictlyCertified);
            C = pLMI(-eye(2));
            t.verifyError(@() C.section(0),'pLMI:constantSectionInfeasible');
        end
        function symbolicScalarMargin(t)
            P = PhasorArray.ndsdpvar(2,2,1); mu = sdpvar(1);
            C = P >= mu;
            residual = phas(C.Symbols{1},0)-(phas(P,0)-mu*ones(2));
            t.verifyLessThan(norm(getbase(residual(:)),'fro'),1e-12);
            residual = phas(C.Symbols{1},1)-phas(P,1);
            t.verifyLessThan(norm(getbase(residual(:)),'fro'),1e-12);
        end
        function lyapunovFactory(t)
            [A,Q] = PLMITest.rotationProblem(1.7);
            [C,P] = pLMI.lyapunov(A,Q,1.7,2,0.1);
            manual = [P >= 0.1*eye(2),1.7*P.d()+A'*P+P*A <= -Q];
            for j = 1:2
                delta = value(C.Symbols{j}-manual.Symbols{j});
                t.verifyLessThan(norm(getbase(reshape(delta,[],1)),'fro'),1e-12);
            end
        end
        function invalidFactoryData(t)
            t.verifyError(@() pLMI.lyapunov(sdpvar(2),eye(2),1,2),'pLMI:numericSystem');
            t.verifyError(@() pLMI.lyapunov(ones(2,3),eye(2),1,2),'pLMI:systemDimensions');
            t.verifyError(@() pLMI.lyapunov(1i*eye(2),eye(2),1,2),'pLMI:realSystem');
            t.verifyError(@() pLMI.riccati(-eye(2),eye(2),eye(2),-eye(2),1,2),'pLMI:riccatiR');
        end
        function riccatiFactoryRecoversTheExactSolution(t)
            [A,B,Q,R,Pexact] = PLMITest.rotatingRiccatiProblem(1.3);
            residual = value(1.3*Pexact.d()+A'*Pexact+Pexact*A+Q-Pexact*B*(R\eye(2))*B'*Pexact);
            t.verifyLessThan(norm(residual(:)),1e-12);
            [C,P] = pLMI.riccati(A,B,Q,R,1.3,2,0.1);
            F = C.gram("real");
            sol = optimize(F,-trace(phas(P,0)),sdpsettings('verbose',0));
            t.verifyEqual(sol.problem,0);
            delta = value(sdpval(P)-Pexact);
            t.verifyLessThan(norm(delta(:)),1e-4);
        end
    end
    methods (Static, Access = private)
        function [A,Q,P] = rotationProblem(omega)
            % omega*P.d() + A'*P + P*A = -Q with A = omega*J + R(t)*diag(-1/2,-1)*R(t)',
            % P = R(t)*diag(2,1)*R(t)', Q = 2I, R(t) the rotation of angle t.
            C = PhasorArray.cos(0,2); S = PhasorArray.sin(0,2);
            Z = [1 0;0 -1]; X = [0 1;1 0]; J = [0 -1;1 0];
            A = PhasorArray(omega*J-0.75*eye(2))+0.25*C*Z+0.25*S*X;
            P = PhasorArray(1.5*eye(2))+0.5*C*Z+0.5*S*X;
            Q = PhasorArray(2*eye(2));
        end
        function [A,B,Q,R,P] = rotatingRiccatiProblem(omega)
            % Same rotating plant, B = R = I, Q = diag(3,8) in rotating
            % coordinates: the stabilizing solution is diagonal there.
            A = PLMITest.rotationProblem(omega);
            B = PhasorArray(eye(2)); R = eye(2);
            roots = [-0.5+sqrt(0.25+3),-1+sqrt(1+8)];
            Z = [1 0;0 -1]; X = [0 1;1 0];
            C = PhasorArray.cos(0,2); S = PhasorArray.sin(0,2);
            P = PhasorArray(mean(roots)*eye(2))+diff(fliplr(roots))/2*(C*Z+S*X);
            Q = PhasorArray(5.5*eye(2))-2.5*(C*Z+S*X);
        end
    end
end
