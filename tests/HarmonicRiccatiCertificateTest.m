classdef HarmonicRiccatiCertificateTest < matlab.unittest.TestCase
    % A posteriori certification of harmonic Riccati candidates (certifyHare):
    % Proposition 7 and Theorem 8 of the corrected Lyapunov/Riccati paper.
    methods (TestClassSetup)
        function sourcePath(testCase)
            source = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(source,'IncludingSubfolders',true));
        end
    end
    methods (Test)
        function exactScalarLti(testCase)
            % x' = -x + u, Q = R = 1: P+ = sqrt(2) - 1.
            c = certifyHare(PhasorArray(-1),PhasorArray(1),PhasorArray(1),1, ...
                PhasorArray(sqrt(2)-1),h=0);
            testCase.verifyTrue(c.accepted);
            testCase.verifyLessThan(c.errorBound,1e-12);
            testCase.verifyEqual(c.beta,1/(2*sqrt(2)),'AbsTol',1e-12);
        end
        function perturbedCandidateBoundsKnownError(testCase)
            c = certifyHare(PhasorArray(-1),PhasorArray(1),PhasorArray(1),1, ...
                PhasorArray(sqrt(2)-1+0.01),h=0);
            testCase.verifyTrue(c.accepted);
            testCase.verifyGreaterThanOrEqual(c.errorBound+1e-13,0.01);
        end
        function unstableClosedLoopRejected(testCase)
            c = certifyHare(PhasorArray(1),PhasorArray(1),PhasorArray(1),1, ...
                PhasorArray(0),h=0);
            testCase.verifyFalse(c.accepted);
            testCase.verifyFalse(c.stabilityEstablished);
        end
        function stableButPoorCandidateRejected(testCase)
            c = certifyHare(PhasorArray(-1),PhasorArray(1),PhasorArray(1),1, ...
                PhasorArray(0),h=0);
            testCase.verifyTrue(c.stabilityEstablished);
            testCase.verifyFalse(c.accepted);   % theta = 1 is not sufficient
            testCase.verifyEqual(c.status,"small_gain_failed");
        end
        function manufacturedPeriodicSolution(testCase)
            % Prescribe P(t) > 0, derive Q exactly, then perturb P.
            A = PhasorArray(-2); B = PhasorArray(1);
            P = PhasorArray(reshape([0.05 1 0.05],1,1,3));
            Q = -P.d(2*pi)-A.'*P-P*A+P*B*B.'*P;
            dP = PhasorArray(reshape([0.0005 0 0.0005],1,1,3));
            c = certifyHare(A,B,Q,1,P+dP,h=12);
            testCase.verifyTrue(c.accepted);
            testCase.verifyGreaterThanOrEqual(c.errorBound+1e-12,0.001);
            testCase.verifyGreaterThanOrEqual(c.residualOrder,2);
        end
        function nonphysicalDataRejected(testCase)
            A = PhasorArray(-1); B = PhasorArray(1); Q = PhasorArray(1);
            testCase.verifyError(@() certifyHare(A,B,Q,-1,Q),'PhasorArray:certifyHare:R');
            testCase.verifyError(@() certifyHare(A,B,Q,1,PhasorArray(1i)), ...
                'PhasorArray:certifyHare:real');
        end
        function negativeBetweenGridPointsIsRejected(testCase)
            % Q(t) dips below zero between the samples of a coarse grid.
            for phase = [pi/256,-pi/256]
                Q = PhasorArray(reshape([-exp(1i*phase)/2,1-1e-6,-exp(-1i*phase)/2],1,1,3));
                sampled = real(evalp(Q,(0:255)*2*pi/256));
                testCase.verifyGreaterThan(min(sampled(:)),0);
                testCase.verifyError(@() certifyHare(PhasorArray(-1),PhasorArray(1), ...
                    Q,1,PhasorArray(0.4),h=4,points=256),'PhasorArray:certifyHare:Q');
            end
        end
        function hareCertifiesItsOwnCandidate(testCase)
            % a(t) = -1 + cos(t)/2, b(t) = 1 + cos(t)/5, Q = R = 1.
            A = PhasorArray(reshape([0.25 -1 0.25],1,1,3));
            B = PhasorArray(reshape([0.1 1 0.1],1,1,3));
            [~,X,info] = hare(A,B,1,1,certify=true);
            testCase.verifyTrue(info.certificate.accepted);
            testCase.verifyLessThan(info.certificate.errorBound,1e-9);
            F = X.d(2*pi)+A.'*X+X*A-X*(B*B.')*X+1;
            testCase.verifyEqual(info.certificate.residualOrder,h(F));
        end
        function hareRejectsCertifyOnDescriptorProblem(testCase)
            testCase.verifyError(@() hare(PhasorArray(-1),PhasorArray(1),1,1, ...
                E=PhasorArray(1),certify=true),'PhasorArray:hare:certifyUnsupported');
        end
    end
end
