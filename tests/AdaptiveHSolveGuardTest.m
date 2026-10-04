classdef AdaptiveHSolveGuardTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function sourcePath(testCase)
            source = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(source,'IncludingSubfolders',true));
            rehash; clear adaptiveHSolve RicHarmonicKlein
            testCase.assertTrue(startsWith(which('adaptiveHSolve'),source));
        end
    end
    methods (Test)
        function plateauBelowBandwidthCanConverge(testCase)
            cfg = testCase.configuration();
            solve = @(h)deal(h,10^-max(0,h-20),10^-max(0,h-20),[]);
            [best,trace] = adaptiveHSolve(solve,2,cfg);
            testCase.verifyEqual(trace.status,0);
            testCase.verifyEqual(best.h,26);
        end
        function wholeWindowMustResolveBandwidth(testCase)
            cfg = testCase.configuration();
            stuck = @(h)deal(h,1,1,[]);
            [~,trace] = adaptiveHSolve(stuck,2,cfg);
            testCase.verifyEqual(trace.status,1);
            testCase.verifyEqual(trace.h_history(end-cfg.stagnationWindow+1),22);
        end
        function orderCapBelowBandwidthRemainsEffective(testCase)
            cfg = testCase.configuration(); cfg.maxh = 10;
            stuck = @(h)deal(h,1,1,[]);
            [~,trace] = adaptiveHSolve(stuck,2,cfg);
            testCase.verifyEqual(trace.status,2);
            testCase.verifyEqual(trace.h_history(end),10);
        end
        function forcedExtrapolationDoesNotEnablePrematureStopping(testCase)
            cfg = testCase.configuration(); cfg.updateMethod = 'adaptive'; cfg.maxUnitSteps = 3;
            solve = @(h)deal(h,10^-max(0,h-20),10^-max(0,h-20),[]);
            [best,trace] = adaptiveHSolve(solve,2,cfg);
            testCase.verifyEqual(trace.status,0);
            testCase.verifyLessThanOrEqual(best.resrelnorm,cfg.thresholdResidual);
        end
        function preBandPowerLawIsNotUnreachableTargetEvidence(testCase)
            cfg = testCase.configuration(); cfg.updateMethod = 'adaptive';
            cfg.stagnationWindow = 100; cfg.stagnationRatio = 0;
            cfg.hOp = 80; cfg.maxh = 160;
            residual = @(h)(h+1)^(-.5)*10^-max(0,h-80);
            solve = @(h)deal(h,residual(h),residual(h),[]);
            [best,trace] = adaptiveHSolve(solve,2,cfg);
            testCase.verifyEqual(trace.status,0);
            testCase.verifyLessThanOrEqual(best.resrelnorm,cfg.thresholdResidual);
        end
        function frozenKleinmanIterateIsNotReportedAsConverged(testCase)
            state = warning('off','RicHarmonicKlein:skipValidate');
            testCase.addTeardown(@()warning(state));
            A = PhasorArray(diag([-1.5 -2.5])); B = PhasorArray([1;.4]);
            [~,~,info] = RicHarmonicKlein(A,B,PhasorArray(diag([3 1])),PhasorArray(1),...
                PhasorArray(zeros(1,2)),2*pi,'h',2,'autoUpdateh',false,...
                'thresholdResidual',1e-14,'relChangeThreshold',1e6);
            testCase.verifyEqual(info.status,1);
            testCase.verifyGreaterThan(info.resRicnorm,1e-14);
            testCase.verifyTrue(contains(info.statusMsg,'frozen without residual convergence'));
        end
        function broadbandKleinmanConvergesWithDefaultGuards(testCase)
            p = 20; a = zeros(2,2,2*p+1); a(:,:,p+1) = [1 2;-1 1];
            for k = 1:p
                ak = zeros(2);
                if mod(k,2)==1, ak(1,1) = -2i/(pi*k); ak(1,2) = 8/(pi^2*k^2); end
                ak(2,1) = (-1)^k*exp(1i*pi/4)/(1i*pi*k);
                if k==1, ak(2,2)=1i; elseif k==3, ak(2,2)=1+1i; elseif k==5, ak(2,2)=1; end
                a(:,:,p+1+k) = ak; a(:,:,p+1-k) = conj(ak);
            end
            b = zeros(2,1,13); b(1,1,7)=1; b(1,1,[5 9])=1;
            b(1,1,1)=2i; b(1,1,13)=-2i;
            A = PhasorArray(a); B = PhasorArray(b); Q = PhasorArray(100*eye(2));
            [~,P,info] = RicHarmonicKlein(A,B,Q,PhasorArray(1),[],1,...
                'h',2,'autoUpdateh',true,'maxh',200,'maxIter',30,'thresholdResidual',1e-6);
            residual = d(P,1)+A.'*P+P*A-P*B*B.'*P+Q;
            testCase.verifyEqual(info.status,0);
            testCase.verifyLessThan(norm(value(residual),'fro')/norm(value(Q),'fro'),1e-6);
        end
    end
    methods (Static)
        function cfg = configuration()
            cfg = struct('thresholdResidual',1e-6,'maxh',40,...
                'stagnationWindow',5,'stagnationRatio',.05,'updateMethod','incremental',...
                'verbose',false,'hOp',20,'hOutFcn',@(h)h,'preamble','','label','h');
        end
    end
end
