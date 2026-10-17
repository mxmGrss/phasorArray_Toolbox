classdef RiccatiReturnContractTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function sourcePath(testCase)
            source = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(source,'IncludingSubfolders',true));
            rehash; clear RicHarmonicKlein
            testCase.assertTrue(startsWith(which('RicHarmonicKlein'),source));
            state = warning('off','RicHarmonicKlein:skipValidate');
            testCase.addTeardown(@()warning(state));
        end
    end
    methods (Test)
        function residualExitReturnsUpdatedGain(testCase)
            [K,S,info] = RicHarmonicKlein(PhasorArray(-1),PhasorArray(1),1,1,0,2*pi,...
                'h',0,'autoUpdateh',false,'maxIter',1,'thresholdResidual',1);
            testCase.verifyEqual(info.status,0);
            testCase.verifyEqual(value(S),.5,'AbsTol',1e-13);
            testCase.verifyEqual(value(K),value(S),'AbsTol',1e-13);
        end
        function frozenExitReturnsUpdatedGain(testCase)
            [K,S,info] = RicHarmonicKlein(PhasorArray(-1),PhasorArray(1),1,1,0,2*pi,...
                'h',0,'autoUpdateh',false,'maxIter',1,...
                'thresholdResidual',1e-14,'relChangeThreshold',2);
            testCase.verifyEqual(info.status,1);
            testCase.verifyEqual(value(K),value(S),'AbsTol',1e-13);
        end
        function iterationLimitDiagnosticsMatchReturnedBest(testCase)
            testCase.checkBestReturn(3,100);
        end
        function stagnationDiagnosticsMatchReturnedBest(testCase)
            testCase.checkBestReturn(2,12);
        end
    end
    methods (Access=private)
        function checkBestReturn(testCase,expectedStatus,window)
            % An unresolved harmonic leaves a residual floor after scalar convergence.
            A = PhasorArray(-1,2,'isreal',true); B = PhasorArray(1);
            [K,S,info] = RicHarmonicKlein(A,B,1,1,0,2*pi,...
                'h',0,'autoUpdateh',false,'maxIter',20,...
                'thresholdResidual',0,'relChangeThreshold',0,...
                'stagnationWindow',window,'stagnationRatio',.05,'storeIterates',true);
            testCase.verifyEqual(info.status,expectedStatus);
            [~,best] = min(info.resRic_history);
            testCase.verifyEqual(value(S),value(info.S_history{best}));
            testCase.verifyEqual(value(K),value(info.K_history{best}));
            testCase.verifyEqual(value(K),value(S),'AbsTol',1e-13);
            residual = d(S,2*pi)+A.'*S+S*A-S*B*B.'*S+1;
            testCase.verifyEqual(info.resRicnorm,norm(value(residual),'fro')/(1+eps),'AbsTol',1e-13);
            testCase.verifyEqual(info.resRicnorm,info.resRic_history(best));
            testCase.verifyEqual(info.relChange,info.relChange_history(best));
            testCase.verifyEqual(info.h,info.h_history(best));
            testCase.verifyEqual(numel(info.resRic_history),info.iter);
        end
    end
end
