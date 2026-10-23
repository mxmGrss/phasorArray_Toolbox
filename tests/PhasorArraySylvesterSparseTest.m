classdef PhasorArraySylvesterSparseTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function sourcePath(testCase)
            source = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Fonctions');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(source,'IncludingSubfolders',true));
            rehash; clear SylvHarmonic
            testCase.assertTrue(startsWith(which('SylvHarmonic'),source));
            state = rng; testCase.addTeardown(@()rng(state)); rng(20261004);
        end
    end
    methods (Test)
        function manufacturedBroadbandSolutions(testCase)
            for temporalReal = [true false]
                for bandwidth = [3 20]
                    A = testCase.coefficients(3,3,bandwidth,temporalReal);
                    B = testCase.coefficients(2,2,bandwidth-1,temporalReal);
                    A(:,:,bandwidth+1) = -20*eye(3); B(:,:,bandwidth) = -15*eye(2);
                    X = testCase.coefficients(3,2,2,temporalReal); h = bandwidth+2; omega = 3.7;
                    for direction = {'backward','forward'}
                        C = testCase.forcing(A,B,X,omega,direction{1});
                        for method = {'square','rectangle'}
                            actual = SylvHarmonic(A,B,C,h,omega,method{1},direction{1});
                            expected = zeros(3,2,2*h+1);
                            expected(:,:,h+1+(-2:2)) = X;
                            testCase.verifyEqual(actual,expected,'AbsTol',1e-11);
                        end
                    end
                end
            end
        end
        function overdeterminedOutputContract(testCase)
            A = testCase.coefficients(2,2,3,false); A(:,:,4) = -10*eye(2);
            B = testCase.coefficients(3,3,2,false); B(:,:,3) = -9*eye(3);
            C = testCase.coefficients(2,3,12,false); h = 4; omega = 4.2;
            for direction = {'backward','forward'}
                [X,M,L,R,q,x] = SylvHarmonic(A,B,C,h,omega,'rectangle',direction{1});
                testCase.verifyEqual(size(M),[6*25,6*9]);
                testCase.verifyTrue(issparse(M) && issparse(L) && issparse(R));
                testCase.verifyEqual(x,reshape(permute(X,[3 1 2]),[],1),'AbsTol',1e-12);
                testCase.verifyEqual(q,reshape(permute(C,[3 1 2]),[],1),'AbsTol',1e-12);
                [AX,XB] = testCase.products(A,B,X,12);
                testCase.verifyEqual(L*x,reshape(permute(AX,[3 1 2]),[],1),'AbsTol',1e-11);
                testCase.verifyEqual(R*x,reshape(permute(XB,[3 1 2]),[],1),'AbsTol',1e-11);
                signN = -1; if strcmp(direction{1},'forward'), signN = 1; end
                testCase.verifyLessThan(norm(M+L+R-signN*N_tb(6,[12 h],[],omega=omega),'fro'),1e-12);
                testCase.verifyLessThan(norm(M'*(M*x-q))/max(norm(M,'fro')*norm(q),eps),1e-12);
                testCase.verifyEqual(SylvHarmonic(A,B,C,h,omega,'rectangle',direction{1}),X,'AbsTol',1e-11);
            end
        end
        function dcAndPhasorInputs(testCase)
            A = -eye(2); B = -2*eye(3); C = ones(2,3);
            testCase.verifyEqual(SylvHarmonic(A,B,C,0,0),C/3,'AbsTol',1e-12);
            testCase.verifyEqual(SylvHarmonic(PhasorArray(A),PhasorArray(B),PhasorArray(C),0,0),C/3,'AbsTol',1e-12);
            raw = cat(3,3i*ones(2),A,-2i*ones(2));
            testCase.verifyEqual(SylvHarmonic(raw,B,C,0,0,'square'),C/3,'AbsTol',1e-12);
        end
        function invalidDimensions(testCase)
            testCase.verifyError(@()SylvHarmonic(ones(2,3),eye(2),eye(2),2,1),...
                'PhasorArray:SylvHarmonic:dimensions');
            testCase.verifyError(@()SylvHarmonic(ones(2,2,2),eye(2),eye(2),2,1),...
                'PhasorArray:SylvHarmonic:dimensions');
        end
    end
    methods (Static)
        function Z = coefficients(n,m,d,temporalReal)
            Z = (randn(n,m,2*d+1)+1i*randn(n,m,2*d+1))/(4*sqrt(max(d,1)));
            if temporalReal
                Z(:,:,d+1) = real(Z(:,:,d+1));
                Z(:,:,1:d) = conj(flip(Z(:,:,d+2:end),3));
            end
        end
        function [AX,XB] = products(A,B,X,hOut)
            hA = (size(A,3)-1)/2; hB = (size(B,3)-1)/2; hX = (size(X,3)-1)/2;
            AX = zeros(size(X,1),size(X,2),2*hOut+1); XB = AX;
            for k = -hOut:hOut
                for j = -hX:hX
                    if abs(k-j) <= hA
                        AX(:,:,hOut+1+k) = AX(:,:,hOut+1+k)+A(:,:,hA+1+k-j)*X(:,:,hX+1+j);
                    end
                    if abs(k-j) <= hB
                        XB(:,:,hOut+1+k) = XB(:,:,hOut+1+k)+X(:,:,hX+1+j)*B(:,:,hB+1+k-j);
                    end
                end
            end
        end
        function C = forcing(A,B,X,omega,direction)
            hX = (size(X,3)-1)/2; hC = hX+max((size(A,3)-1)/2,(size(B,3)-1)/2);
            [AX,XB] = PhasorArraySylvesterSparseTest.products(A,B,X,hC);
            C = -AX-XB;
            signN = -1; if strcmp(direction,'forward'), signN = 1; end
            for k = -hX:hX
                C(:,:,hC+1+k) = C(:,:,hC+1+k)+signN*1i*k*omega*X(:,:,hX+1+k);
            end
        end
    end
end
