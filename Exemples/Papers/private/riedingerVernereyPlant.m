function [A, B, T] = riedingerVernereyPlant(p, inputHarmonic, triangleFactor)
%RIEDINGERVERNEREYPLANT Analytic Fourier coefficients of the paper plant.
%   [A,B,T] = riedingerVernereyPlant(p,inputHarmonic) retains |k| <= p in A.
%   INPUTS: p >= 5; inputHarmonic = 3 for the corrected paper examples.
%   OUTPUTS: real-valued PhasorArrays A (2x2), B (2x1), period T = 1.
%   EXAMPLE: [A,B,T] = riedingerVernereyPlant(20,3);
%   Corrected manuscripts: Lyapunov/Riccati Eq.35; TBLMI Eq.43.
%   Corrected papers use inputHarmonic=3 and no phase offset in a21.
%   triangleFactor=0.5 for Pole Placement Eq.21 (8/pi^2, not 16/pi^2).
arguments
    p (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(p,5)}
    inputHarmonic (1,1) double {mustBeMember(inputHarmonic,[3 6])}
    triangleFactor (1,1) double {mustBePositive} = 1
end
T = 1;
coefficients = zeros(2,2,2*p+1);
coefficients(:,:,p+1) = [1 2; -1 1];
for k = 1:p
    Ak = zeros(2);
    if mod(k,2) == 1
        Ak(1,1) = -2i/(pi*k);
        Ak(1,2) = triangleFactor*8/(pi*k)^2;
    end
    Ak(2,1) = (-1)^k/(1i*pi*k);
    if k == 1, Ak(2,2) = 1i; end
    if k == 3, Ak(2,2) = 1+1i; end
    if k == 5, Ak(2,2) = 1; end
    coefficients(:,:,p+1+k) = Ak;
    coefficients(:,:,p+1-k) = conj(Ak);
end
A = PhasorArray(coefficients);
B = [1 + 2*PhasorArray.cos(0,2) + 4*PhasorArray.sin(0,inputHarmonic); 0];
end
