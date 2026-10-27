function [residual,xe,se,Ae] = blinRectifierResidual(coefficients)
%BLINRECTIFIERRESIDUAL Weighted equilibrium problem, corrected paper Sec.6.1.
% INPUT: DC, real and imaginary positive Fourier coefficients of s, h=5.
% OUTPUT: real weighted residual, equilibrium, control and state matrix.
% EXAMPLE: residual = blinRectifierResidual([0;zeros(5,1);-0.25;zeros(4,1)]);
arguments
    coefficients (11,1) double
end
R = 1e-3; L = 1e-3; C = 4e-3; RL = 10; T = 1/50;
se = ScalarPhasorArray([coefficients(1);coefficients(2:6)+1i*coefficients(7:11)].','isreal',true);
Ae = [-R/L,-se/L;se/C,-1/(RL*C)];
forcing = [100*PhasorArray.sin()/L;0];
xe = lyap(-Ae,0,-forcing,T=T,h=5,systemType="square");
current = reshape(phas(xe{1,1},-5:5),[],1);
voltage = reshape(phas(xe{2,1},-5:5),[],1);
v = [sqrt(1e3)*(voltage(6)-200);voltage([1:5,7:11]);sqrt(1e3)*current([1:4,8:11]);real(current)];
residual = [real(v);imag(v)];
end
