% Solve a scalar periodic Lyapunov equation and measure its full residual.
%[text] # Riedinger: scalar least-squares Lyapunov
clear;
T = 2*pi; orders = [2 4 6 8 10 12 14 16 18 20 24];
A = -1-PhasorArray.cos()+2*PhasorArray.sin()+PhasorArray.cos(0,2);
residualNorm = zeros(size(orders));
for j = 1:numel(orders)
    P = lyap(A,1,T=T,h=orders(j),autoUpdateh=false);
    residualNorm(j) = sqrt(energy(P.d(T)+A'*P+P*A+1));
end
disp(table(orders.',residualNorm.','VariableNames',{'h','FullResidual'}))
assert(residualNorm(end)<1e-9,'Paper:Lyapunov','Full residual exceeds tolerance.');
%[appendix]{"version":"1.0"}
