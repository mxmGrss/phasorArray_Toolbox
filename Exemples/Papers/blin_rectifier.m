% Blin et al., corrected formulation, Sec.6.1-6.2, Eqs.49 and 35.
% Requires Optimization Toolbox. Outputs: xe, se, P and feedback law s(t,x).
%[text] # Blin: single-phase rectifier equilibrium and Lyapunov feedback
clear;
T = 1/50; gamma = 1e-4; h = 5;
c0 = [0;zeros(h,1);-0.25;zeros(h-1,1)];
[c,~,~,exitflag] = lsqnonlin(@blinRectifierResidual,c0,[],[],optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',4000));
assert(exitflag>0,'Paper:Rectifier','Equilibrium optimization did not converge.');
[objectiveResidual,xe,se,Ae] = blinRectifierResidual(c);
[P,info] = lyap(Ae,diag([1 0.01]),T=T,h=20,autoUpdateh=true,maxh=100);
g = [0 -1e3;250 0];
feedback = @(t,x) min(1,max(-1,real(evalp(se,2*pi*t/T)-gamma*(g*x)'*evalp(P,2*pi*t/T)*(x-evalp(xe,2*pi*t/T)))));
[t,x] = ode15s(@(t,x) (diag([-1 -25])+feedback(t,x)*g)*x+[1e5*sin(2*pi*t/T);0],[0 0.2],[0;0],odeset('RelTol',1e-7,'AbsTol',1e-8));
equilibriumResidual = sqrt(energy(xe.d(T)-Ae*xe-[1e5*PhasorArray.sin();0]));
assert(info.resrelnorm<1e-6,'Paper:Rectifier','Lyapunov residual too large.');
disp(struct('meanVoltage',real(phas(xe{2,1},0)),'objective',norm(objectiveResidual)^2,'equilibriumResidual',equilibriumResidual,'lyapunovResidual',info.resrelnorm,'startupFinalVoltage',x(end,2)))
% h=5 square equilibrium reproduces the paper truncation; report its full residual.
%[appendix]{"version":"1.0"}
