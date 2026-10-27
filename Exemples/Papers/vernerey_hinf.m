% Bound the disturbance-to-output L2 gain with feedback u=K*x.
% Requires YALMIP/MOSEK. Outputs: P, Y, K (u=K*x), certificates.
%[text] # Vernerey: H-infinity synthesis with real Gram constraints
clear;
p = 15; q = 20; gamma = 24.922; margin = 1e-3; D = diag([1 10]);
[A,B,T] = riedingerVernereyPlant(p,3); A = D*A/D; B = D*B;
Bw = PhasorArray(D); Cz = PhasorArray([1 0;0 100;0 0]/D); Dzu = PhasorArray([0;0;1]);
P = PhasorArray.ndsdpvar(2,2,q); Y = PhasorArray.ndsdpvar(1,2,q,'symmetry','real');
L = A*P+P*A'-P.d(T)+B*Y+Y'*B'; H = Cz*P+Dzu*Y;
M = -[L,Bw,H'; Bw',-gamma*eye(2),zeros(2,3); H,zeros(3,2),-gamma*eye(3)];
[F1,G1,c1] = M.gramConstraint(margin,"real"); [F2,G2,c2] = P.gramConstraint(margin,"real");
sol = optimize([F1,F2],-trace(D\phas(P,0)/D'),sdpsettings('solver','mosek','verbose',0,'dualize',1));
assert(sol.problem==0,'Paper:Hinf','SDP failed: %s',sol.info);
certificates = [sdpval(M).checkGramCertificate(margin,value(G1),c1), sdpval(P).checkGramCertificate(margin,value(G2),c2)];
assert(all([certificates.accepted]),'Paper:Hinf','Gram validation failed; increase margin.');
P = sdpval(P); Y = sdpval(Y);
K = @(theta) evalp(Y,theta)/evalp(P,theta)*D; % exact rational gain, u=K(theta)*x
disp(struct('gamma',gamma,'lowerBounds',[certificates.lower],'gainAtZero',K(0)))
%[appendix]{"version":"1.0"}
