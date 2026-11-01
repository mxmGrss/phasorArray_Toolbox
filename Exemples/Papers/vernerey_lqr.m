% Solve a periodic Riccati inequality with real Gram constraints.
% Requires YALMIP/MOSEK. Outputs: P, K (u=-K*x), certificates, meanTrace.
%[text] # Vernerey: direct-P LQR with real Gram constraints
clear;
p = 15; q = 40; margin = 1e-4; D = diag([1 100]);
[A,B,T] = riedingerVernereyPlant(p,3); A = D*A/D; B = D*B;
Q = PhasorArray(D'\diag([1 1e4])/D); R = 1;
P = PhasorArray.ndsdpvar(2,2,q);
M = [P.d(T)+A'*P+P*A+Q, P*B; B'*P, R];
[F1,G1,c1] = M.gramConstraint(margin,"real");
[F2,G2,c2] = P.gramConstraint(margin,"real");
sol = optimize([F1,F2],-trace(D'*phas(P,0)*D)/1e4,sdpsettings('solver','mosek','verbose',0,'dualize',1));
assert(sol.problem==0,'Paper:LQR','SDP failed: %s',sol.info);
certificates = [sdpval(M).checkGramCertificate(margin,value(G1),c1), sdpval(P).checkGramCertificate(margin,value(G2),c2)];
assert(all([certificates.accepted]),'Paper:LQR','Gram validation failed; increase margin.');
P = D'*sdpval(P)*D; [~,B] = riedingerVernereyPlant(p,3); K = B'*P;
meanTrace = real(trace(phas(P,0))); disp(struct('q',q,'meanTrace',meanTrace,'lowerBounds',[certificates.lower]))
%[appendix]{"version":"1.0"}
