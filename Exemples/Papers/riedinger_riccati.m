% Compute periodic LQR feedback u=-K*x and certify the Riccati solution.
%[text] # Riedinger: least-squares Riccati and independent certificate
clear;
p = 15; hG = 120; D = diag([1 64]);
[A,B,T] = riedingerVernereyPlant(p,3);
Q = PhasorArray(diag([1 1e4])); R = 1;
[K,P,info] = hare(A,B,Q,R,T=T,h=hG,autoUpdateh=false,thresholdResidual=1e-12);
certificate = certifyHare(D*A/D,D*B,D'\Q/D,R,D'\P/D,T=T,h=2*hG+30);
assert(certificate.accepted,'Paper:Riccati','Certificate rejected: %s',certificate.status);
disp(struct('h',info.h,'meanTrace',real(trace(phas(P,0))),'certificate',certificate))
%[appendix]{"version":"1.0"}
