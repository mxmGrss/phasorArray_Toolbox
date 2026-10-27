% Riedinger-Daafouz, Lyapunov/Riccati, corrected formulation, Sec.7, Eq.35 and Table 2.
% No workspace inputs. Outputs: K (u=-K*x), P, info, certificate.
%[text] # Riedinger: least-squares Riccati and independent certificate
clear;
p = 15; hG = 120; D = diag([1 64]);
[A,B,T] = riedingerVernereyPlant(p,3);
Q = PhasorArray(diag([1 1e4])); R = 1;
[K,P,info] = hare(A,B,Q,R,T=T,h=hG,autoUpdateh=false,thresholdResidual=1e-12);
certificate = certifyHare(D*A/D,D*B,D'\Q/D,R,D'\P/D,T=T,h=2*hG+30);
assert(certificate.accepted,'Paper:Riccati','Certificate rejected: %s',certificate.status);
disp(struct('h',info.h,'meanTrace',real(trace(phas(P,0))),'certificate',certificate))
% Bounds concern p=15 in coordinates z=D*x, not the original infinite series.
%[appendix]{"version":"1.0"}
