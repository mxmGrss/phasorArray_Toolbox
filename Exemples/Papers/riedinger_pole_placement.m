% Riedinger-Daafouz, Harmonic Pole Placement, corrected formulation, Sec.IV-A, Eq.21.
% Requires Symbolic Math Toolbox (FloquetDec). Outputs: K, P, info, poles.
%[text] # Riedinger: harmonic pole placement
clear;
[A,B,T] = riedingerVernereyPlant(20,3,0.5);
[V,J] = FloquetDec(A,T,FixedStepTransPow=12,FixedStepInitProbPow=12);
V = trunc(PhasorArray(V),40); % explicit Floquet approximation order
G = B'/V'; Lambda = -J'-eye(2);
[K,P,res,info] = place(A,B,Lambda,G=G,T=T,hLyap=20,maxh=200,thresholdResidual=1e-8);
poles = HmqNEig(A-B*K,80,T,'fundamental');
[~,ip] = sort(imag(poles)); target = eig(Lambda); [~,it] = sort(imag(target));
poleError = max(abs(poles(ip)-target(it)));
assert(info.resrelnorm<1e-8 && poleError<1e-3,'Paper:Placement','Placement check failed.');
disp(struct('poles',poles,'poleError',poleError))
%[appendix]{"version":"1.0"}
