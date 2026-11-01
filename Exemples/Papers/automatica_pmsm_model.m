% Build a variable-speed PMSM and its output-integrator augmentation.
%[text] # Variable-frequency PMSM model
clear;
r = 0.5; L = 1.5e-3; psi = 0.14; J = 0.03; Bf = 0.02; p = 4;
omega = [10 200]; % mechanical rad/s; phase theta is mechanical
phi = [PhasorArray.sin(0,p);PhasorArray.sin(-2*pi/3,p);PhasorArray.sin(2*pi/3,p)];
A = [-r/L*eye(3),psi*p/L*phi;-p*psi/J*phi',-Bf/J];
Bu = PhasorArray([eye(3)/L;zeros(1,3)]);
E = diag([L L L J]);
assert(sqrt(energy(E*A+A'*E+2*diag([r r r Bf])))<1e-12,'Paper:PMSM','Energy balance failed.');
Tp = PhasorArray.Park(0,p,false);
Cnegative = [zeros(1,3),1;Tp{1,:},0]; Cpositive = Cnegative;
for k = [0 2 6 8]
    if k==0, Tk = PhasorArray([2/3 -1/3 -1/3;0 1/sqrt(3) -1/sqrt(3)]); else, Tk = PhasorArray.Park(0,k,false); end
    Cnegative = [Cnegative;Tk,zeros(2,1)]; Cpositive = [Cpositive;diag([1 -1])*Tk,zeros(2,1)]; %#ok<AGROW> Four fixed harmonic blocks.
end
vertexNegative = cell(1,2); vertexPositive = cell(1,2);
for j = 1:2
    vertexNegative{j} = [A,zeros(4,10);omega(j)*Cnegative,zeros(10)];
    vertexPositive{j} = [A,zeros(4,10);omega(j)*Cpositive,zeros(10)];
end
Baug = [Bu;zeros(10,3)]; Qfull = eye(14); Qintegrators = diag([zeros(1,4),ones(1,10)]); R = 100*eye(3);
assert(isequal(size(A,[1 2]),[4 4]) && isequal(size(Cnegative,[1 2]),[10 4]),'Paper:PMSM','Unexpected dimensions.');
disp(struct('plantOrder',h(A),'outputOrder',h(Cnegative),'states',size(vertexNegative{1},1),'omegaMechanical',omega))
% At vertex omega(j): operator = T(vertex{j})-omega(j)*N, N uses period 2*pi.
%[appendix]{"version":"1.0"}
