function periodic_defective_demo()
% PERIODIC_DEFECTIVE_DEMO  Periodic system with a prescribed defective Floquet
%   exponent, built by the Floquet-Lyapunov construction
%       A(t) = Rdot(t) R(t)^-1 + R(t) J R(t)^-1,
%   with R(t) a periodic rotation and J a constant Jordan block. The Floquet
%   exponents equal eig(J) (defective), yet A(t) is a genuine full periodic
%   matrix, so HmqNEig produces real truncation pollution -- unlike a constant
%   matrix, whose Hill operator is exactly block-diagonal and trivially clean.
%
%   Period note: omega must match PhasorArray's fundamental (default 2*pi here).
%   A mismatch scatters the whole HSS spectrum (difficulty D7).

omega = 1; Tper = 2*pi/omega;
J = [-1 1; 0 -1];                 % Jordan block: lambda=-1, alg mult 2, geom mult 1
a = 0.6;                          % rotation amplitude
Nt = 512; t = (0:Nt-1)/Nt*Tper;

At = zeros(2,2,Nt);
for i = 1:Nt
    th  = a*sin(omega*t(i));  thd = a*omega*cos(omega*t(i));
    R   = [cos(th) -sin(th); sin(th) cos(th)];
    Rd  = thd*[-sin(th) -cos(th); cos(th) -sin(th)];
    At(:,:,i) = Rd/R + R*J/R;
end
assert(max(abs(squeeze(At(2,1,:)))) > 0.1, 'construction collapsed to triangular');

Aper = PhasorArray(TimeArray2Phasors(At, 1, t));

fprintf('periodic defective (lambda=-1, mult 2):\n');
for h = [20 40 80]
    E  = HmqNEig(Aper, h); E = E(:);
    mu = findTruelmV3(E, Tper, 2);
    assert(max(abs(mu+1)) < 1e-4, 'h=%d: recovered %s', h, mat2str(mu.'));
    fprintf('  h=%2d | mu = %s | max|mu+1| = %.1e\n', h, mat2str(round(sort(mu),4).'), max(abs(mu+1)));
end
fprintf('PASS\n');
end
