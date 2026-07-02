function [Aper, info] = randomWithNPole(J_or_V, h, opts)
% RANDOMWITHNPOLE  Generate a random periodic matrix with prescribed Floquet exponents.
%
%   Returns a PhasorArray representing a random nx×nx time-periodic matrix A(t)
%   whose Floquet exponents are the eigenvalues of J_or_V.
%
%   Floquet exponents play the same role as eigenvalues for LTI systems:
%   stability (Re < 0 → stable), oscillatory modes (Im part).
%
%   Two underlying constructions, three options (opts.Method):
%
%   --- Construction 1 : triangular P  (Method = 'structured') ---------------
%
%     P(t) = I + Q(t),  Q strictly lower-triangular, K = floor(h/2) harmonics
%     per entry.  P^{-1} computed exactly via Neumann series (Q^{nx} = 0).
%     A(t) = Q̇·P⁻¹ + P·J·P⁻¹,  then trunc(h).
%     Everything stays in the harmonic domain — no time sampling, no aliasing.
%
%     Result: Floquet exponents = eig(J) to floating-point precision.
%             Bandwidth = exactly h.
%             A(t) has lower-triangular structure (full generic coupling only
%             when J is non-diagonal, e.g. Jordan blocks).
%
%   --- Construction 2 : rotation expm  (Method = 'generic' | 'truncated') ---
%
%     G = random unit skew-symmetric matrix.
%     θ(t) = Σ_{k=1}^{hRot} amp/k · sin(k·t + φk)
%     R(t) = expm(θ(t)·G)   sampled on Nt points, then FFT → PhasorArray.
%     A(t) = θ̇·G + R·J·R⁻¹  (computed sample by sample).
%
%     R(t) has infinite Fourier harmonics (Jacobi-Anger expansion).  The FFT
%     aliases harmonics beyond Nt/2; aliasing error ≈ J_{Nt/2}(amp) < 1e-8
%     for default amp.  A(t) is fully generic (no structural bias).
%
%     'generic'   — keep all FFT harmonics (~Nt/2).  Use when the wide-band
%                   content is intentional (convergence studies, stress tests).
%     'truncated' — apply trunc(h) after FFT.  Bandwidth = exactly h.
%                   Adds a Bessel truncation error |J_{h+1}(amp)| on top of
%                   the aliasing error; < 1e-11 for amp≤0.4, h≥8.
%                   Use when you need controlled bandwidth and generic coupling.
%
%   Summary:
%
%     Method         Floquet error         h exact   A(t) generic
%     'structured'   ~0 for real/well-      yes          no
%                     separated modes;
%                     up to ~1e-4 for
%                     complex modes close
%                     in (Re,Im) -- solver
%                     extraction limit,
%                     see findTruelm
%     'generic'      < 1e-8 (alias)         no           yes
%     'truncated'    < 1e-11 (Bessel)       yes          yes
%
%   Which to pick: 'structured' (default) covers almost every case. Its only
%   drawback is a triangular/sparse A(t) coupling -- fix that with
%   Densify=true: post-conjugation by a random constant orthogonal S
%   densifies A(t) to full generic coupling while leaving the Floquet
%   exponents EXACTLY unchanged (constant similarity, zero added error) and
%   bandwidth still exactly h. structured+Densify strictly dominates
%   'truncated' (same guarantees, zero error instead of <1e-11) -- prefer it.
%   Reserve 'generic'/'truncated' for the one case they still cover:
%   intentional wide-band content beyond h (testing a solver's behaviour
%   under real aliasing/truncation, not just missing high harmonics).
%
% Syntax
%   Aper = PhasorArray.randomWithNPole(V, h)
%   Aper = PhasorArray.randomWithNPole(J, h)
%   Aper = PhasorArray.randomWithNPole(J_or_V, h, Method='generic')
%   Aper = PhasorArray.randomWithNPole(J_or_V, h, seed=42)
%
% Inputs
%   J_or_V : prescribed Floquet exponents — either:
%              nx×1 real/complex vector  (diagonal J built automatically)
%              nx×nx real matrix         (use Jordan blocks for defective
%                                         eigenvalues, companion blocks for
%                                         complex conjugate pairs)
%              Example for exponents -1 ± 0.4i : J = [-1 0.4; -0.4 -1]
%              nx is inferred from the size of J_or_V.
%   h      : target harmonic order
%
% Options (name-value)
%   Method    'structured' (default) | 'generic' | 'truncated'
%   Densify   true|false (default false). Post-conjugate by a random constant
%             orthogonal matrix S: A -> S*A*S'. Densifies the coupling to
%             full generic structure at ZERO cost in Floquet accuracy and
%             bandwidth (constant similarity). S and the pre-conjugation
%             matrix L are returned in info.
%   seed      RNG seed for reproducibility (default: unseeded)
%   amp       Rotation amplitude, 'generic'/'truncated' only (default 0.4).
%             Larger amp → richer harmonic content, larger Bessel error.
%   hRot      Harmonics in rotation angle, 'generic'/'truncated' only (default 4).
%
% Output
%   Aper : PhasorArray, nx×nx, omega = 1 (T = 2π)
%   info : struct with the full chain of similarity-transform matrices:
%          J -> Q -> P -> L -> (S) -> Aper
%            .J    prescribed Jordan/diagonal form (nx×nx double)
%            .Q    'structured' only: Q(t), PhasorArray
%            .P    'structured' only: P(t) = I+Q(t), PhasorArray
%            .Pinv 'structured' only: P(t)^{-1}, PhasorArray (exact)
%            .L    Densify=true only: pre-conjugation matrix
%                  L(t) = d(P)*Pinv + P*J*Pinv (exponents already = eig(J))
%            .S    Densify=true only: constant orthogonal mixing matrix,
%                  Aper = S*L*S'
%          Rebuild check: Aper == trunc(d(P)*Pinv + P*J*Pinv, h) (Densify
%          off) -- exact to machine precision.
%
% Known limitation
%   Exponents with Im(μ) = ±ω/2 lie on the Brillouin zone boundary and
%   cause aliasing ambiguity in Floquet extraction.  Keep |Im(μ)| < 0.45.
%
% Examples
%   % Stable system, two real exponents
%   Ap = PhasorArray.randomWithNPole([-1; -0.5], 20);
%
%   % Defective eigenvalue  (algebraic mult 2, geometric mult 1)
%   Ap = PhasorArray.randomWithNPole([-1 1; 0 -1], 20);
%
%   % Complex conjugate pair  -1 ± 0.4i
%   Ap = PhasorArray.randomWithNPole([-1 0.4; -0.4 -1], 20);
%
%   % Generic coupling, spectrally rich
%   Ap = PhasorArray.randomWithNPole([-1; -0.5], 20, Method='generic', hRot=6);
%
%   % Controlled bandwidth, generic structure, approx Floquet (Bessel < 1e-11)
%   Ap = PhasorArray.randomWithNPole([-1; -0.5], 20, Method='truncated');
%
%   % Dense coupling, EXACT Floquet exponents, exact bandwidth h (preferred
%   % over 'truncated' -- see "Which to pick" above)
%   [Ap, info] = PhasorArray.randomWithNPole([-1; -0.5], 20, Densify=true);
%   % info.S = mixing matrix, info.L = pre-conjugation matrix
%
% See also: findTruelm, findTrueFloquet, HmqNEig

arguments
    J_or_V  double
    h       (1,1) {mustBeInteger,mustBePositive}
    opts.Method (1,1) string ...
        {mustBeMember(opts.Method, ["structured","generic","truncated"])} = "structured"
    opts.seed   = []
    opts.amp    (1,1) double = 0.4
    opts.hRot   (1,1) double = 4
    opts.acAmp  (1,1) double = 0.3   % 'structured' only: Q harmonic amplitude (AC/DC ratio knob)
    opts.Densify (1,1) logical = false  % conjugate by a random constant orthogonal S (exact)
end

if ~isempty(opts.seed), rng(opts.seed); end

% --- Build J and infer nx
if isvector(J_or_V)
    nx = numel(J_or_V);
    J  = diag(J_or_V(:));
else
    assert(ismatrix(J_or_V) && size(J_or_V,1)==size(J_or_V,2), ...
        'J_or_V must be a vector (nx×1) or square matrix (nx×nx)');
    nx = size(J_or_V, 1);
    J  = J_or_V;
end

% --- Dispatch
info = struct('J', J, 'P', [], 'Pinv', [], 'Q', [], 'L', [], 'S', []);
switch opts.Method
    case "structured"
        [Aper, info.P, info.Pinv, info.Q] = buildStructured(nx, J, h, opts.acAmp);
    case "generic"
        Aper = buildGeneric(nx, J, h, opts.amp, opts.hRot);
    case "truncated"
        Aper = buildGeneric(nx, J, h, opts.amp, opts.hRot);
        Aper = Aper.trunc(h);
end

% --- Optional densification: constant orthogonal similarity (exact) ---------
if opts.Densify
    info.L = Aper;                 % pre-conjugation matrix, exponents = eig(J)
    info.S = orth(randn(nx));
    Aper   = info.S * Aper * info.S';
end
end

% =========================================================================
function [Aper, P_pa, Pinv, Q_pa] = buildStructured(nx, J, h, acAmp)
% Floquet-Lyapunov with lower-triangular unipotent P = I + Q.
% P^{-1} exact via Neumann series (Q^nx = 0).
% Floquet exponents = eig(J) exact.  Output bandwidth = exactly h.

T  = 2*pi;
K  = max(1, floor(h/2));   % harmonics per Q entry
Nt = 2^nextpow2(max(64, 8*K));
t  = (0:Nt-1)/Nt * T;

% Q(t): strictly lower triangular, exactly K harmonics per entry
Q_t = zeros(nx, nx, Nt);
for i = 2:nx
    for j = 1:i-1
        amp    = acAmp * randn(K,1) ./ (1:K)';
        phases = 2*pi * rand(K,1);
        Q_t(i,j,:) = sum(amp .* cos((1:K)' .* t + phases), 1);
    end
end

Q_pa = PhasorArray(TimeArray2Phasors(Q_t, 1, t));
I_pa = PhasorArray(eye(nx));
J_pa = PhasorArray(J);
P_pa = I_pa + Q_pa;

% P^{-1} = sum_{k=0}^{nx-1} (-Q)^k  (exact, Q^nx = 0)
Pinv = I_pa;
Qpow = Q_pa;
for k = 1:nx-1
    if mod(k,2)==1, Pinv = Pinv - Qpow; else, Pinv = Pinv + Qpow; end
    if k < nx-1,    Qpow = Qpow * Q_pa; end
end

A_pa = d(Q_pa) * Pinv  +  P_pa * J_pa * Pinv;
Aper = A_pa.trunc(h);
end

% =========================================================================
function Aper = buildGeneric(nx, J, h, amp, hRot)
% Floquet-Lyapunov with R(t) = expm(theta(t)*G), G random skew-symmetric.
% Floquet exponents = eig(J) exact.
% Output bandwidth ≈ Nt/2 (Jacobi-Anger expansion, not controlled to h).

T     = 2*pi;
omega = 1;
Nt    = 2^nextpow2(8 * (2*h + 1));
t     = (0:Nt-1)/Nt * T;

% Random unit skew-symmetric generator
S = randn(nx);
G = (S - S') / 2;
G = G / norm(G,'fro');

% Multi-harmonic rotation angle: theta(t) = sum amp/k * sin(k*omega*t + phi)
amps   = amp ./ (1:hRot);
phases = 2*pi * rand(1, hRot);

At = zeros(nx, nx, Nt);
for i = 1:Nt
    kvec = 1:hRot;
    th   = sum(amps .* sin(kvec * omega * t(i) + phases));
    thd  = sum(amps .* (kvec * omega) .* cos(kvec * omega * t(i) + phases));
    R    = expm(th * G);
    At(:,:,i) = thd * G  +  R * J / R;
end

is_real_At = isreal(At) || (max(abs(imag(At(:)))) < 1e-12);
Aper = PhasorArray(TimeArray2Phasors(At, 1, t, 'isReal', is_real_At));
end
