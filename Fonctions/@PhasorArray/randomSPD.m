function P = randomSPD(nx, h, nvp)
% RANDOMSPD  Generate a random symmetric positive definite periodic matrix.
%
%   P(t) = L(t)·L(t)ᵀ  where L(t) is lower-triangular with positive diagonal.
%   Guarantees P(t) ≻ 0 for all t by construction: x'P(t)x = ‖L(t)'x‖² > 0.
%
%   P(t) is generic (no diagonal dominance): off-diagonal entries of L have
%   a nonzero DC offset plus oscillations, giving P significant off-diagonal
%   energy at all frequencies including DC.
%
% Syntax
%   P = PhasorArray.randomSPD(nx, h)
%   P = PhasorArray.randomSPD(nx, h, seed=42, amp=0.5)
%
% Inputs
%   nx  : matrix dimension
%   h   : harmonic order of output PhasorArray
%
% Options
%   seed   : RNG seed (default: unseeded)
%   amp    : relative oscillation amplitude on diagonal entries (default 0.4).
%            Diagonal of L: c_i·(1 + amp·f(t)) with |f| < 1 → always positive.
%            Off-diagonal amplitude scales with amp too.
%   alphaMin : prescribed lower bound, P(t) ⪰ alphaMin·I for all t (default 0).
%            Adding alphaMin·I to L(t)L(t)' gives the bound pointwise and exactly,
%            since L(t)L(t)' ⪰ 0 already. It is a pure DC term, so the bandwidth
%            and the lossless truncation below are unaffected.
%
%            Note this is a *pointwise* guarantee, which is the strong one.
%            T_h(P) ⪰ alphaMin·I would not be: a finite Toeplitz block is the
%            compression of the infinite operator, so lambda_min(T_h(P)) only
%            decreases towards min_t lambda_min(P(t)) as h grows. The equivalence
%            P(t) ⪰ alpha·I  <=>  T(P) ⪰ alpha·I holds at h = inf alone.
%
% Output
%   P : PhasorArray, nx×nx, symmetric PD, bandwidth = h
%
% Period
%   None to give: P(t) = L(t)L(t)' is positive definite for whatever period
%   the caller later evaluates at, since positivity is a property of the
%   coefficients and no derivative enters the construction. Measured identical
%   at T = 2*pi, 1 and 0.1. Contrast randomWithNPole, whose rotation methods
%   do differentiate and therefore take a T.

arguments
    nx   (1,1) {mustBeInteger,mustBePositive}
    h    (1,1) {mustBeInteger,mustBePositive}
    nvp.seed = []
    nvp.amp  (1,1) double = 0.4
    nvp.alphaMin (1,1) double {mustBeNonnegative} = 0
end

assert(nvp.amp < 1, 'amp must be < 1 to guarantee positive diagonal');
if ~isempty(nvp.seed), rng(nvp.seed); end

T  = 2*pi;
K  = max(1, floor(h/2));
Nt = 2^nextpow2(max(64, 8*K));
t  = (0:Nt-1)/Nt * T;

% --- Build L(t) sample by sample
L_t = zeros(nx, nx, Nt);
for i = 1:nx
    % Diagonal: c_i > 0, oscillations bounded by amp*c_i → always positive
    c_i    = 0.5 + rand();                        % base in [0.5, 1.5]
    amps_d = nvp.amp * c_i * rand(K,1) ./ (1:K)';  % |sum| < amp*c_i
    phi_d  = 2*pi * rand(K,1);
    L_t(i,i,:) = c_i * (1 + sum(amps_d./c_i .* cos((1:K)'.*t + phi_d), 1));

    % Off-diagonal (strictly lower triangular): DC offset + oscillations
    % DC offset gives P_{ij} a nonzero mean → significant off-diagonal energy
    for j = 1:i-1
        c_off  = nvp.amp * c_i * randn();          % nonzero constant
        amps_o = nvp.amp * c_i * randn(K,1) ./ (1:K)';
        phi_o  = 2*pi * rand(K,1);
        L_t(i,j,:) = c_off + sum(amps_o .* cos((1:K)'.*t + phi_o), 1);
    end
end

% --- P(t) = L(t)*L(t)' + alphaMin*I sample by sample
P_t = zeros(nx, nx, Nt);
aI  = nvp.alphaMin * eye(nx);
for k = 1:Nt
    Lk       = L_t(:,:,k);
    P_t(:,:,k) = Lk * Lk' + aI;
end

P = PhasorArray(TimeArray2Phasors(P_t, 1, t));
P = P.trunc(h);
end
