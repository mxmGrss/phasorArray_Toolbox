function P = randomSPD(nx, h, opts)
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
%
% Output
%   P : PhasorArray, nx×nx, symmetric PD, omega=1 (T=2π), bandwidth = h

arguments
    nx   (1,1) {mustBeInteger,mustBePositive}
    h    (1,1) {mustBeInteger,mustBePositive}
    opts.seed = []
    opts.amp  (1,1) double = 0.4
end

assert(opts.amp < 1, 'amp must be < 1 to guarantee positive diagonal');
if ~isempty(opts.seed), rng(opts.seed); end

T  = 2*pi;
K  = max(1, floor(h/2));
Nt = 2^nextpow2(max(64, 8*K));
t  = (0:Nt-1)/Nt * T;

% --- Build L(t) sample by sample
L_t = zeros(nx, nx, Nt);
for i = 1:nx
    % Diagonal: c_i > 0, oscillations bounded by amp*c_i → always positive
    c_i    = 0.5 + rand();                        % base in [0.5, 1.5]
    amps_d = opts.amp * c_i * rand(K,1) ./ (1:K)';  % |sum| < amp*c_i
    phi_d  = 2*pi * rand(K,1);
    L_t(i,i,:) = c_i * (1 + sum(amps_d./c_i .* cos((1:K)'.*t + phi_d), 1));

    % Off-diagonal (strictly lower triangular): DC offset + oscillations
    % DC offset gives P_{ij} a nonzero mean → significant off-diagonal energy
    for j = 1:i-1
        c_off  = opts.amp * c_i * randn();          % nonzero constant
        amps_o = opts.amp * c_i * randn(K,1) ./ (1:K)';
        phi_o  = 2*pi * rand(K,1);
        L_t(i,j,:) = c_off + sum(amps_o .* cos((1:K)'.*t + phi_o), 1);
    end
end

% --- P(t) = L(t)*L(t)' sample by sample
P_t = zeros(nx, nx, Nt);
for k = 1:Nt
    Lk       = L_t(:,:,k);
    P_t(:,:,k) = Lk * Lk';
end

P = PhasorArray(TimeArray2Phasors(P_t, 1, t));
P = P.trunc(h);
end
