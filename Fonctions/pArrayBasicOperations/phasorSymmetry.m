function [A, held] = phasorSymmetry(A, req)
%PHASORSYMMETRY  Project a PhasorArray onto a symmetry class.
%   B = PHASORSYMMETRY(A, REQ) projects A onto the set of PhasorArrays
%   satisfying every property named in the string array REQ.
%
%   [B, HELD] = PHASORSYMMETRY(...) also returns the properties B satisfies.
%   HELD contains REQ and is usually larger: the satisfied set is closed
%   under composition, so REAL and SYMMETRIC together force HERMITIAN.
%
%   Valid names (14, exhaustive):
%       symmetric      skewSymmetric        A(t).' = +-A(t)
%       real           imaginary            conj(A(t)) = +-A(t)
%       even           odd                  A(-t) = +-A(t)
%       hermitian      skewHermitian        A(t)' = +-A(t)
%       paraSymmetric  skewParaSymmetric    A(-t).' = +-A(t)
%       paraConjugate  skewParaConjugate    conj(A(-t)) = +-A(t)
%       paraHermitian  skewParaHermitian    A(-t)' = +-A(t)
%
%   These are every property built from transposition, conjugation and time
%   reversal; the three generate (Z/2)^3. A PhasorArray therefore satisfies
%   0, 1, 3 or 7 of them at once, never 2, 4, 5 or 6, and the 14 names
%   combine into 51 distinct classes.
%
%   REQ defaults to "real", the usual case for a physical periodic system.
%   Pass [], "" or [""] to lift it and obtain a full complex array: there is
%   no name for unconstrained, it is the empty request.
%
%   Contradictory requests error out; the only PhasorArray satisfying
%   SYMMETRIC and SKEWSYMMETRIC together is zero.
%
%   Example
%       A = PhasorArray.random(3, 3, 3, "symmetry", []);
%       [B, held] = phasorSymmetry(A, ["real" "symmetric"]);
%       % held = ["symmetric" "real" "hermitian"]
%
%   See also symmetryClosure, mreal, mimag, mtranspose, mconj, retro,
%   trretro, ctrretro, PhasorArray.random, PhasorArray.ndsdpvar.

arguments
    A   PhasorArray
    req string = "real"
end

[held, eps_] = symmetryClosure(req);
if any(bitget(find(~isnan(eps_)) - 1, 1)) && size(A,1) ~= size(A,2)
    error('PhasorArray:symmetry:notSquare', ...
        'Requested symmetry transposes and needs a square array; got %dx%d.', ...
        size(A,1), size(A,2));
end

V  = pvalue(A);
el = find(~isnan(eps_));
P  = V;                                       % identity term, eps = +1 always
for j = el(2:end)
    P = P + eps_(j) * act(V, j-1);
end
P = P * (1/numel(el));                        % * not / : YALMIP has no ndsdpvar mrdivide

if isnumeric(P) && ~any(P(:))
    % Combinatorially consistent, yet empty: A had no component in the class.
    % A collapsed ndsdpvar lands here too, YALMIP returning a plain double 0.
    warning('PhasorArray:symmetry:emptyProjection', ...
        'Projection onto [%s] is zero; the input has no component in that class.', ...
        strjoin(req, ', '));
end
A = PhasorArray(P);
end

% -------------------------------------------------------------------------
function Y = act(X, code)
% bit 1 = transpose, bit 2 = conjugation, bit 3 = time reversal (k -> -k)
Y = X;
if bitget(code, 1), Y = permute(Y, [2 1 3]); end
if bitget(code, 2), Y = cj(flip(Y, 3));      end
if bitget(code, 3), Y = flip(Y, 3);          end
end

function Y = cj(Y)
if isobject(Y) && ~isa(Y, 'sym')
    Y = real(Y) - 1i*imag(Y);                 % YALMIP has no conj for ndsdpvar
else
    Y = conj(Y);
end
end
