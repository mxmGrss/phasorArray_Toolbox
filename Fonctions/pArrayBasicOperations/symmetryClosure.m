function [held, eps_, names] = symmetryClosure(req)
%SYMMETRYCLOSURE  Properties implied by a set of requested symmetries.
%   HELD = SYMMETRYCLOSURE(REQ) returns every property a PhasorArray built to
%   satisfy REQ also satisfies. HELD contains REQ and is usually larger:
%   "real" and "symmetric" together force "hermitian".
%
%   [HELD, EPS] = SYMMETRYCLOSURE(REQ) also returns the sign of each of the 8
%   involutions of (Z/2)^3, indexed by 1+code with code = bit1 transpose,
%   bit2 conjugation, bit3 time reversal; NaN where the involution is
%   unconstrained.
%
%   [HELD, EPS, NAMES] = SYMMETRYCLOSURE() returns the 14 valid names.
%
%   Because the satisfied set is a subgroup with a character, HELD always has
%   0, 1, 3 or 7 entries, never 2, 4, 5 or 6.
%
%   Contradictory requests error out: SYMMETRIC and SKEWSYMMETRIC together
%   admit only the zero array.
%
%   See also phasorSymmetry, PhasorArray.random, PhasorArray.ndsdpvar.

arguments
    req string = strings(1,0)
end
req(req == "" | ismissing(req)) = [];

names = ["symmetric"     "skewSymmetric"     "real"          "imaginary" ...
         "even"          "odd"               "hermitian"     "skewHermitian" ...
         "paraSymmetric" "skewParaSymmetric" "paraConjugate" "skewParaConjugate" ...
         "paraHermitian" "skewParaHermitian"];
code  = [1 1  2 2  4 4  3 3  5 5  6 6  7 7];
sign_ = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1];

[tf, loc] = ismember(req, names);
if ~all(tf)
    error('PhasorArray:symmetry:unknownName', ...
        'Unknown symmetry "%s". Valid names: %s.', req(find(~tf,1)), strjoin(names, ', '));
end

eps_ = nan(1,8);
eps_(1) = 1;                                  % identity
for g = loc
    for j = find(~isnan(eps_))
        k = 1 + bitxor(j-1, code(g));
        v = eps_(j) * sign_(g);
        if ~isnan(eps_(k)) && eps_(k) ~= v
            error('PhasorArray:symmetry:contradiction', ...
                'Requested symmetries [%s] admit only the zero array.', strjoin(req, ', '));
        end
        eps_(k) = v;
    end
end

held = names(arrayfun(@(k) eps_(code(k)+1) == sign_(k), 1:numel(names)));
end
