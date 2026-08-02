function A = setHarmonic(A, idx, val, nvp)
%SETHARMONIC  Write harmonic coefficients, keeping A(t) real if it was.
%
%   A = SETHARMONIC(A, IDX, VAL) writes VAL at the harmonics IDX selects. When
%   A(t) is real the opposite harmonics receive conj(VAL) in the same call, so
%   the conjugate symmetry that makes A(t) real is not broken by an edit.
%
%   IDX accepts, from terse to explicit:
%
%     k              every entry of harmonic k          (i and j implied ':')
%     [i j k]        one entry
%     {i, j, k}      any mix, each of i, j, k being ':' , a scalar or a vector
%     {':', j, k}    column j of harmonic k
%
%   A three-element numeric index is always read as [i j k]. To write several
%   harmonics at once, use the cell form: {':', ':', [1 3 5]}. Harmonics beyond
%   the current order are padded in, as brace assignment does.
%
%   Name-Value Arguments:
%     isReal : mirror onto -k with the conjugate. Default: isreal(A), so an
%              array that is real stays real and one that is not is left
%              alone. Pass false to write k only -- deliberately breaking the
%              symmetry is legitimate and stays possible.
%
%   Harmonic 0 is its own mirror, so a real A(t) needs a real DC. A complex
%   value there is written as given and warned about rather than silently
%   coerced.
%
%   Example
%       A = PhasorArray.random(2,2,4);
%       A = setHarmonic(A, 2, 0.3+0.1i);          % harmonic +-2, stays real
%       A = setHarmonic(A, {1,':',3}, [1 2]);     % row 1 of harmonic +-3
%       A = setHarmonic(A, 3, 1i, isReal=false);  % +3 only, A(t) now complex
%
%   See also braceAssign, pad, isreal, mreal.

arguments
    A
    idx
    val
    nvp.isReal (1,1) logical = isreal(A)
end

[n1, n2, k] = normaliseIndex(idx);

kmax = max(abs(k));
if kmax > A.h
    A = A.pad(kmax - A.h);
end

v  = pvalue(A);
h  = A.h;
kp = k + h + 1;                       % array index of harmonic k

v(n1, n2, kp) = val;

if nvp.isReal
    % Mirroring a k that already holds both signs would overwrite what was
    % just written, and which half won would depend on the order.
    if any(ismember(k(k~=0), -k))
        error('PhasorArray:setHarmonic:mirroredTwice', ...
            ['IDX names both +k and -k, so mirroring would overwrite the values ' ...
             'just written. Give one sign and let isReal fill the other, or pass ' ...
             'isReal=false and write both halves yourself.']);
    end
    km = -k + h + 1;                  % array index of harmonic -k
    v(n1, n2, km) = conj(val);
    if any(k == 0) && ~isreal(val)
        warning('PhasorArray:setHarmonic:complexDC', ...
            ['Harmonic 0 is its own mirror, so a real A(t) needs a real DC. ' ...
             'The value was written unchanged; take real(val) or pass isReal=false.']);
    end
end

A = PhasorArray(v);
end

%% =========================================================================
function [n1, n2, k] = normaliseIndex(idx)
%NORMALISEINDEX  Bring the three accepted index forms to {rows, cols, harmonics}.
n1 = ':'; n2 = ':';
if iscell(idx)
    switch numel(idx)
        case 1, k = idx{1};
        case 2, error('PhasorArray:setHarmonic:badIndex', ...
                    'A two-element index is ambiguous. Use k, [i j k] or {i, j, k}.');
        case 3, n1 = idx{1}; n2 = idx{2}; k = idx{3};
        otherwise, error('PhasorArray:setHarmonic:badIndex', ...
                    'IDX must hold at most three subscripts.');
    end
elseif isnumeric(idx) && isscalar(idx)
    k = idx;
elseif isnumeric(idx) && numel(idx) == 3
    n1 = idx(1); n2 = idx(2); k = idx(3);
else
    error('PhasorArray:setHarmonic:badIndex', ...
        ['IDX must be a scalar k, a triple [i j k], or a cell {i, j, k}. ' ...
         'Use the cell form for '':'' or vector subscripts.']);
end
if ~isnumeric(k) || ~isreal(k) || any(mod(k,1) ~= 0)
    error('PhasorArray:setHarmonic:badHarmonic', ...
        'The harmonic index must be an integer or a vector of integers.');
end
end
