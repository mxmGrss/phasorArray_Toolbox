function v = solvedPayload(v, caller, alternative)
%SOLVEDPAYLOAD  Numeric coefficients of a phasor array, or a clear refusal.
%
%   V = SOLVEDPAYLOAD(V, CALLER) returns V unchanged when it is numeric. When
%   it holds YALMIP decision variables it returns their current value, which
%   is what an operation like inv or det needs, but only if that value exists:
%   before the problem is solved YALMIP answers NaN, and the caller would
%   otherwise return a numeric array of NaN with no indication that the
%   question could not be answered. A sym payload is refused outright.
%
%   V = SOLVEDPAYLOAD(V, CALLER, ALTERNATIVE) names a function that does work
%   symbolically, so the message can point at it.
%
%   These operations have no meaning on an unsolved variable. inv(P) is not an
%   expression YALMIP can carry, det(P) would be a polynomial in every entry,
%   and expm(P) is not expressible at all. Solve first, then call sdpval and
%   operate on the result.
%
%   See also sdpval, value, PhasorInv, PhasorDet, PhasorArray.expm.

arguments
    v
    caller (1,1) string
    alternative (1,1) string = ""
end

if isnumeric(v)
    return
end

hint = "";
if alternative ~= ""
    hint = " Use " + alternative + " to work symbolically.";
end

% The identifier has to be a literal for error(), so it is built first;
% passing a format string there leaves the identifier empty.
if isa(v, 'sym')
    error(sprintf('PhasorArray:%s:symbolicPayload', caller), ...
        '%s needs numeric coefficients; this array holds sym.%s', caller, hint);
end

if isa(v, 'ndsdpvar') || isa(v, 'sdpvar')
    v = value(v);
    if any(isnan(v(:)))
        error(sprintf('PhasorArray:%s:unsolvedVariable', caller), ...
            ['%s needs numeric coefficients and this array holds unsolved ' ...
             'decision variables. Solve the problem first, then take sdpval ' ...
             'of the variable and call %s on that.%s'], caller, caller, hint);
    end
end
end
