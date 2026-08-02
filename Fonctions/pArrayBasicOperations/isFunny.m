function tf = isFunny(x)
%ISFUNNY True for "funny" element types that cannot use the vectorized double
%   fast path (sym, sdpvar, ndsdpvar) and must take the generic cat-based route.
%
%   tf = ISFUNNY(x) returns true when x is a symbolic (sym) or YALMIP
%   decision (sdpvar / ndsdpvar) array. These types do not support the
%   page-wise builtins (pagemtimes, tensorprod) nor plain typed allocation,
%   so the toolbox dispatches them onto concatenation-based code paths.
%
%   This factors the triple isa(...) test that is repeated across the
%   pArrayBasicOperations dispatch sites. It intentionally does NOT cover the
%   sdpvar/ndsdpvar-only checks (use isa directly when sym must be excluded).
%
%   See also: isa, PhasorArrayTimes, sparray2TBlocks.
arguments
    x
end
tf = isa(x, 'sym') || isa(x, 'sdpvar') || isa(x, 'ndsdpvar');
end
