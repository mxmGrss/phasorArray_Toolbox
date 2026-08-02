function warnIfNotConverged(info, methodCall)
%WARNIFNOTCONVERGED  Point at the method form when an operator solve falls short.
%
%   WARNIFNOTCONVERGED(INFO, METHODCALL) warns when INFO.status says the
%   harmonic-order refinement did not converge, and names the call that accepts
%   the knobs. The operators \ and / take no name-value pairs, so a caller who
%   hits a stagnated or truncated solve has no way to act on it from the
%   operator syntax alone.
%
%   Statuses: 0 converged, 1 stagnated, 2 maxh reached, 3 fixed h,
%   4 convergence judged unreachable.
%
%   See also mlHmcDivide, mrHmcDivide, adaptiveHSolve.

arguments
    info
    methodCall (1,1) string
end

if ~isstruct(info) || ~isfield(info, 'status') || info.status == 0 || info.status == 3
    return
end

why = "did not converge";
switch info.status
    case 1, why = "stagnated";
    case 2, why = "hit the harmonic-order ceiling";
    case 4, why = "judged the target order unreachable";
end

res = NaN;
if isfield(info, 'resrelnorm'), res = info.resrelnorm; end
hEnd = NaN;
if isfield(info, 'h'), hEnd = info.h; end

warning('PhasorArray:divide:notConverged', ...
    ['The harmonic solve %s at h = %d (relative residual %.2e).\n' ...
     '  The operator form takes no options. Call the method to tune it:\n' ...
     '    %s\n' ...
     '  Options that usually help, in order:\n' ...
     '    ''maxh'', <larger>              raise the ceiling\n' ...
     '    ''thresholdResidual'', <looser> accept a coarser answer\n' ...
     '    ''stagnationRatio'', <smaller>  keep going through a flat stretch\n' ...
     '    ''verbose'', true               print the refinement table\n' ...
     '  Suppress: warning(''off'',''PhasorArray:divide:notConverged'')'], ...
    why, hEnd, res, methodCall);
end
