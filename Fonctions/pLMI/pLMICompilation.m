classdef pLMICompilation
    %PLMICOMPILATION  Result of PLMI.GRAM or PLMI.SECTION, with its validation.
    %
    %   REPORTS = COMPILED.CHECK(SOL) validates the current YALMIP values after
    %   SOL = optimize(...), one report per inequality, with fields
    %     method, solverProblem, lowerBound, sectionMinimum, sectionFeasible,
    %     strictlyCertified, status
    %   status is "strictly_feasible", "not_certified", "section_only",
    %   "solver_failed" or "nonfinite". A Gram report bounds min_t of the
    %   smallest eigenvalue of the symbol (CHECKGRAMCERTIFICATE); a section
    %   report is never a certificate.
    %   Call check right after the corresponding solve: a later solve can change
    %   the values of shared decision variables.
    %
    %   Name-value:
    %     sectionTolerance  Tolerance on the section minimum (default 1e-8).
    %
    %   See also: pLMI, checkGramCertificate.
    properties (SetAccess = private)
        Method string
        Symbols cell
        Matrices cell
        Metadata cell
        Order
        Guarantee string
    end
    methods
        function obj = pLMICompilation(method,symbols,matrices,metadata,order)
            arguments
                method (1,1) string {mustBeMember(method,["gram","section"])}
                symbols cell
                matrices cell
                metadata cell
                order
            end
            obj.Method = method; obj.Symbols = symbols; obj.Matrices = matrices;
            obj.Metadata = metadata; obj.Order = order;
            if method == "gram"
                obj.Guarantee = "Exact polynomial representation; numerical validation required";
            else
                obj.Guarantee = "Necessary finite-section relaxation only";
            end
        end
        function reports = check(obj,sol,nvp)
            arguments
                obj (1,1) pLMICompilation
                sol (1,1) struct
                nvp.sectionTolerance (1,1) double {mustBeNonnegative,mustBeFinite} = 1e-8
            end
            if ~isfield(sol,'problem') || ~isnumeric(sol.problem) || ...
                    ~isscalar(sol.problem) || ~isfinite(sol.problem)
                error('pLMI:diagnostics','Supply the diagnostics returned by optimize.');
            end
            template = struct('method',obj.Method,'solverProblem',sol.problem, ...
                'lowerBound',NaN,'sectionMinimum',NaN,'sectionFeasible',false, ...
                'strictlyCertified',false,'status',"not_checked");
            reports = repmat(template,1,numel(obj.Symbols));
            for j = 1:numel(reports)
                if sol.problem ~= 0
                    reports(j).status = "solver_failed";
                    continue
                end
                matrix = obj.Matrices{j};
                if ~isnumeric(matrix), matrix = value(matrix); end
                if any(~isfinite(matrix(:)))
                    reports(j).status = "nonfinite";
                    continue
                end
                if obj.Method == "gram"
                    M = obj.Symbols{j};
                    if pLMI.hasDecisionPayload(M), M = sdpval(M); end
                    coefficients = value(M);
                    if any(~isfinite(coefficients(:)))
                        reports(j).status = "nonfinite";
                        continue
                    end
                    % M already contains the right-hand side, margin included.
                    bound = M.checkGramCertificate(0,matrix,obj.Metadata{j});
                    reports(j).lowerBound = bound.lower;
                    reports(j).strictlyCertified = bound.lower > 0;
                    reports(j).status = "not_certified";
                    if reports(j).strictlyCertified
                        reports(j).status = "strictly_feasible";
                    end
                else
                    reports(j).sectionMinimum = min(real(eig((matrix+matrix')/2)));
                    reports(j).sectionFeasible = reports(j).sectionMinimum >= -nvp.sectionTolerance;
                    reports(j).status = "section_only";
                end
            end
        end
        function disp(obj)
            fprintf('pLMI compilation: %s, %d inequalities.\n',obj.Method,numel(obj.Symbols));
            fprintf('  %s\n',obj.Guarantee);
        end
    end
end
