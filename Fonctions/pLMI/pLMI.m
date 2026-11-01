classdef pLMI
    %PLMI  Collection of periodic matrix inequalities M(theta) >= 0.
    %
    %   C = pLMI(M) represents M(theta) >= 0 for all theta, with M a square
    %   PhasorArray, Hermitian in time and affine in YALMIP decision variables.
    %   Comparisons A >= B and A <= B between PhasorArray operands holding
    %   decision variables return a pLMI; strict inequalities are refused.
    %   [C1, C2] or [C1; C2] concatenates collections.
    %
    %   Compilation to YALMIP constraints:
    %     [F, compiled] = C.gram()        exact Gram certificates (GRAMCONSTRAINT)
    %     [F, compiled] = C.gram("real")  real Gram basis, real symmetric symbols
    %     [F, compiled] = C.section(h)    finite sections T_tb(M,h) >= 0, a
    %                                     necessary condition only
    %   After optimize, compiled.check(sol) returns one report per inequality
    %   (see PLMICOMPILATION).
    %
    %   Factories:
    %     [C, P] = pLMI.lyapunov(A, Q, omega, degree, margin)
    %     [C, P] = pLMI.riccati(A, B, Q, R, omega, degree, margin)
    %
    %   A scalar right-hand side follows the YALMIP convention and applies to
    %   every entry: P >= 1 means P - ones(n) >= 0. Use P >= mu*eye(n) for a
    %   spectral margin. P.d() differentiates in phase; omega*P.d() in time.
    %
    %   EXAMPLE
    %     A = PhasorArray([-1 1;-1 -2]) + 0.2*PhasorArray.cos()*[0 1;1 0];
    %     P = PhasorArray.ndsdpvar(2, 2, 1);
    %     C = [P >= 0.1*eye(2), P.d() + A'*P + P*A <= -0.1*eye(2)];
    %     [F, compiled] = C.gram("real");
    %     sol = optimize(F, [], sdpsettings('solver', 'mosek'));
    %     report = compiled.check(sol);
    %
    %   Requires YALMIP and an SDP solver.
    %
    %   See also: pLMICompilation, gramConstraint, checkGramCertificate, ndsdpvar.
    properties (SetAccess = private)
        Symbols cell = {}
    end
    methods
        function obj = pLMI(M)
            arguments
                M = []
            end
            if isempty(M), return; end
            if ~isa(M,'PhasorArray'), M = PhasorArray(M); end
            plmiValidateSymbol(M);
            obj.Symbols = {M};
        end
        function n = count(obj)
            %COUNT Number of inequalities in the collection.
            arguments
                obj (1,1) pLMI
            end
            n = numel(obj.Symbols);
        end
        function combined = horzcat(varargin)
            combined = pLMI();
            for j = 1:nargin
                item = varargin{j};
                if isnumeric(item) && isempty(item), continue; end
                if ~isa(item,'pLMI') || ~isscalar(item)
                    error('pLMI:collection','Combine only scalar pLMI collections.');
                end
                combined.Symbols = [combined.Symbols,item.Symbols];
            end
        end
        function combined = vertcat(varargin)
            combined = horzcat(varargin{:});
        end
        function [F, compiled] = gram(obj,basis)
            %GRAM Impose every inequality exactly by a Gram certificate.
            %   basis "hermitian" (default) or "real" (real symmetric symbols).
            arguments
                obj (1,1) pLMI
                basis (1,1) string {mustBeMember(basis,["real","hermitian"])} = "hermitian"
            end
            F = []; grams = cell(1,obj.count()); metadata = grams;
            for j = 1:obj.count()
                M = obj.Symbols{j};
                if basis == "real", plmiValidateSymbol(M,true); end
                [Fj,grams{j},metadata{j}] = M.gramConstraint(0,basis);
                F = [F,Fj]; %#ok<AGROW>
            end
            compiled = pLMICompilation("gram",obj.Symbols,grams,metadata,[]);
        end
        function [F, compiled] = section(obj,h)
            %SECTION Impose the finite sections T_tb(M,h) >= 0 (necessary only).
            arguments
                obj (1,1) pLMI
                h (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite}
            end
            F = []; sections = cell(1,obj.count());
            for j = 1:obj.count()
                H = obj.Symbols{j}.T_tb(h);
                H = (H+H')/2;
                sections{j} = H;
                if isnumeric(H)
                    if min(eig(H)) < 0
                        error('pLMI:constantSectionInfeasible','The constant section is not positive semidefinite.');
                    end
                else
                    F = [F,H >= 0]; %#ok<AGROW>
                end
            end
            compiled = pLMICompilation("section",obj.Symbols,sections,{},h);
        end
        function disp(obj)
            fprintf('pLMI: %d periodic matrix inequalities, not yet compiled.\n',obj.count());
            for j = 1:obj.count()
                M = obj.Symbols{j};
                fprintf('  %d: %d-by-%d Hermitian symbol, degree %d, M(theta) >= 0\n', ...
                    j,size(M,1),size(M,2),M.h);
            end
        end
    end
    methods (Static)
        function [obj,P] = lyapunov(A,Q,omega,degree,margin)
            %LYAPUNOV P >= margin*I and omega*P.d() + A'*P + P*A <= -Q.
            %   A and Q are real periodic data, Q symmetric; P is a new
            %   ndsdpvar of the given degree.
            arguments
                A
                Q
                omega (1,1) double {mustBePositive,mustBeFinite}
                degree (1,1) double {mustBeInteger,mustBeNonnegative}
                margin (1,1) double {mustBeNonnegative,mustBeFinite} = 0
            end
            [A,Q] = validatePeriodicData(A,Q);
            n = size(A,1); P = PhasorArray.ndsdpvar(n,n,degree);
            obj = [P >= margin*eye(n),omega*P.d()+A'*P+P*A <= -Q];
        end
        function [obj,P] = riccati(A,B,Q,R,omega,degree,margin)
            %RICCATI Riccati inequality in Schur form, with P >= margin*I:
            %   [omega*P.d() + A'*P + P*A + Q, P*B; B'*P, R] >= 0.
            %   Maximizing the mean trace of P approaches the stabilizing
            %   solution from below. R is constant, symmetric positive definite.
            arguments
                A
                B
                Q
                R (:,:) double
                omega (1,1) double {mustBePositive,mustBeFinite}
                degree (1,1) double {mustBeInteger,mustBeNonnegative}
                margin (1,1) double {mustBeNonnegative,mustBeFinite} = 0
            end
            [A,Q] = validatePeriodicData(A,Q);
            B = validatePeriodicReal(B);
            n = size(A,1); m = size(B,2);
            if size(B,1) ~= n || ~isequal(size(R),[m m]) || ...
                    any(~isfinite(R(:))) || ~isreal(R) || norm(R-R','fro') > 1e-12
                error('pLMI:riccatiData','B and R must have compatible dimensions and R must be finite real symmetric.');
            end
            [~,flag] = chol(R);
            if flag ~= 0, error('pLMI:riccatiR','R must be positive definite.'); end
            P = PhasorArray.ndsdpvar(n,n,degree);
            top = omega*P.d()+A'*P+P*A+Q;
            block = [top,P*B;B'*P,PhasorArray(R)];
            obj = [P >= margin*eye(n),block >= 0];
        end
        function tf = hasDecisionPayload(varargin)
            %HASDECISIONPAYLOAD True if any operand holds YALMIP variables.
            tf = false;
            for j = 1:nargin
                payload = varargin{j};
                if isa(payload,'PhasorArray'), payload = value(payload); end
                tf = tf || isa(payload,'sdpvar') || isa(payload,'ndsdpvar');
            end
        end
        function obj = fromComparison(left,right,sense)
            %FROMCOMPARISON pLMI of left >= right ("ge") or left <= right ("le").
            arguments
                left
                right
                sense (1,1) string {mustBeMember(sense,["ge","le"])}
            end
            left = PhasorArray(left); right = PhasorArray(right);
            leftSize = [size(left,1),size(left,2)];
            rightSize = [size(right,1),size(right,2)];
            if ~isequal(leftSize,rightSize) && prod(leftSize) ~= 1 && prod(rightSize) ~= 1
                error('pLMI:dimensions','Operands must have equal matrix dimensions or be scalar.');
            end
            if sense == "ge", M = left-right; else, M = right-left; end
            obj = pLMI(M);
        end
    end
end
