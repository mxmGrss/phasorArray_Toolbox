function [P, held] = ndsdpvar(n1,n2,h,nvp)
    % NDSDPVAR Construct an `sdpvar`-based PhasorArray of specified size and symmetry.
    %
    %   P = NDSDPVAR(N1, N2, H, <name-value arguments>) generates an SDP variable-based
    %   PhasorArray suitable for optimization in YALMIP.
    %
    %   Inputs:
    %     N1 - (integer) First dimension size.
    %     N2 - (integer, optional) Second dimension size (default: N1). May differ
    %          from N1, down to N2 = 1 for a column of decision variables.
    %     H  - (integer, optional) Number of harmonics (default: 0).
    %
    %   A rectangular P cannot be symmetric, so the default symmetry drops that
    %   request with a PhasorArray:ndsdpvar:symmetryFallback warning and HELD
    %   reports what survived.
    %
    %   Name-Value Arguments:
    %     'symmetry' (string array) - Symmetry class of P(t) (default:
    %                           ["real" "symmetric"], the Lyapunov variable).
    %                           Pass [] for a full complex P. The 14 names are:
    %                             symmetric      skewSymmetric
    %                             real           imaginary
    %                             even           odd
    %                             hermitian      skewHermitian
    %                             paraSymmetric  skewParaSymmetric
    %                             paraConjugate  skewParaConjugate
    %                             paraHermitian  skewParaHermitian
    %                           PHASORSYMMETRY states what each one means.
    %
    %   Outputs:
    %     P    - (PhasorArray) PhasorArray containing `sdpvar` elements.
    %     HELD - (string array) properties P satisfies, usually more than requested.
    %
    %   Symmetries are stated on P(t) in the time domain, never on the coefficients.
    %   The coefficient-level condition P_k = P_k' is a property in its own right,
    %   named "paraHermitian", and it is NOT time-domain hermitian:
    %
    %       symmetry = "hermitian"       ->  P(t)' = P(t)
    %       symmetry = "paraHermitian"   ->  P(-t)' = P(t), i.e. every P_k hermitian
    %
    %   Whenever the request maps onto a native YALMIP declaration the variable is
    %   declared with it, so no degree of freedom is created only to be projected out.
    %
    %   There is no "diagonal" or "triangular" symmetry: those are sparsity
    %   patterns, not relations between P(t) and a transformed copy of itself, so
    %   they do not compose in the closure the other names live in. Build a
    %   diagonal variable from a column instead, which already costs exactly the
    %   column: diag(PhasorArray.ndsdpvar(3,1,h)) is 3-by-3 with 3*(2h+1) free
    %   entries, against 9*(2h+1) for the full square.
    %
    %   Example:
    %     P = PhasorArray.ndsdpvar(4,4,5);                             % real symmetric
    %     P = PhasorArray.ndsdpvar(4,4,5, "symmetry", "hermitian");     % P(t)' = P(t)
    %     A = PhasorArray.ndsdpvar(4,4,5, "symmetry", "real");          % real, no more
    %     A = PhasorArray.ndsdpvar(4,4,5, "symmetry", []);              % full complex
    %
    %   See also: sdpvar, phasorSymmetry, symmetryClosure, PhasorArray.random,
    %             PosPart2PhasorArray.
    arguments
        n1
        n2=n1
        h=0
        nvp.symmetry string = ["real" "symmetric"]
    end
    if nargin ==1
        switch numel(n1)
            case 1
            case 2
                n2 = n1(2);
                n1 = n1(1);
                h = 0;
            case 3
                n2 = n1(2);
                h = n1(3);
                n1 = n1(1);
        end
    end

    req = nvp.symmetry;
    if n1 ~= n2
        transposing = intersect(req, ["symmetric" "skewSymmetric" "hermitian" ...
            "skewHermitian" "paraSymmetric" "skewParaSymmetric" ...
            "paraHermitian" "skewParaHermitian"]);
        if ~isempty(transposing)
            warning('PhasorArray:ndsdpvar:symmetryFallback', ...
                'Non-square %dx%d: dropping [%s], which needs a square array.', ...
                n1, n2, strjoin(transposing, ', '));
            req = setdiff(req, transposing, 'stable');
        end
    end

    [pType, useReal, native] = nativeDeclaration(req);

    if useReal
        P1 = ndsdpvar(n1,n2,1,pType,'real');
        if h > 0
            P2 = ndsdpvar(n1,n2,h,pType,'complex');
            P  = PosPart2PhasorArray(P1,P2);
        else
            P  = PhasorArray(P1);
        end
    else
        P = PhasorArray(ndsdpvar(n1,n2,2*h+1,pType,'complex'));
    end

    if ~native
        P = phasorSymmetry(P, req);
    end
    if nargout > 1
        held = symmetryClosure(req);
    end
end

%% =========================================================================
function [pType, useReal, native] = nativeDeclaration(req)
%NATIVEDECLARATION  Largest part of REQ that YALMIP can declare directly.
%   YALMIP structures act page by page, so they express the per-coefficient
%   involutions only: 'symmetric' is P_k = P_k.', 'hermitian' is P_k = P_k',
%   which is paraHermitian on P(t). Conjugate symmetry across +-k is the
%   separate 'real' construction. NATIVE is true when nothing else is left.

pType   = 'full';
useReal = any(req == "real");
covered = req(req == "real");

pages = ["symmetric" "symmetric" ; "skewSymmetric" "skew" ; "paraHermitian" "hermitian"];
for k = 1:size(pages,1)
    if any(req == pages(k,1))
        pType   = char(pages(k,2));
        covered = [covered pages(k,1)];  %#ok<AGROW>
        break
    end
end
native = isempty(setdiff(req, covered));
end
