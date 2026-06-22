function P = ndsdpvar(n1,n2,h,varg)
    % NDSPDPVAR Construct an `sdpvar`-based PhasorArray of specified size and structure.
    %
    %   P = NDSPDPVAR(N1, N2, H, <name-value arguments>) generates an SDP variable-based
    %   PhasorArray suitable for optimization in YALMIP.
    %
    %   Inputs:
    %     N1 - (integer) First dimension size.
    %     N2 - (integer, optional) Second dimension size (default: N1).
    %     H  - (integer, optional) Number of harmonics (default: 0).
    %
    %   Name-Value Arguments:
    %     'PhasorType' (char) - Defines the structure of the phasor (default: 'symmetric').
    %                           Options: 'symmetric', 'full', 'diagonal', etc.
    %                           See YALMIP documentation for the complete list.
    %     'real' (logical)    - If true, ensures conjugate symmetry for real-valued signals (default: true).
    %
    %   Outputs:
    %     P - (PhasorArray) PhasorArray containing `sdpvar` elements.
    %
    %   Notes:
    %     - The `PhasorType` argument defines **phasor structure**, not time-domain structure.
    %     - Example: A **Hermitian** phasor structure enforces:
    %          conj(A_(ij)(t)) = A_(ji)(-t)
    %       However, for A(t) to be Hermitian in the **time domain**, phasors must satisfy:
    %          A_k = ctrans(A_-k)
    %     - If `real=true`, ensures the phasor array represents a real-valued periodic matrix.
    %     - If `h>0`, higher-order harmonics are created as additional `sdpvar` variables.
    %
    %   Example:
    %     P = PhasorArray.ndsdpvar(4,4,5, 'PhasorType', 'symmetric', 'real', true);
    %       -> Produces a real-valued P(t), with 5 harmonics (size 11 along the third dimension),
    %          and enforces symmetry (i.e., P_ij(t) = P_ji(t)).
    %
    %     A = PhasorArray.ndsdpvar(4,4,5, 'PhasorType', 'full', 'real', true);
    %       -> Produces a real-valued A(t) with no additional structure constraints.
    %
    %   See also: sdpvar, PhasorArray, PosPart2PhasorArray.
    arguments
        n1
        n2=n1
        h=0
        varg.PhasorType='symmetric'
        varg.real=true
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
    if n1~=n2
        if ismember(varg.PhasorType,{'symmetric'})
            warning('PhasorArray:randomPhasorArray:phasorTypeFallback', 'Non-square matrix: PhasorType changed from ''symmetric'' to ''full''.')
            varg.PhasorType='full';
        end
    end
    if varg.real
        P1=(ndsdpvar(n1,n2,1,varg.PhasorType,'real'));
        if h>0
            P2=(ndsdpvar(n1,n2,h,varg.PhasorType,'complex'));
            P = PosPart2PhasorArray(P1,P2);
        else
            P=PhasorArray(P1);
        end
    else
        P=PhasorArray(ndsdpvar(n1,n2,2*h+1,varg.PhasorType,'complex'));
    end
end
