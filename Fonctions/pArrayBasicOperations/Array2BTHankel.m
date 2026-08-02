function [HpJ,JHm,Hp,Hm] = Array2BTHankel(Aph,m,nvp)
% ARRAY2BTHANKEL Construct Block Toeplitz Hankel matrices from a 3D array.
%
%   Same outputs as ARRAY2TBHANKEL but in the Block-Toeplitz (BT) layout:
%   harmonic-outer blocks of size N×N, instead of the entrywise Toeplitz
%   blocks of the TB layout. Obtained by perfect-shuffle conjugation of the
%   TB-layout matrices (pure permutation, exact).
%
%   Inputs:
%     Aph     - (3D array or PhasorArray) Fourier coefficient array.
%     m       - (integer) Truncation order for harmonics.
%     nvp    - 'method': 'cell2mat' (default) or 'cat' (symbolic/SDP).
%
%   Outputs:
%     HpJ, JHm, Hp, Hm - ((m+1)N × (m+1)N) matrices, BT layout.
%
%   See also: Array2TBHankel, spArray2BTHankel, BTHankel, TB_2_BT.
arguments
    Aph
    m (1,1) {mustBeInteger, mustBeNonnegative}
    nvp.method {mustBeMember(nvp.method,{'cell2mat','cat'})} = 'cell2mat'
end
[HpJ,JHm,Hp,Hm] = Array2TBHankel(Aph, m, "method", nvp.method);

% state-outer -> harmonic-outer: same permutation as TB_2_BT, applied by
% indexing so that sparse and symbolic inputs are preserved.
nb = m + 1;                                     % Hankel block size per entry
vr = shufflePerm(size(HpJ,1)/nb, nb);
vc = shufflePerm(size(HpJ,2)/nb, nb);
HpJ = HpJ(vr, vc);
JHm = JHm(vr, vc);
Hp  = Hp(vr, vc);
Hm  = Hm(vr, vc);
end

function v = shufflePerm(p, q)
% Row-selection order of shuffle_matrix(p,q): rows qi:q:p*q for qi = 1..q.
v = reshape(reshape(1:p*q, q, p)', [], 1);
end
