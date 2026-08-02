function [HpJ,JHm,Hp,Hm] = spArray2BTHankel(Aph,m,nvp)
% SPARRAY2BTHANKEL Sparse Block Toeplitz Hankel matrices from a 3D array.
%
%   Sparse counterpart of ARRAY2BTHANKEL: same outputs as SPARRAY2TBHANKEL
%   but in the Block-Toeplitz (BT, harmonic-outer) layout, obtained by
%   perfect-shuffle index permutation (sparsity preserved).
%
%   See also: Array2BTHankel, spArray2TBHankel, spBTHankel.
arguments
    Aph
    m (1,1) {mustBeInteger, mustBeNonnegative}
    nvp.method {mustBeMember(nvp.method,{'cell2mat','cat'})} = 'cell2mat'
end
[HpJ,JHm,Hp,Hm] = spArray2TBHankel(Aph, m, "method", nvp.method);

nb = m + 1;
vr = shufflePerm(size(HpJ,1)/nb, nb);
vc = shufflePerm(size(HpJ,2)/nb, nb);
HpJ = HpJ(vr, vc);
JHm = JHm(vr, vc);
Hp  = Hp(vr, vc);
Hm  = Hm(vr, vc);
end

function v = shufflePerm(p, q)
v = reshape(reshape(1:p*q, q, p)', [], 1);
end
