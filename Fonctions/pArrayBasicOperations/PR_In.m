function [outM] = PR_In(A,nxB,h)
%PR_IN Produit rond (round product) I_nB o A, block-wise Kronecker on a TB matrix.
%where A is a block matrix of type 'TB'
% B o A = [ B kron A11  | B kron A12  | ... | B kron A1(ny)   ]
%         [ B kron A21  |                                     ]
%         [                                                   ]
%         [ B kron Anx1 | B kron Anx2 | ... | B kron A(nx)(ny)]
%
%h is the number of harmonics in each bloc of A. ie 2*h+1 is the length of
%the blocs
%
%    A can be provided as a 3D array of phasors, then A=array2TB(A, 2*h) is
%    used
%
arguments
    A
    nxB (1,1) double
    h (1,1) double
end

if isscalar(h)
    hOut = h;
    hIn = h;
else
    hOut = h(1);
    hIn = h(2);
end

if isa(A,'PhasorArray')
    A=A.Value;
    A=sparray2TBlocks(A,[hOut,hIn]);
end

if ndims(A)==3
    A=sparray2TBlocks(A,[hOut,hIn]);
end

nxA=size(A,1)/(2*hOut+1);
nyA=size(A,2)/(2*hIn+1);
% outM = sparse(size(A).*[nxB nxB]);
outM = spalloc(size(A,1)*nxB,size(A,2)*nxB,numel(A)*nxB);
nyB=nxB;
xM=nxB*(2*hOut+1);
yM=nyB*(2*hIn+1);

for nxii=1:nxA
    for nyjj=1:nyA
        BlocijA = A(((nxii-1)*(2*hOut+1))+(1:(2*hOut+1)),((nyjj-1)*(2*hIn+1))+(1:(2*hIn+1)));
        tmp = repmat({BlocijA},nxB,1);
        dB = blkdiag(tmp{:});
        outM((nxii-1)*xM+(1:xM),(nyjj-1)*yM+(1:yM))=dB;
    end
end
