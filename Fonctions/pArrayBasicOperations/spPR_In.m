function [outM] = spPR_In(A,nxB,h)
%produit rond  I_nB o A 
%where A is a block matrix of type 'TB'
% B o A = [ B kron A11  | B kron A12  | ... | B kron A1(ny)   ]
%         [ B kron A21  |                                     ]
%         [                                                   ]
%         [ B kron Anx1 | B kron Anx2 | ... | B kron A(nx)(ny)]
%
%h is the harmonics number of the truncature ie in each bloc of A.  2*h+1 is the length of
%the blocs
%
%    A can be provided as a 3D array of phasors, then A=array2TB(A, 2*h) is
%    used
%

if isa(A,'PhasorArray')
    A=A.Value;
    A=sparray2TBlocks(A,2*h);
end

if ndims(A)==3
    A=sparray2TBlocks(A,2*h);
end

nxA=size(A,1)/(2*h+1);
nyA=size(A,2)/(2*h+1);
% outM = sparse(size(A).*[nxB nxB]);
% outM = spalloc(size(A,1)*nxB,size(A,2)*nxB,numel(A)*nxB);
% nyB=nxB;
% xM=nxB*(2*h+1);
% yM=nyB*(2*h+1);

c=cell(nxA,nyA);
for nxii=1:nxA
    for nyjj=1:nyA
        BlocijA = A(((nxii-1)*(2*h+1))+(1:(2*h+1)),((nyjj-1)*(2*h+1))+(1:(2*h+1)));
        tmp = repmat({BlocijA},nxB,1);
        dB = blkdiag(tmp{:});
        c{nxii,nyjj}=dB;
%         outM((nxii-1)*xM+(1:xM),(nyjj-1)*yM+(1:yM))=dB;
    end
outM=cell2mat(c);
end