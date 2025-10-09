function [outM] = pr(A,B,h)
%produit rond  B o A 
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

if ndims(A)==3
    A=array2TBlocks(A,2*h);
end

nx=size(A,1)/(2*h+1);
ny=size(A,2)/(2*h+1);
% outM = zeros(size(A).*size(B));

outM=[]
for nxii=1:nx
    ddB=[];
    for nyjj=1:ny
        BlocijA = A(((nxii-1)*(2*h+1))+(1:(2*h+1)),((nyjj-1)*(2*h+1))+(1:(2*h+1)));
        dB=kron(B,BlocijA);
        ddB=[ddB dB];
    end
    outM=[outM ;ddB];
end