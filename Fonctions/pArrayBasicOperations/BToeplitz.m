function [Tmat,Mat] = BToeplitz(Lmat1,Lmat2)

%BToeplitz Summary of this function goes here
%   Tmat =BToeplitz(LMat1,LMat2) build the toeplitz block matrix of dim ni*nbloc x nj*nbloc, where Lmat1 is a 3D (ni,nj,nbloc) array of nbloc (ni,nj) matrix and def the first
%   colomn of the resulting toeplitz matrix, with its first element being the upper left element of the
%   block matrix. if specified Lmat2 (ni,nj,nbloc) define the first row of matrix
%   if Lmat2 is left blank, Lmat2 is the block transconjugate of Lmat1, so the resulting
%   matrix is block hermitian off the main diagonal (if Lmat(1) is hermitian, the
%   resulting matrix is hermitian), to force régularisation (ie (Lmat(:,:,1) +
%   Lmat(:,:,1)')/2 is used as diagonal element), specify (1) as a
%   second argument.
%   if Lmati is 0, the resulting matrix is down triangular (i=2) or upper
%   triangular (i=1)
%   first element of Lmat2 is ignored and considered equal to the first of
%   Lmat1
%
%   Hence, for a,b two 3D array of size (ni,nj,nbloc) :
%   - BBToeplitz(a,b) is a block toeplitz matrix determined by a and b
%   -BToeplitz(a) is hermitian off diagonal
%   -BToeplitz(a,1) is hermitian and take the hermitian part of a(:,:,1) on
%   the diagonal
%   -BToeplitz(a,a) is symetric but not hermitian (no conjugation) if a is
%   complex valued, 
%   if a is real,  BToeplitz(a)=BToeplitz(a,1)=BToeplitz(a,a)
%   

    [nx,nj,nk] = size(Lmat1);
    if nargin==1
        Lmat2 = conj(Lmat1(:,:,2:end));
    elseif nargin==2
        if Lmat2==0
            Lmat2=zeros(nx,nj,nk-1);
        elseif Lmat1==0
            [nx,nj,nk] = size(Lmat2);
            Lmat1=zeros(nx,nj,nk);
            Lmat1(:,:,1)=Lmat2(:,:,1);
            Lmat2=Lmat2(:,:,2:end);
        elseif Lmat2==1
            Lmat1(:,:,1)=(Lmat1(:,:,1)+Lmat1(:,:,1)')/2;
        else
            if Lmat2(:,:,1) ~= Lmat1(:,:,1)
                disp('Attention : première valeures différentes, premier élément de Lmat2 ignoré')
            end
            Lmat2=Lmat2(:,:,2:end);
        end
    end
    Lmat2 = flip(Lmat2,3);

    Mat=zeros(nx,nj,2*nk-1);
    Mat(:,:,1:(nk-1))=Lmat2;
    Mat(:,:,nk:end)=Lmat1(:,:,:);

    listMat2=mat2cell(Mat,nx,nj,ones(2*nk-1,1));
    Toe=toeplitz((nk):(2*nk-1),(nk):-1:1);
    Tmat = cell2mat(listMat2(Toe));
    
    end