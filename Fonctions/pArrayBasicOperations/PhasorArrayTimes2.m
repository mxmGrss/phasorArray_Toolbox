function D = PhasorArrayTimes2(A,B,m)
%Takes 2 3D array as input, seen as a Matrix sequence indexed on Z and
%truncated, stored along the third dimension. 
%Compute the matrix convolution
% It is equivalent of computing the phasor of C(t) = A(t)B(t) from the
% phasor of A(t) and B(t)
%
% C = A conv B 
% C(:,:,mC+1+kk) = sum_ll A(:,:,mA+1+ll) * B(:,:,mB+1+kk-ll) 
%
% A is a 3D array of size m x n x mA
% B is a 3D array of size n x p x mB
%
% This function can be used with ndsdpvar arrays

arguments
    A
    B
    m=[]
end

if isa(A,'PhasorArray')
    A=A.Value;
end

if isa(B,'PhasorArray')
    B=B.Value;
end


nA=size(A);
nB=size(B);
nA1=nA(1);
nA2=nA(2);
nB1=nB(1);
nB2=nB(2);

    mA=(size(A,3)-1)/2;
    mB=(size(B,3)-1)/2;


if isscalar(A)
    A=A*eye(size(B,1));
end

if isscalar(B)
    B=B*eye(size(A,2));
end

if nA1==1 && nA2==1 && mA>0
    if isa(A,'sdpvar') || isa(A,'ndsdpvar')
        %Ad=ndsdpvar(size(B,1),size(B,1),2*mA+1);
        Ad =  zero(size(B,1),size(B,1),2*mA+1);
        for ii=-mA:mA
            Ad(:,:,ii+mA+1)= Ad(:,:,ii+mA+1) + eye(size(B,1))*A(1,1,ii+mA+1);
        end
        A=Ad;
    elseif isa(A,'sym') 
        Ad=sym('Ad',[size(B,1),size(B,1),2*mA+1]);
        for ii=-mA:mA
            Ad(:,:,ii+mA+1)=eye(size(B,1))*A(1,1,ii+mA+1);
        end
        A=Ad;
    else
    A=pagemtimes(A,eye(size(B,1)));
    end
end

if nB1==1 && nB2==1 && mB>0
    if isa(B,'sdpvar') || isa(B,'ndsdpvar')
       Bd=ndsdpvar(size(A,2),size(A,2),1+2*mB)*0;
        % Bd =  zero(size(A,2),size(A,2),1+2*mB);
        for ii=-mB:mB
            Bd(:,:,ii+mB+1)=Bd(:,:,ii+mB+1)+eye(size(A,2))*B(1,1,ii+mB+1);
        end
        B=Bd;
    elseif isa(B,'sym') 
        Bd=sym('Bd',[size(A,2),size(A,2),1+2*mB]);
        for ii=-mB:mB
            Bd(:,:,ii+mB+1)=eye(size(A,2))*B(1,1,ii+mB+1);
        end
        B=Bd;
    else
    
    B=pagemtimes(B,eye(size(A,2)));
    end
end



    if isempty(m)
        m=mA+mB ;
    end

    
% % D=sdpvar(size(A,1),size(B,2),2*m+1);;
Dd=[];
    %Convolution
    for k=(-m:m)
        l1=max(k-mB,-mA);
        l2=min(k+mB,mA);
        V=A(:,:,(mA+1+l1):(mA+1+l2));
        U=B(:,:,(k+mB+1-l1):-1:(k+mB+1-l2));
        VV=reshape(V,size(V,1), size(V,2)*size(V,3));
        UU=reshape(permute(U,[1 3 2]),size(V,2)*size(V,3),[]);
        DD=VV*UU;
%         D(:,:,m+1+k)=DD;
        Dd=[Dd, DD];
    end

    D=reshape(Dd,size(V,1),size(U,2),[]);
%     D=ReduceArray(D);

% %2 Initialize D as an zero array with the same size as the expected output
% D2 = sdpvar(size(A,1), size(B,2), 2*m+1)*0;
% 
% % Loop over the third dimension of D
% for k = -m:m
%     % Calculate the limits of the convolution sum for the current value of k
%     l1 = max(k-mB, -mA);
%     l2 = min(k+mB, mA);
% 
%     % Initialize the current slice of D as a zero matrix
%     % D2(:,:,m+1+k) = zeros(size(A,1), size(B,2));
% 
%     % Loop over the third dimension of A and B within the calculated limits
%     for i = (mA+1+l1):(mA+1+l2)
%         j = k+mB+1-(i-mA-1);
% 
%         % Extract the current slices of A and B
%         V = A(:,:,i);
%         U = B(:,:,j);
% 
% 
%         % Perform the matrix multiplication and add the result to the current slice of D
%         D2(:,:,m+1+k) = D2(:,:,m+1+k) + V*U;
%     end
% end
% D2;
% D=D2;
end
