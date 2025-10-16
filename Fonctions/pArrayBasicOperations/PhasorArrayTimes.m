function D = PhasorArrayTimes(A,B,m,varg)
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

arguments
    A
    B
    m=[]
    varg.reduce = 0
    varg.output {mustBeMember(varg.output,["Array","PhasorArray","Herited"])} = "Herited"
end
if ~isMATLABReleaseOlderThan("R2022a")
    if isa(A,'PhasorArray')
        A=A.Value;
        if strcmp(varg.output,"Herited")
            varg.output="PhasorArray";
        end
    else
        if strcmp(varg.output,"Herited")
            varg.output="Array";
        end
    end

    if isa(B,'PhasorArray')
        B=B.Value;
    end

    %switch to function adapted to sdpvar
    if isa(A,"sdpvar") || isa(A,"ndsdpvar") || isa(B,"sdpvar") || isa(B,"ndsdpvar") || isa(A,"sym")  || isa(B,"sym")
        D=PhasorArrayTimes2(A,B,m);
        return
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
        D=A*B;
        return
        A=A*eye(size(B,1));
    end

    if isscalar(B)
        D=A*B;
        return
        B=B*eye(size(A,2));
    end

    if nA1==1 && nA2==1 && mA>0
        A=pagemtimes(A,eye(size(B,1)));
        %update the size of A
        nA=size(A);
        nA1=nA(1);
        nA2=nA(2);
    end

    if nB1==1 && nB2==1 && mB>0
        B=pagemtimes(B,eye(size(A,2)));
        %update the size of B
        nB=size(B);
        nB1=nB(1);
        nB2=nB(2);
    end

    %finally assert that the size of the two arrays are compatible
    if nA2~=nB1
        error('The two arrays must have compatible sizes. Left array %s is %d x %d x %d and right array is  %s %d x %d x %d',inputname(1),nA1,nA2,size(A,3),inputname(2),nB1,nB2,size(B,3));
    end

    if isempty(m)
        m=mA+mB ;
    end
    D=zeros(size(A,1),size(B,2),2*m+1);

    %Convolution using tensorprod
    for k=(-m:m)
        l1=max(k-mB,-mA);
        l2=min(k+mB,mA);
        V=A(:,:,(mA+1+l1):(mA+1+l2));
        U=B(:,:,(k+mB+1-l1):-1:(k+mB+1-l2));
        D(:,:,m+1+k)=tensorprod(V,U,[2,3],[1,3]);
    end
    if varg.reduce
        D=ReduceArray(D);
    end

else
    warning('Error : Switching to legacy matricial convolution')
    D=PhasorArrayTimes2(A,B,m);
end

end
