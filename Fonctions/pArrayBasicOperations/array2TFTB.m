function r= array2TFTB(o1,m)
%F_tb Compute the Fourier representation of A in a form compatible with T_tb(A, m).
%
%   This function computes the **Fourier representation** of A(t) up to order `m`, 
%   but instead of forming a **square Toeplitz Block matrix** (T_tb(A, m)), it **stacks
%   the phasors of each element of A into a structured vectorized form**.
%
%   **Definition:**
%   - Let `a(t)` be a scalar function with Fourier coefficients `(a_k)_{k∈ℤ}`.
%   - Then, `F_tb(a, m) = [a_(-m); a_(-m+1); ... a_0; a_1; ... a_m]`.
%   - If `A(t)` is an `N×M` matrix, its Fourier representation is given by:
%       F_tb(A, m) = [F_tb(A_11), F_tb(A_12); F_tb(A_21), F_tb(A_22)].
%
%   **Key Property (TB Compatibility):**
%   If `y(t) = A(t)x(t)`, then:
%       F_TB(y) = T_tb(A, m) ⋅ F_TB(x),
%   where `T_tb(A, m)` is the **Toeplitz Block representation** of `A` (see `TB` method).
%
%   **Dimension of the Output:**
%   - If `A` is an `N×M` matrix, then `F_tb(A, m)` is a `((2m+1)N) × M` matrix.
%
%   Syntax:
%   r = F_tb(A, m)  
%       Computes the Fourier representation of `A` up to order `m`, compatible with T_tb(A, m).
%
%   Input Arguments:
%   - o1 (PhasorArray) : The input PhasorArray representing A(t).
%   - m (integer, optional) : The highest harmonic order to retain. Default: `A.h`.
%
%   Output:
%   - r ((2m+1)N × M matrix) : The Fourier representation of A in a form compatible with T_tb(A, m).
%
%   Example:
%   % Compute Fourier representation of A in TB form
%   A = PhasorArray(rand(4,4,11));
%   F_tb_A = F_tb(A, 5);
%
%   See also: F_bt, TB.
arguments
    o1
    m=(size(o1,3)-1)/2;
end


if isa(o1,'PhasorArray')
    o1=o1.Value;
end

nho=(size(o1,3)-1)/2;


if nho<m
    toto = PhasorArrayPad(o1,m-nho);
elseif nho>m
    toto = o1(:,:,(nho+1)+(-m:m));
end

titi=permute(toto,[3 1 2]);
r=reshape(titi,[],size(o1,2),1);
end