function [Xph,M,M1,M2,colQ,colX] = Sylv_harmonique(Ahm,Bhm,Chm,h,omeg,varargin)
%SYLV_HARMONIQUE  Solve Harmonic Sylvester Equations.
%
%   X = Sylv_harmonique(Ahm, Bhm, Chm, h, omeg) solves the Harmonic Sylvester equation:
%   
%       (Ahm + Nh) * X + X * (Bhm - Nh) + Chm = 0
%
%   Where:
%       - Ahm, Bhm, and Chm are 3D arrays representing the harmonic components of the 
%         time-periodic matrices A(t), B(t), and C(t) respectively.
%       - Nh is a diagonal matrix defined by Nh = diag( jk * omeg ) where |k| <= h.
%       - omeg represents the pulsation in Nh.
%
%   The function implements the algorithm from Riedinger (2022), "Solving Infinite-Dimensional 
%   Harmonic Lyapunov and Riccati Equations," by solving the equivalent system in the 
%   harmonic domain.
%
%   Theoretical Background:
%   The Sylvester equation in the time domain is:
%       dx + A(t) * x + x * B(t) + C(t) = 0
%   
%   The associated harmonic Sylvester equation is:
%       dX + A * X + X * B + C + N * X - X * N = 0
%   
%   For a time-periodic solution, dX = 0, and the homogeneous Sylvester equation holds:
%       (Ahm + Nh) * X + X * (Bhm - Nh) + Chm = 0
%
%   Alternatively, it can be written as:
%       (Ahm - N*) * X + X * (Bhm - Nh) + Chm = 0
%
%   The algorithm used for solving this equation involves the following steps:
%       - Time Domain Expression: 
%           dXvec = -[ B(t)^T ⊗ A(t) ] * Xvec + Cvec
%       - Harmonic Expression:
%           dXv = -[I_nb ⊗ Ahm + I_na ⊗ Bhm^* + N_nanb] * Xv - Chmv
%           Xv =  [-I_nb ⊗ Ahm - I_na ⊗ Bhm^* - N_nanb]^-1 * Chmv
%
%   Here, the Kronecker product of I_na and Bhm^* is represented as:
%       Toep(B(t)^T ⊗ I_na)
%
%   The final system requires solving for X by matrix inversion.
%
%   Input Arguments:
%   - Ahm (3D array): The phasors of the matrix A(t), of size n x n x p1, where p1 is odd.
%   - Bhm (3D array): The phasors of the matrix B(t), of size m x m x p2.
%   - Chm (3D array): The phasors of the matrix C(t), of size n x m x p3.
%   - h (scalar): The harmonic truncation order.
%   - omeg (scalar): The pulsation (frequency) of the system.
%
%   Output Arguments:
%   - Xph (3D array): The harmonic phasors of the matrix X(t), stored along the third dimension.
%       - Xph(:,:,h+1+d) is the dth phasor of X(t), where X(t) is the time matrix associated with the harmonic matrix Xhm.
%       - The size of Xph is n x m x (2*h+1).
% - M (matrix): The matrix used in the calculation.
%   - M1, M2 (matrices): Intermediate matrices used for efficient matrix inversion.
%   - colQ (vector): The column vector of Q, extracted from Chm.
%   - colX (vector): The column vector of X, the solution to the Sylvester equation.
%
%   Example Usage:
%   % Solve the harmonic Sylvester equation for given matrices
%   Xph = Sylv_harmonique(Ahm, Bhm, Chm, h, omeg);
%
%   % Solve the harmonic Lyapunov equation for given matrices
%   Xph = Sylv_harmonique(Ahm.', Ahm, Qhm, h, omeg);
%
%   Notes:
%   - The matrices Ahm, Bhm, and Chm must be provided as 3D arrays, where the third dimension
%     contains the phasors of the time-periodic matrices.
%   - The harmonic truncation order (h) must be selected carefully based on the frequency content
%     of the matrices involved.
%
%   See also: kron, blkdiag, sparse, omeg.

        narginchk(4,6)
        % dX + A(t) X + X B(t) + Q =0
        %Sylvester (A+N) X + X (B-N) + Q = 0
        %col(dX) = (I kron A(t) + B(t)^T kron I ) col(P) + col(Q(t)) 
        %(I kron A + I o B* - N_{nAnB}) col X + colQ = 0
        
       Q = Chm;
       
       A = sparray2TBlocks(Ahm,2*h);
       B = sparray2TBlocks(Bhm,2*h);
       
       if isa(Q,'PhasorArray')
       Q=value(Q);
       end
       nxa = size(Ahm(:,:,1),1);
       nxB = size(Bhm(:,:,1),1);
       
       N = kron(speye(nxa*nxB), sparse(1:2*h+1,1:2*h+1,1i*(-h:h)'*omeg));
       

    %colQ 
       hQ=(size(Q,3)-1)/2; %troncature order of provided Q
       if hQ < h 
           dQ=padarray(Q,[0 0 h-hQ],0,'both'); %pad Q with 0 phasor to match h truncatur
       elseif hQ > h 
           dQ= Q(:,:,hQ+1+(-h:h));
       else
           dQ=Q;
       end
       
       dQ=permute(dQ,[3,1,2]);
       colQ=squeeze(reshape(dQ,[],1,1));
       % colQ=sparse(colQ);

    %more effcient I kron A by using blkdiag and repmat
        tmp = repmat({A},nxB,1);
        M1 = blkdiag(tmp{:});
        M1=sparse(M1);

    % I o B*
        M2=PR_In(B',nxa,h);
        M2=sparse(M2);

    % E = -I kron A - I o B* - N 
        M=-M1-M2-N;

    % E\colQ
        colX=(M\colQ);

    %col2Mat
       dX=reshape(full(colX),[],nxa,nxB);
       Xph=permute(dX,[2,3,1]);

end