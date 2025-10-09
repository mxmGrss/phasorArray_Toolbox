function [Xph,M,M1,M2,colQ,colX] = Lyap_Harmonique(Ahm,Bhm,h,omeg,varargin)
    %LYAP_HARMONIQUE Solve Harmonic Lyapunov, Sylvester, and Riccati equations.
    %
    %   Xph = Lyap_Harmonique(Ahm, Qhm, h, omeg) solves the Harmonic Lyapunov matrix equation:
    %       (Ahm - Nh)' * X + X * (Ahm - Nh) + Qhm = 0
    %
    %   Xph = Lyap_Harmonique(Ahm, Bhm, Chm, h, omeg) solves the Harmonic Sylvester equation:
    %       (Ahm - Nh) * X + X * (Bhm - Nh) + Chm = 0   
    %
    %   Xph = Lyap_Harmonique(Ahm, Bhm, Rhm, Qhm, h, omeg) solves the Harmonic Riccati equation:
    %       (Ahm - Nh)' * X + X * (Ahm - Nh) + X * Bhm * Rhm^-1 * Bhm' * X + Qhm = 0
    %
    %   The input matrices Ahm, Bhm, and Qhm should be provided as 3D arrays of size n x m x p,
    %   where each block represents a different harmonic component of the system, with p being odd
    %   and the central element representing the 0th harmonic.
    %
    %   Inputs:
    %       - Ahm: A 3D array representing the harmonic components of matrix A(t).
    %       - Bhm: A 3D array representing the harmonic components of matrix B(t).
    %       - h: The harmonic truncation order (positive integer).
    %       - omeg: The frequency corresponding to the diagonal elements of the matrix N.
    %       - varargin: Optional arguments specifying the type of equation to solve:
    %           - Qhm: The matrix Q for Lyapunov and Sylvester equations.
    %           - Rhm: The matrix R for the Riccati equation (in addition to Qhm).
    %
    %   Outputs:
    %       - Xph: The solution matrix X in phasor form, represented as a 3D array of phasors.
    %       - M, M1, M2: Intermediate matrices used in the computation.
    %       - colQ: The column vector of Q values, adjusted to the harmonic order.
    %       - colX: The column vector of the solution X values.
    %
    %   See also: kron, blkdiag, sparse, toeplitz, evalTimeCmplx

narginchk(4,6)
ni = nargin;
if ni<6 %4 or 5 argument is lyapunov or sylvester equation, common resolution
    t0=tic;
    if ni<5 %Lyapunov ( A-N )* P + P ( A-N ) + Q = 0
       Q=Bhm;
       A = array2TBlocks(Ahm,2*h);
       nxa = size(Ahm(:,:,1),1);
       nxB=nxa;
       NhA = kron(eye(nxa), diag(1i*(-h:h)*omeg));
       AmN=sparse((A-NhA)');
       B=sparse(A) ;
    elseif ni<6  %Sylvester (A-N) X + X (B-N) + Q = 0
        % <=>colX =  -(idna kron (A-N) - idn o B')^-1 colQ
       Q = -h;
       h = omeg; %Sylvester
       omeg = varargin{1}; %Sylvester
       
       A = array2TBlocks(Ahm,2*h);
       B = array2TBlocks(Bhm,2*h);

       nxa = size(Ahm(:,:,1),1);
       NhA = kron(eye(nxa), diag(1i*(-h:h)'*omeg));
       nxB = size(Bhm(:,:,1),1);
       
       
       AmN=sparse((A-NhA)');
       
    end
    tprep=toc(t0);
    Q=value(Q);
       hQ=(size(Q,3)-1)/2; %troncature order of provided Q
       if hQ < h 
           dQ=padarray(Q,[0 0 h-hQ],0,'both'); %pad Q with 0 phasor to match h truncatur
       elseif hQ > h 
           dQ= Q(:,:,hQ+(-h:h));
       else
           dQ=Q;
       end
       
       dQ=permute(dQ,[3,1,2]);
       colQ=squeeze(reshape(dQ,[],1,1));
       colQ=sparse(colQ);
    %more effcient kron nxB Amb by using blkdiag
        t1=tic;
        tmp = repmat({AmN},nxB,1);
        M1 = blkdiag(tmp{:});
        M1=sparse(M1);
        term1=toc(t1);
        t2=tic;
        M2=PR_In(B',nxa,h);
        M2=sparse(M2);
        term2=toc(t2);
        M=M1+M2;
        t3=tic;
        colX=-(M\colQ);
        Tinv=toc(t3);
        t4=tic;

       dX=reshape(full(colX),[],nxa,nxB);
       Xph=permute(dX,[2,3,1]);
       treshape=toc(t4);

else %Ricatti, TO DO
   Q= h;
   h = omeg; %Sylvester
   omeg = varargin{1}; %Sylvester
   Qhm= h;
   Rhm = varargin{1}; 
   h = varargin{2}; %Riccati
end

% lqr(a,b,q,r,)
end