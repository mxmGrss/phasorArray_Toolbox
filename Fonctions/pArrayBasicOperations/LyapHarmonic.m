function [Xph, M, M1, M2, colQ, colX] = LyapHarmonic(Ahm, Qhm, h, omega, nvp)
%LYAPHARMONIC Solve Harmonic Lyapunov equations using Toeplitz blocks
%
%   X = LyapHarmonic(Ahm, Qhm, h, omega) solves the Continuous-Time
%   Harmonic Lyapunov equation:
%       (Ahm - Nh)' * X + X * Ahm + Qhm = 0
%   where Nh is the harmonic shift operator. The derivative operator Nh is
%   folded into the closed-loop harmonic matrix (T(A) - Nh), so this is the
%   'backward' (cost-to-go) equation:  dP/dt + A'P + PA + Q = 0.
%
%   X = LyapHarmonic(..., direction='forward') flips the sign of Nh and
%   solves dP/dt = A'P + PA + Q instead (covariance-type integration; for
%   the usual covariance equation dP/dt = AP + PA' + Q, pass A' as input).
%
%   SYNTAX:
%       X = LyapHarmonic(Ahm, Qhm, h, omega)
%       X = LyapHarmonic(Ahm, Qhm, h, omega, direction='forward')
%       [X, M, M1, M2, colQ, colX] = LyapHarmonic(...)
%
%   INPUTS:
%       Ahm     - 3D array of harmonic components for the system matrix
%       Qhm     - 3D array of harmonic components for the weight matrix
%       h       - Harmonic truncation order
%       omega   - Fundamental frequency (rad/s)
%       direction - 'backward' (default) or 'forward': sign of dP/dt
%
%   OUTPUTS:
%       Xph     - Solution 3D array (harmonic components)
%       M       - Global Kronecker-like matrix used for resolution
%       colQ    - Vectorized version of the weighted input
%
%   PhasorArray Toolbox
%   Author: Maxime GROSSO
%   Date: March 2026

    arguments
        Ahm
        Qhm
        h (1,1) {mustBeInteger, mustBeNonnegative}
        omega (1,1) {mustBeNumeric}
        nvp.B = []
        nvp.R = []
        nvp.direction {mustBeMember(nvp.direction, {'backward','forward'})} = 'backward'
    end

    if isempty(nvp.B)
        % --- Lyapunov Mode ---
        % 1. Convert to Toeplitz Block structure
        A_tb = array2TBlocks(Ahm, h);
        nxA = size(Ahm, 1);
        
        % 2. Harmonic Shift Matrix: Nh = kron(I, diag(jk*omega))
        NhA = kron(eye(nxA), diag(1i*(-h:h)*omega));

        % 3. Build LHS operator: M1 = I ⊗ (A -+ Nh)'
        % 'backward' folds -Nh (solves dP/dt + A'P + PA + Q = 0),
        % 'forward' folds +Nh (solves dP/dt = A'P + PA + Q).
        if strcmp(nvp.direction, 'forward')
            signN = +1;
        else
            signN = -1;
        end
        AmN = sparse((A_tb + signN * NhA)');
        tmp = repmat({AmN}, nxA, 1);
        M1 = sparse(blkdiag(tmp{:}));
        
        % 4. Build RHS operator: M2 = A' ⊗ I
        M2 = sparse(PR_In(A_tb', nxA, h));
        
        % Global operator
        M = M1 + M2;

        % 5. Handle Vectorized Input Q
        Q = Qhm;
        if isa(Q, 'PhasorArray'), Q = Q.value; end
        hQ = (size(Q, 3) - 1) / 2;
        if hQ < h
            dQ = phasorPad(Q, [0 0 h-hQ], 0, 'both');
        elseif hQ > h
            dQ = Q(:, :, hQ+1+(-h:h));
        else
            dQ = Q;
        end
        dQ = permute(dQ, [3, 1, 2]);
        colQ = sparse(reshape(dQ, [], 1));

        % 6. Solve Linear System
        colX = -(M \ colQ);
        dX = reshape(full(colX), [], nxA, nxA);
        Xph = permute(dX, [2, 3, 1]);
    else
        % --- Riccati Mode ---
        error('LyapHarmonic:NotImplemented', 'Riccati solver is not yet implemented.');
    end
end
