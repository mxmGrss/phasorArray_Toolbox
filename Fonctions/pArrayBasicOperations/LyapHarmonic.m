function [Xph,M,M1,M2,colQ,colX] = LyapHarmonic(Ahm, Qhm, h, omega, options)
% LYAPHARMORIC Solve Harmonic Lyapunov and Riccati equations.
%
%   X = LyapHarmonic(Ahm, Qhm, h, omega) solves the Harmonic Lyapunov equation:
%       (Ahm - Nh)' * X + X * Ahm + Qhm = 0
%
%   This matches the internal logic of the legacy Lyap_Harmonique.m.
%
%   Inputs:
%       - Ahm, Qhm: 3D arrays of harmonic components.
%       - h: Harmonic truncation order.
%       - omega: Fundamental frequency.

    arguments
        Ahm
        Qhm
        h (1,1) {mustBeInteger, mustBePositive}
        omega (1,1) {mustBeNumeric}
        options.B = []
        options.R = []
    end

    if isempty(options.B)
        % Lyapunov Mode
        A_tb = array2TBlocks(Ahm, 2*h);
        nxA = size(Ahm, 1);
        
        % Nh = kron(I, diag(jk*omega))
        NhA = kron(eye(nxA), diag(1i*(-h:h)*omega));
        
        % M1 = I ⊗ (A - Nh)'
        AmN = sparse((A_tb - NhA)');
        tmp = repmat({AmN}, nxA, 1);
        M1 = sparse(blkdiag(tmp{:}));
        
        % M2 = A' ⊗ I
        M2 = sparse(PR_In(A_tb', nxA, h));
        
        M = M1 + M2;

        % Handle Q
        Q = Qhm;
        if isa(Q, 'PhasorArray'), Q = Q.value; end
        hQ = (size(Q, 3) - 1) / 2;
        if hQ < h
            dQ = padarray(Q, [0 0 h-hQ], 0, 'both');
        elseif hQ > h
            dQ = Q(:, :, hQ+1+(-h:h));
        else
            dQ = Q;
        end
        dQ = permute(dQ, [3, 1, 2]);
        colQ = sparse(reshape(dQ, [], 1));

        % Solve
        colX = -(M \ colQ);
        dX = reshape(full(colX), [], nxA, nxA);
        Xph = permute(dX, [2, 3, 1]);
    else
        error('LyapHarmonic:NotImplemented', 'Riccati solver is not yet implemented.');
    end
end
