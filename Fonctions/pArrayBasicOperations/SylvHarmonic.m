function [Xph,M,M1,M2,colQ,colX] = SylvHarmonic(Ahmad, Bhmad, Chmad, h, omega)
% SYLVHARMONIC Solve Harmonic Sylvester Equations.
%
%   X = SylvHarmonic(Ahmad, Bhmad, Chmad, h, omega) solves:
%       (Ahmad + Nh) * X + X * Bhmad + Chmad = 0
%
%   This matches the internal logic of the legacy Sylv_harmonique.m.
%
%   Inputs:
%       - Ahmad, Bhmad, Chmad: 3D arrays of harmonic components.
%       - h: Harmonic truncation order.
%       - omega: Fundamental frequency.

    arguments
        Ahmad
        Bhmad
        Chmad
        h (1,1) double {mustBeInteger, mustBeNonnegative}
        omega (1,1) 
    end

    % 1. Extract raw values from PhasorArray objects (Bug fix: do this for all inputs)
    % This ensures all matrix operations (permute, reshape, cat) happen on double arrays.
    Ahm = Ahmad; if isa(Ahm, 'PhasorArray'), Ahm = Ahm.value; end
    Bhm = Bhmad; if isa(Bhm, 'PhasorArray'), Bhm = Bhm.value; end
    Chm = Chmad; if isa(Chm, 'PhasorArray'), Chm = Chm.value; end

    % 2. Get dimensions and handle harmonic orders
    hA = (size(Ahm, 3) - 1) / 2;
    hB = (size(Bhm, 3) - 1) / 2;
    hC = (size(Chm, 3) - 1) / 2;

    % Ensure h is at least max(hA, hB, hC) if not specified (though h is required here)
    if h < hC
        h = hC;
    end

    % 1. Convert to Tensor-Block representations (Corrected DEFINITION)
    A_tb = array2TBlocks(Ahm, 2*h);
    B_tb = array2TBlocks(Bhm, 2*h);

    % 3. Pad Chm to target h using standard MATLAB concatenation (removes padarray dependency)
    nxA = size(Ahm, 1);
    nxB = size(Bhm, 1);
    if h > hC
        padSize = h - hC;
        dQraw = cat(3, zeros(nxA, nxB, padSize), Chm, zeros(nxA, nxB, padSize));
    elseif hC > h
        dQraw = Chm(:, :, hC+1+(-h:h));
    else % h == hC
        dQraw = Chm;
    end
    dQraw = permute(dQraw, [3, 1, 2]);
    colQ = squeeze(reshape(dQraw, [], 1));

    % 3. Build Operator Matrix M to match legacy: M = -M1 - M2 - N
    % M1 = I ⊗ A_tb
    tmp = repmat({A_tb}, nxB, 1);
    M1 = sparse(blkdiag(tmp{:}));
    
    % M2 = B_tb' ⊗ I
    M2 = sparse(PR_In(B_tb', nxA, h));
    
    % N = I_spatial ⊗ diag(jk*omega)
    N_diag = 1i * (-h:h)' * omega;
    N = kron(speye(nxA * nxB), sparse(diag(N_diag)));
    
    M = -M1 - M2 - N;

    % 4. Solve and Reshape
    colX = (M \ colQ);
    dX = reshape(full(colX), [], nxA, nxB);
    Xph = permute(dX, [2, 3, 1]);
end
