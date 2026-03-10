function [Xph,M,M1,M2,colQ,colX] = SylvHarmonic(Ahmad, Bhmad, Chmad, h, omega)
% SYLVHARMONIC Solve harmonic Sylvester equations.
%
%   X = SylvHarmonic(Ahmad, Bhmad, Chmad, h, omega) solves:
%       dX/dt + Ahmad*X + X*Bhmad + Chmad = 0
%
%   Inputs:
%       - Ahmad, Bhmad, Chmad: 3D arrays of harmonic components.
%       - h: Harmonic truncation order for unknown X (hIn).
%       - omega: Fundamental pulsation.

    arguments
        Ahmad
        Bhmad
        Chmad
        h (1,1) double {mustBeInteger, mustBeNonnegative}
        omega (1,1)
    end

    % Extract raw values from PhasorArray objects if needed.
    Ahm = Ahmad; if isa(Ahm, 'PhasorArray'), Ahm = Ahm.value; end
    Bhm = Bhmad; if isa(Bhm, 'PhasorArray'), Bhm = Bhm.value; end
    Chm = Chmad; if isa(Chm, 'PhasorArray'), Chm = Chm.value; end

    % Harmonic orders of inputs.
    hA = (size(Ahm, 3) - 1) / 2;
    hB = (size(Bhm, 3) - 1) / 2;
    hC = (size(Chm, 3) - 1) / 2;

    % Backward compatibility: requested h must at least represent C.
    if h < hC
        h = hC;
    end

    % Rectangular truncation: unknown on hIn, equations enforced on hOut.
    hIn = h;
    hOut = max([hIn + hA, hIn + hB, hC]);

    % Rectangular Toeplitz-Blocks for left/right multiplication operators.
    A_tb = array2TBlocks(Ahm, [hOut, hIn]);
    B_tb = array2TBlocks(Bhm, [hOut, hIn]);

    nxA = size(Ahm, 1);
    nxB = size(Bhm, 1);

    % Pad/truncate C to equation order hOut.
    if hOut > hC
        padSize = hOut - hC;
        dQraw = cat(3, zeros(nxA, nxB, padSize), Chm, zeros(nxA, nxB, padSize));
    elseif hC > hOut
        dQraw = Chm(:, :, hC + 1 + (-hOut:hOut));
    else
        dQraw = Chm;
    end
    dQraw = permute(dQraw, [3, 1, 2]);
    colQ = squeeze(reshape(dQraw, [], 1));

    % Build operator matrix M = -M1 - M2 - N with rectangular sizes.
    tmp = repmat({A_tb}, nxB, 1);
    M1 = sparse(blkdiag(tmp{:}));

    % Right-multiplication operator for vec convention used in this file.
    % Keep exact legacy path for square truncation.
    if hOut == hIn
        M2 = sparse(PR_In(B_tb', nxA, hIn));
    else
        nRowBlock = 2*hOut + 1;
        nColBlock = 2*hIn + 1;
        nRows = nxA * nxB * nRowBlock;
        nCols = nxA * nxB * nColBlock;
        M2 = spalloc(nRows, nCols, nxA * nxB * nxB * nRowBlock * nColBlock);

        for j = 1:nxB
            for l = 1:nxB
                Bij = reshape(Bhm(l,j,:), 1, 1, []);
                Tlj = array2TBlocks(Bij, [hOut, hIn]);
                for i = 1:nxA
                    row0 = ((j-1)*nxA + (i-1)) * nRowBlock;
                    col0 = ((l-1)*nxA + (i-1)) * nColBlock;
                    M2(row0 + (1:nRowBlock), col0 + (1:nColBlock)) = ...
                        M2(row0 + (1:nRowBlock), col0 + (1:nColBlock)) + Tlj;
                end
            end
        end
    end

    % Derivative term: I_spatial kron D(hOut,hIn), where D maps X_k -> j*k*omega*X_k.
    k = (-hIn:hIn).';
    rowIdx = k + hOut + 1;
    colIdx = k + hIn + 1;
    Drect = sparse(rowIdx, colIdx, 1i * k * omega, 2*hOut + 1, 2*hIn + 1);
    N = kron(speye(nxA * nxB), Drect);

    M = -M1 - M2 - N;

    % Least-squares solve of rectangular harmonic system.
    colX = (M \ colQ);
    dX = reshape(full(colX), [], nxA, nxB);
    Xph = permute(dX, [2, 3, 1]);
end
