function [Xph, M, M1, M2, colQ, colX] = SylvHarmonic(A, B, C, h, omega, method, direction)
%SYLVHARMONIC Solve a harmonic Sylvester equation.
%   X = SylvHarmonic(A,B,C,h,omega) solves dX/dt + AX + XB + C = 0.
%   SylvHarmonic(...,method,'forward') solves dX/dt = AX + XB + C.
%   A and B are square coefficient arrays, C has compatible spatial
%   dimensions. Inputs may also be PhasorArray objects. Harmonics are stored
%   in ascending order, with an odd number of coefficient slices.
%   omega is the fundamental angular frequency.
%
%   h sets the solution order. 'rectangle' (default) retains all product
%   and forcing rows: hOut = max(h+hA,h+hB,hC). 'square' uses hOut = h.
%   Rectangular systems are solved in the least-squares sense.
%
%   [X,M,M1,M2,colQ,colX] also returns the sparse system M*colX = colQ,
%   where M = -M1-M2-N (backward) or -M1-M2+N (forward). These outputs
%   use harmonics-fastest ordering: colX = reshape(permute(X,[3 1 2]),[],1).
%   M1 and M2 represent left and right multiplication. Right multiplication
%   uses B_k.' without conjugation, including for complex-valued signals.
%
%   See also PhasorArray/lyap.
    arguments
        A
        B
        C
        h (1,1) double {mustBeInteger, mustBeNonnegative}
        omega (1,1)
        method {mustBeMember(method, {'rectangle','square'})} = 'rectangle'
        direction {mustBeMember(direction, {'backward','forward'})} = 'backward'
    end

    if isa(A, 'PhasorArray'), A = value(A); end
    if isa(B, 'PhasorArray'), B = value(B); end
    if isa(C, 'PhasorArray'), C = value(C); end
    if ~(isnumeric(A) && isnumeric(B) && isnumeric(C))
        error('PhasorArray:SylvHarmonic:dataType', 'Numeric harmonic coefficients required.');
    end
    n = size(A, 1); m = size(B, 1);
    if size(A, 2) ~= n || size(B, 2) ~= m || size(C, 1) ~= n || size(C, 2) ~= m ...
            || any(mod([size(A,3), size(B,3), size(C,3)], 2) ~= 1)
        error('PhasorArray:SylvHarmonic:dimensions', ...
            'Square A/B, compatible C and odd coefficient counts required.');
    end

    wantParts = nargout >= 3;
    [K, q, L, R, orderOut, orderIn] = assembleSystem(A, B, C, h, omega, method, direction, wantParts);
    x = K \ q;
    Xph = reshape(full(x), n, m, 2*h+1);
    if nargout >= 2
        inverseOut = zeros(1, numel(orderOut));
        inverseIn = zeros(1, numel(orderIn));
        inverseOut(orderOut) = 1:numel(orderOut);
        inverseIn(orderIn) = 1:numel(orderIn);
        M = K(inverseOut, inverseIn);
    end
    if wantParts
        M1 = L(inverseOut, inverseIn); M2 = R(inverseOut, inverseIn);
    end
    if nargout >= 5, colQ = q(inverseOut); end
    if nargout >= 6, colX = x(inverseIn); end
end

function [K, q, L, R, orderOut, orderIn] = assembleSystem(A, B, C, h, omega, method, direction, wantParts)
    n = size(A,1); m = size(B,1); spatialSize = n*m;
    hA = (size(A,3)-1)/2; hB = (size(B,3)-1)/2; hC = (size(C,3)-1)/2;
    hOut = h;
    if strcmp(method, 'rectangle'), hOut = max([h+hA, h+hB, hC]); end
    numIn = 2*h+1; numOut = 2*hOut+1;
    maxShift = min(max(hA,hB), h+hOut);
    shifts = -maxShift:maxShift;
    entries = cell(numel(shifts)+1, 3); L = []; R = [];
    if wantParts
        leftEntries = cell(numel(shifts), 3); rightEntries = cell(numel(shifts), 3);
    end

    % In harmonic-major order, each shift contributes S_k kron
    % (I_m kron A_k + B_k.' kron I_n). Build sparse matrices once from triplets.
    for index = 1:numel(shifts)
        k = shifts(index); left = sparse(spatialSize,spatialSize); right = left;
        if abs(k) <= hA, left = kron(speye(m), sparse(A(:,:,hA+1+k))); end
        if abs(k) <= hB, right = kron(sparse(B(:,:,hB+1+k).'), speye(n)); end
        rows = (-h:h)' + k + hOut + 1;
        valid = rows >= 1 & rows <= numOut;
        shift = sparse(rows(valid), find(valid), 1, numOut, numIn);
        [entries{index,1}, entries{index,2}, entries{index,3}] = find(-kron(shift, left+right));
        if wantParts
            [leftEntries{index,1}, leftEntries{index,2}, leftEntries{index,3}] = find(kron(shift,left));
            [rightEntries{index,1}, rightEntries{index,2}, rightEntries{index,3}] = find(kron(shift,right));
        end
    end
    harmonics = (-h:h)';
    derivative = sparse(harmonics+hOut+1, (1:numIn)', 1i*harmonics*omega, numOut, numIn);
    signN = -1;
    if strcmp(direction, 'forward'), signN = 1; end
    index = numel(shifts)+1;
    [entries{index,1}, entries{index,2}, entries{index,3}] = find(signN*kron(derivative,speye(spatialSize)));
    K = sparse(vertcat(entries{:,1}), vertcat(entries{:,2}), vertcat(entries{:,3}), ...
        spatialSize*numOut, spatialSize*numIn);
    if wantParts
        L = sparse(vertcat(leftEntries{:,1}), vertcat(leftEntries{:,2}), vertcat(leftEntries{:,3}), ...
            spatialSize*numOut, spatialSize*numIn);
        R = sparse(vertcat(rightEntries{:,1}), vertcat(rightEntries{:,2}), vertcat(rightEntries{:,3}), ...
            spatialSize*numOut, spatialSize*numIn);
    end
    forcing = zeros(n,m,numOut); keep = -min(hOut,hC):min(hOut,hC);
    forcing(:,:,hOut+1+keep) = C(:,:,hC+1+keep); q = forcing(:);
    idsIn = reshape(1:spatialSize*numIn, numIn, spatialSize);
    idsOut = reshape(1:spatialSize*numOut, numOut, spatialSize);
    orderIn = reshape(idsIn.',[],1); orderOut = reshape(idsOut.',[],1);
end
