function isT = isToeplitz(A,opt)
    %ISTOEPLITZ Check if the input matrix A is Toeplitz.
    % A Toeplitz matrix has constant diagonals, i.e., 
    % A(i, j) == A(i+1, j+1) for all valid i, j.
    arguments
        A
        opt (1,1) double = 0
    end

    % Validate input
    if ~ismatrix(A)
        error('isToeplitz:invalidInput', 'Input must be a 2-D matrix; got an array of size %s.', mat2str(size(A)))
    end

    switch opt 
            case 0
        % Check if A equals toeplitz(A(:,1), A(1,:))
        isT = isequal(A, toeplitz(A(:,1), A(1,:)));
            case 1 %isequaln 
        % Check if A equals toeplitz(A(:,1), A(1,:)) ignoring NaNs
        isT = isequaln(A, toeplitz(A(:,1), A(1,:)));
            case 2 %ismembertol
            % Check if A equals toeplitz(A(:,1), A(1,:)) with tolerance
            tol = 1e-10; % Set a tolerance value
            isT = ismembertol(A, toeplitz(A(:,1), A(1,:)), tol);
        case 3 %manual tolerance
            tol = 1e-10; % Set a tolerance value
            % Check if A equals toeplitz(A(:,1), A(1,:)) with manual tolerance % Use a tolerance-based comparison:
            isT = all(abs(A - toeplitz(A(:,1), A(1,:))) < tol);
        otherwise
        error('isToeplitz:invalidOption', 'Invalid option %d. Use 0 (exact), 1 (NaN-tolerant), 2 (ismembertol), or 3 (manual tolerance).', opt)
    end


end
