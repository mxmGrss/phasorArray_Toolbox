function test_RectangularToeplitz()
    % Test script for rectangular Toeplitz truncation in array2TBlocks
    
    fprintf('Running Rectangular Toeplitz Verification...\n');
    
    %% 1. Square Case (Backward Compatibility)
    h = 2;
    val = zeros(1,1,2*h+1);
    val(1,1,:) = -h:h; % Harmonics match their index
    A = PhasorArray(val);
    
    TB_square = array2TBlocks(A, h);
    expected_size = [2*h+1, 2*h+1];
    assert(all(size(TB_square) == expected_size), 'Square size mismatch');
    % Harmonic 0 should be on the diagonal
    assert(all(diag(TB_square) == 0), 'Square diagonal (harmonic 0) mismatch');
    fprintf('  ✓ Square case (h=2) passed.\n');
    
    %% 2. Rectangular Case (h1 > h2)
    h1 = 3; h2 = 1;
    TB_rect = array2TBlocks(A, [h1, h2]);
    % Rows: 2h1+1 = 7, Cols: 2h2+1 = 3
    assert(all(size(TB_rect) == [7, 3]), 'Rectangular size mismatch (h1>h2)');
    % Element (i, j) = k = i-j (centered)
    % Center row of 7 is 4 (-h1 to h1 => -3, -2, -1, 0, 1, 2, 3)
    % Center col of 3 is 2 (-h2 to h2 => -1, 0, 1)
    % TB(4, 2) should be harmonic 0 (index i=0, j=0 => k=0)
    assert(TB_rect(4, 2) == 0, 'Rectangular center (harmonic 0) mismatch');
    % TB(1, 1) should be k = -3 - (-1) = -2
    assert(TB_rect(1, 1) == -2, 'Rectangular extreme k mismatch');
    fprintf('  ✓ Rectangular case (h1=3, h2=1) passed.\n');
    
    %% 3. Rectangular Case (h1 < h2)
    h1 = 1; h2 = 2;
    TB_rect2 = array2TBlocks(A, [h1, h2]);
    assert(all(size(TB_rect2) == [3, 5]), 'Rectangular size mismatch (h1<h2)');
    assert(TB_rect2(2, 3) == 0, 'Rectangular center (harmonic 0) mismatch (h1<h2)');
    fprintf('  ✓ Rectangular case (h1=1, h2=2) passed.\n');

    %% 4. Automatic Padding
    % Current A has h=2. Requesting h1+h2=4.
    h1 = 2; h2 = 2;
    TB_pad = array2TBlocks(A, [h1, h2]); % req h_req = 4.
    assert(all(size(TB_pad) == [5, 5]), 'Square size mismatch with padding');
    % Harmonic 4 and -4 should be 0 because they were padded
    % TB(1, 5) => k = -2 - 2 = -4
    assert(TB_pad(1, 5) == 0, 'Padding verification failed');
    %% 5. Symbolic Compatibility
    if exist('sym', 'file')
        fprintf('  Checking Symbolic compatibility...\n');
        h1 = 1; h2 = 1;
        Asym = PhasorArray(sym('A', [1, 1, 3]));
        TB_sym = array2TBlocks(Asym, [h1, h2]);
        assert(isa(TB_sym, 'sym'), 'Symbolic type lost in array2TBlocks');
        assert(all(size(TB_sym) == [3, 3]), 'Symbolic size mismatch');
        fprintf('    ✓ Symbolic variables preserved.\n');
    end

    %% 6. YALMIP (sdpvar) Compatibility
    if exist('sdpvar', 'file')
        fprintf('  Checking YALMIP (sdpvar) compatibility...\n');
        h1 = 1; h2 = 1;
        % Create a 1x1x3 sdpvar array manually to avoid PhonArray.random for now
        v_sdp = sdpvar(1, 1, 3);
        TB_sdp = array2TBlocks(v_sdp, [h1, h2]);
        assert(isa(TB_sdp, 'sdpvar'), 'sdpvar type lost in array2TBlocks');
        assert(all(size(TB_sdp) == [3, 3]), 'sdpvar size mismatch');
        fprintf('    ✓ YALMIP variables preserved.\n');
    end

    fprintf('All Rectangular Toeplitz tests (including types) passed!\n');
end
