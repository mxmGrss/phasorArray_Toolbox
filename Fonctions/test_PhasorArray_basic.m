function results = test_PhasorArray_basic()
%TEST_PHASORARRAY_BASIC Basic test suite for the PhasorArray Toolbox.
%
%   results = TEST_PHASORARRAY_BASIC() runs a comprehensive set of tests
%   covering constructors, arithmetic, indexing, concatenation, energy,
%   derivatives, and the Toeplitz/Fourier representations.
%
%   The function returns a struct array with fields:
%       .name    - Test name (string)
%       .passed  - Logical (true if passed)
%       .message - Error message if failed, '' if passed
%
%   USAGE:
%       results = test_PhasorArray_basic();
%       % Display summary
%
%   REQUIREMENTS:
%       - PhasorArray Toolbox installed and on the MATLAB path
%       - No additional toolboxes required (YALMIP tests are skipped if absent)
%
%   See also: PHASORARRAY, INSTALLALLTOOLBOX.

fprintf('\n========================================\n');
fprintf('  PhasorArray Toolbox — Basic Test Suite\n');
fprintf('========================================\n\n');

results = struct('name', {}, 'passed', {}, 'message', {});
tol = 1e-10;   % numerical tolerance for floating-point comparisons

% In tests we intentionally exercise inv(); hide advisory noise.
warnId = 'phasorArray:PhasorInv:useHmcDivide';
warnState = warning('query', warnId);
cleanupWarn = onCleanup(@() warning(warnState.state, warnId)); %#ok<NASGU>
warning('off', warnId);

%% ========================================================================
%  1. CONSTRUCTORS
%  ========================================================================
results(end+1) = runTest('Constructor: 3D array', @() test_constructor_3D(tol));
results(end+1) = runTest('Constructor: zero matrix', @() test_constructor_zero());
results(end+1) = runTest('Constructor: eye', @() test_constructor_eye(tol));
results(end+1) = runTest('Constructor: random', @() test_constructor_random());
results(end+1) = runTest('Constructor: cos / sin', @() test_constructor_trig(tol));
results(end+1) = runTest('Constructor: isreal flag', @() test_constructor_isreal(tol));
results(end+1) = runTest('Constructor: ScalarPhasorArray', @() test_constructor_scalar(tol));

%% ========================================================================
%  2. PROPERTIES & SIZE
%  ========================================================================
results(end+1) = runTest('Properties: size / h / dim', @() test_properties());

%% ========================================================================
%  3. INDEXING
%  ========================================================================
results(end+1) = runTest('Indexing: parentheses ()', @() test_indexing_parens(tol));
results(end+1) = runTest('Indexing: braces {} read', @() test_indexing_braces_read(tol));
results(end+1) = runTest('Indexing: braces {} absent harmonic', @() test_indexing_braces_absent(tol));
results(end+1) = runTest('Indexing: braces {} write', @() test_indexing_braces_write(tol));

%% ========================================================================
%  4. ARITHMETIC
%  ========================================================================
results(end+1) = runTest('Arithmetic: addition', @() test_add(tol));
results(end+1) = runTest('Arithmetic: addition (different h)', @() test_add_pad(tol));
results(end+1) = runTest('Arithmetic: multiplication (Cauchy)', @() test_mul(tol));
results(end+1) = runTest('Arithmetic: scalar multiply', @() test_scalar_mul(tol));
results(end+1) = runTest('Arithmetic: inverse', @() test_inv(tol));
results(end+1) = runTest('Arithmetic: A\\B vs inv(A)*B (square, vector+matrix)', @() test_left_divide_vs_inv_mul_square());
results(end+1) = runTest('Arithmetic: B/A vs B*inv(A) (square, vector+matrix)', @() test_right_divide_vs_mul_inv_square());
results(end+1) = runTest('Arithmetic: A\\B vs inv(A)*B (rectangular, vector+matrix)', @() test_left_divide_vs_inv_mul_rect());
results(end+1) = runTest('Arithmetic: B/A vs B*inv(A) (rectangular, vector+matrix)', @() test_right_divide_vs_mul_inv_rect());
results(end+1) = runTest('Arithmetic: transpose', @() test_transpose(tol));
results(end+1) = runTest('Arithmetic: Hermitian', @() test_hermitian(tol));

%% ========================================================================
%  5. DERIVATIVE & ANTI-DERIVATIVE
%  ========================================================================
results(end+1) = runTest('Calculus: derivative d()', @() test_derivative(tol));
results(end+1) = runTest('Calculus: antiD() + d() roundtrip', @() test_antiD_roundtrip(tol));
results(end+1) = runTest('Calculus: antiD() DC forced to zero', @() test_antiD_DC(tol));

%% ========================================================================
%  6. CONCATENATION
%  ========================================================================
results(end+1) = runTest('Concatenation: horzcat', @() test_horzcat(tol));
results(end+1) = runTest('Concatenation: vertcat', @() test_vertcat(tol));
results(end+1) = runTest('Concatenation: blkdiag', @() test_blkdiag(tol));

%% ========================================================================
%  7. ENERGY & DIAGNOSTICS
%  ========================================================================
results(end+1) = runTest('Energy: Parseval (total = DC + AC)', @() test_energy_parseval(tol));
results(end+1) = runTest('Energy: real + imag = total', @() test_energy_decomp(tol));
results(end+1) = runTest('Diagnostics: isreal for real signal', @() test_isreal());

%% ========================================================================
%  8. TOEPLITZ & OPERATOR FUNCTIONS
%  ========================================================================
results(end+1) = runTest('Toeplitz: T_tb round-trip', @() test_T_tb_roundtrip(tol));
results(end+1) = runTest('Toeplitz: T_tb multiplication consistency', @() test_T_tb_mul(tol));
results(end+1) = runTest('Operator: N_tb structure', @() test_N_tb(tol));

%% ========================================================================
%  9. EVAL & RECONSTRUCTION
%  ========================================================================
results(end+1) = runTest('Eval: evalp consistency', @() test_evalp(tol));

%% ========================================================================
%  10. DETERMINANT
%  ========================================================================
results(end+1) = runTest('Det: det() runs without error', @() test_det_runs());
results(end+1) = runTest('Det: detLeibnizHmc 2x2', @() test_detLeibniz_2x2(tol));
results(end+1) = runTest('Det: detLeibnizHmc 3x3', @() test_detLeibniz_3x3(tol));

%% ========================================================================
%  11. FREQUENCY BASE
%  ========================================================================
results(end+1) = runTest('Base: expandBase / squishBase round-trip', @() test_expandSquish(tol));

%% ========================================================================
%  SUMMARY
%  ========================================================================
nPass = sum([results.passed]);
nFail = sum(~[results.passed]);
nTotal = numel(results);

fprintf('\n========================================\n');
fprintf('  RESULTS: %d/%d passed', nPass, nTotal);
if nFail > 0
    fprintf(', %d FAILED', nFail);
end
fprintf('\n========================================\n\n');

if nFail > 0
    fprintf('FAILED TESTS:\n');
    for k = 1:nTotal
        if ~results(k).passed
            fprintf('  ✗ %s: %s\n', results(k).name, results(k).message);
        end
    end
    fprintf('\n');
end

end

%% ========================================================================
%  TEST RUNNER
%  ========================================================================
function result = runTest(name, testFcn)
    result.name = name;
    try
        testFcn();
        result.passed = true;
        result.message = '';
        fprintf('  ✓ %s\n', name);
    catch ME
        result.passed = false;
        result.message = ME.message;
        fprintf('  ✗ %s — %s\n', name, ME.message);
    end
end

%% ========================================================================
%  1. CONSTRUCTOR TESTS
%  ========================================================================
function test_constructor_3D(tol)
    V = randn(3, 3, 11);   % h = 5
    A = PhasorArray(V);
    assert(isequal(size(A), [3, 3, 11]), 'Size mismatch');
    assert(A.h == 5, 'Harmonic order mismatch');
    Aval = A.value;
    assert(max(abs(Aval(:) - V(:))) < tol, 'Value mismatch');
end

function test_constructor_zero()
    A = PhasorArray.zeros(3, 4);
    assert(isequal(size(A, 1), 3), 'Row count');
    assert(isequal(size(A, 2), 4), 'Col count');
    assert(A.h == 0, 'Zero matrix should have h=0');
end

function test_constructor_eye(tol)
    A = PhasorArray.eye(3);
    assert(isequal(size(A,1), 3), 'Size mismatch');
    val = A{1:3, 1:3, 0};  % DC component
    err = max(abs(val - eye(3)), [], 'all');
    assert(err < tol, 'DC should be identity');
end

function test_constructor_random()
    A = PhasorArray.random(4, 3, 7);
    assert(isequal(size(A,1), 4), 'Row mismatch');
    assert(isequal(size(A,2), 3), 'Col mismatch');
    assert(A.h == 7, 'Harmonic order mismatch');
end

function test_constructor_trig(tol)
    c = PhasorArray.cos();
    s = PhasorArray.sin();
    assert(c.h == 1, 'cos should have h=1');
    assert(s.h == 1, 'sin should have h=1');
    % Evaluate at t=0: cos(0)=1, sin(0)=0
    c0 = evalp(c, 0);
    s0 = evalp(s, 0);
    assert(abs(c0 - 1) < tol, 'cos(0) should be 1');
    assert(abs(s0) < tol, 'sin(0) should be 0');
end

function test_constructor_isreal(tol)
    V = randn(2, 2, 4);  % h will be interpreted as positive-only phasors
    A = PhasorArray(V, 'isreal', true);
    assert(isreal(A), 'Should be flagged as real');
    % Evaluate: should produce real values
    At = evalp(A, pi/4);
    assert(max(abs(imag(At(:)))) < tol, 'Evaluation should be real');
end

function test_constructor_scalar(tol)
    v = [0.5, 1, 0.5];   % DC=1, h1=0.5
    A = ScalarPhasorArray(v);
    assert(A.h == 1, 'h should be 1');
    % DC coefficient should be 1
    dc = A{1,1,0};
    assert(abs(dc - 1) < tol, 'DC mismatch');
end

%% ========================================================================
%  2. PROPERTIES TESTS
%  ========================================================================
function test_properties()
    A = PhasorArray.random(3, 4, 5);
    assert(size(A, 1) == 3, 'size(A,1)');
    assert(size(A, 2) == 4, 'size(A,2)');
    assert(A.h == 5, 'h');
    assert(dim(A) == 12, 'dim = n1*n2');
end

%% ========================================================================
%  3. INDEXING TESTS
%  ========================================================================
function test_indexing_parens(tol)
    V = randn(3, 3, 11);
    A = PhasorArray(V);
    % Parentheses access the raw 3D array
    slice = A(1, 1, 6);   % (1,1, h+1) = DC for h=5
    assert(abs(slice - V(1, 1, 6)) < tol, 'Parentheses indexing mismatch');
end

function test_indexing_braces_read(tol)
    V = randn(3, 3, 11);
    A = PhasorArray(V);
    % Brace access uses harmonic order: k=0 -> DC -> slice index h+1=6
    dc = A{1, 1, 0};
    assert(abs(dc - V(1, 1, 6)) < tol, 'Brace DC mismatch');
    % k=1 -> slice index h+2=7
    h1 = A{1, 1, 1};
    assert(abs(h1 - V(1, 1, 7)) < tol, 'Brace h=1 mismatch');
end

function test_indexing_braces_absent(tol)
    A = PhasorArray.random(2, 2, 3);  % h=3
    % Request harmonic far beyond stored range
    val = A{1, 1, 100};
    assert(abs(val) < tol, 'Absent harmonic should return 0');
end

function test_indexing_braces_write(tol)
    A = PhasorArray.zeros(2, 2);   % h=0
    A{1, 1, 2} = 5;               % assign to harmonic k=2 (auto-expands)
    assert(A.h >= 2, 'h should have grown');
    assert(abs(A{1, 1, 2} - 5) < tol, 'Assignment mismatch');
end

%% ========================================================================
%  4. ARITHMETIC TESTS
%  ========================================================================
function test_add(tol)
    A = PhasorArray.random(2, 2, 3);
    B = PhasorArray.random(2, 2, 3);
    C = A + B;
    % Check DC: C{:,:,0} = A{:,:,0} + B{:,:,0}
    assert(max(abs(C{:,:,0} - A{:,:,0} - B{:,:,0}), [], 'all') < tol, 'DC addition');
end

function test_add_pad(tol)
    A = PhasorArray.random(2, 2, 3);  % h=3
    B = PhasorArray.random(2, 2, 5);  % h=5
    C = A + B;
    assert(C.h == 5, 'Result should have h=max(3,5)=5');
    % A's h=4,5 harmonics are zero, so C{:,:,4} = B{:,:,4}
    assert(max(abs(C{:,:,4} - B{:,:,4}), [], 'all') < tol, 'Padding mismatch');
end

function test_mul(tol)
    A = PhasorArray.random(2, 2, 2);  % h=2
    B = PhasorArray.random(2, 2, 3);  % h=3
    C = A * B;
    assert(C.h == 5, 'h_C should be h_A + h_B = 5');
    % Verify via time-domain sampling
    t = linspace(0, 2*pi, 200);
    for k = 1:5
        tk = t(k*30);
        At = evalp(A, tk);
        Bt = evalp(B, tk);
        Ct = evalp(C, tk);
        err = max(abs(At * Bt - Ct), [], 'all');
        assert(err < 1e-6, sprintf('Time-domain verification failed at t=%g', tk));
    end
end

function test_scalar_mul(tol)
    A = PhasorArray.random(2, 2, 3);
    C = 3 * A;
    assert(max(abs(C{:,:,0} - 3 * A{:,:,0}), [], 'all') < tol, 'Scalar multiply DC');
end

function test_inv(tol)
    A = PhasorArray.eye(2) + 0.1 * PhasorArray.random(2, 2, 2);
    Ainv = inv(A);
    I_approx = A * Ainv;  % should be ≈ identity
    % Check DC of result is close to I
    I_dc = I_approx{:,:,0};
    assert(max(abs(I_dc - eye(2)), [], 'all') < 0.05, 'A*inv(A) DC should be ~I');
end

function test_left_divide_vs_inv_mul_square()
    [A, Bsquare, BvecLeft] = make_pd_case_3x3_h15();
    Brect = PhasorArray.random(3, 5, 15); % n x m

    % A square (n x n), b vector (n x 1)
    check_left_divide_equivalence(A, BvecLeft, 3e-3, 3e-3, 'square: b vector (n x 1)', true);
    % A square (n x n), b square (n x n)
    check_left_divide_equivalence(A, Bsquare, 3e-3, 3e-3, 'square: b square (n x n)', true);
    % A square (n x n), b rectangle (n x m)
    check_left_divide_equivalence(A, Brect, 3e-3, 3e-3, 'square: b rectangle (n x m)', true);
end

function test_right_divide_vs_mul_inv_square()
    [A, Bsquare, ~, BvecRight] = make_pd_case_3x3_h15();
    Brect = PhasorArray.random(5, 3, 15); % l x n

    % xA = b with A square (n x n), b vector (1 x n)
    check_right_divide_equivalence(BvecRight, A, 3e-3, 3e-3, 'square: b vector (1 x n)', true);
    % xA = b with A square (n x n), b square (n x n)
    check_right_divide_equivalence(Bsquare, A, 3e-3, 3e-3, 'square: b square (n x n)', true);
    % xA = b with A square (n x n), b rectangle (l x n)
    check_right_divide_equivalence(Brect, A, 3e-3, 3e-3, 'square: b rectangle (l x n)', true);
end

function test_left_divide_vs_inv_mul_rect()
    [A, BmatLeft, BvecLeft] = make_rect_case_4x3_h15();

    % A rectangle (n x m), b vector (n x 1)
    check_left_divide_equivalence(A, BvecLeft, 1e-2, 1e-2, 'rectangular: b vector (n x 1)', false);
    % A rectangle (n x m), b rectangle (n x l)
    check_left_divide_equivalence(A, BmatLeft, 1e-2, 1e-2, 'rectangular: b rectangle (n x l)', false);
end

function test_right_divide_vs_mul_inv_rect()
    [A, ~, ~, BmatRight, BvecRight] = make_rect_case_4x3_h15();

    % xA = b with A rectangle (n x m), b vector (1 x m)
    check_right_divide_equivalence(BvecRight, A, 1e-2, 1e-2, 'rectangular: b vector (1 x m)', false);
    % xA = b with A rectangle (n x m), b rectangle (l x m)
    check_right_divide_equivalence(BmatRight, A, 1e-2, 1e-2, 'rectangular: b rectangle (l x m)', false);
end

function test_transpose(tol)
    A = PhasorArray.random(3, 2, 5);
    At = A.';
    assert(size(At, 1) == 2 && size(At, 2) == 3, 'Transposed size');
    % Check DC transpose
    assert(max(abs(At{:,:,0} - A{:,:,0}.'), [], 'all') < tol, 'DC transpose');
end

function test_hermitian(tol)
    A = PhasorArray.random(3, 2, 5);
    Ah = A';
    assert(size(Ah, 1) == 2 && size(Ah, 2) == 3, 'Hermitian size');
    % For Hermitian: (A')_k = conj(A_{-k})'
    % At k=0: conj(A_0)' 
    assert(max(abs(Ah{:,:,0} - conj(A{:,:,0}).'), [], 'all') < tol, 'Hermitian DC');
end

%% ========================================================================
%  5. CALCULUS TESTS
%  ========================================================================
function test_derivative(tol)
    % d/dt of cos(ωt) = -ω sin(ωt)
    c = PhasorArray.cos();
    T = 2*pi;
    dc = d(c, T);
    % Evaluate at t=0: d/dt cos(t)|_0 = -sin(0) = 0
    dc0 = evalp(dc, 0);
    assert(abs(dc0) < tol, 'd/dt cos(0) should be 0');
    % Evaluate at t=pi/2: d/dt cos(t)|_{pi/2} = -sin(pi/2) = -1
    dcpi2 = evalp(dc, pi/2);
    assert(abs(dcpi2 - (-1)) < tol, 'd/dt cos(pi/2) should be -1');
end

function test_antiD_roundtrip(tol)
    c = PhasorArray.cos();
    T = 2*pi;
    % antiD(d(cos)) should give back cos (up to DC)
    dc = d(c, T);
    adc = antiD(dc, T);
    % Compare AC part: should be cos(t) (DC may differ)
    err_h1 = abs(adc{1,1,1} - c{1,1,1});
    assert(err_h1 < tol, 'antiD(d(cos)) h=1 mismatch');
end

function test_antiD_DC(tol)
    A = PhasorArray.random(2, 2, 3);
    intA = antiD(A, 2*pi);
    dc = intA{:,:,0};
    assert(max(abs(dc(:))) < tol, 'antiD DC should be forced to zero');
end

%% ========================================================================
%  6. CONCATENATION TESTS
%  ========================================================================
function test_horzcat(tol)
    A = PhasorArray.random(2, 2, 3);
    B = PhasorArray.random(2, 3, 5);
    C = [A, B];
    assert(size(C, 1) == 2, 'Rows');
    assert(size(C, 2) == 5, 'Cols = 2+3');
    assert(C.h == 5, 'h = max(3,5)');
end

function test_vertcat(tol)
    A = PhasorArray.random(2, 3, 3);
    B = PhasorArray.random(4, 3, 5);
    C = [A; B];
    assert(size(C, 1) == 6, 'Rows = 2+4');
    assert(size(C, 2) == 3, 'Cols');
    assert(C.h == 5, 'h = max(3,5)');
end

function test_blkdiag(tol)
    A = PhasorArray.random(2, 2, 3);
    B = PhasorArray.random(3, 3, 5);
    C = blkdiag(A, B);
    assert(size(C, 1) == 5, 'Rows = 2+3');
    assert(size(C, 2) == 5, 'Cols = 2+3');
end

%% ========================================================================
%  7. ENERGY & DIAGNOSTICS TESTS
%  ========================================================================
function test_energy_parseval(tol)
    A = PhasorArray.random(3, 3, 5);
    E = energy(A);
    Edc = DCenergy(A);
    Eac = ACenergy(A);
    assert(abs(E - Edc - Eac) < tol * E, 'Parseval: total ≠ DC + AC');
end

function test_energy_decomp(tol)
    A = PhasorArray.random(3, 3, 5);
    E = energy(A);
    Er = realEnergy(A);
    Ei = imagEnergy(A);
    assert(abs(E - Er - Ei) < tol * max(E, 1), 'Energy: total ≠ real + imag');
end

function test_isreal()
    v = randn(2, 2, 4);
    A = PhasorArray(v, 'isreal', true);
    assert(isreal(A), 'Should be real');
end

%% ========================================================================
%  8. TOEPLITZ TESTS
%  ========================================================================
function test_T_tb_roundtrip(tol)
    A = PhasorArray.random(2, 2, 3);
    hTrunc = 10;
    T = A.T_tb(hTrunc);
    % Reconstruct: the central column blocks should contain A's phasors
    n = size(A, 1);
    Nt = 2*hTrunc + 1;
    expectedSize = n * Nt;
    assert(isequal(size(T), [expectedSize, expectedSize]), 'T_tb size mismatch');
end

function test_T_tb_mul(tol)
    % A*B in PhasorArray should be consistent with T_tb(A) * F_tb(B)
    A = PhasorArray.random(2, 2, 2);
    B = PhasorArray.random(2, 1, 2);
    C = A * B;
    
    hTrunc = 10;
    T_A = A.T_tb(hTrunc);
    F_B = B.F_tb(hTrunc);
    F_C_expected = T_A * F_B;
    F_C_actual = C.F_tb(hTrunc);
    
    err = max(abs(F_C_expected - F_C_actual));
    assert(err < 1e-8, 'T_tb * F_tb inconsistency');
end

function test_N_tb(tol)
    n = 3;
    h = 5;
    T = 2*pi;
    N = N_tb(n, h, T);
    expectedSize = n * (2*h + 1);
    assert(isequal(size(N), [expectedSize, expectedSize]), 'N_tb size');
    % In TB formalism: N = I_n ⊗ diag(jkω)
    % The k=0 block is at row/col indices (h*n+1):(h*n+n) in each n×n sub-block
    % But in TB formalism, the diagonal is organized as n blocks of (2h+1)×(2h+1)
    % The (h+1)-th diagonal entry in each block corresponds to k=0 and should be 0
    % Simpler check: the diagonal of N at positions corresponding to k=0 should be 0
    Ndiag = diag(N);
    % In TB: I_n ⊗ diag(jkω), so diagonal = repmat(jkω vector, n, 1)
    kvals = -h:h;
    omega = 2*pi/T;
    expected_diag = repmat(1i * kvals(:) * omega, n, 1);
    assert(max(abs(Ndiag - expected_diag)) < tol, 'N_tb diagonal mismatch');
end

%% ========================================================================
%  9. EVAL TESTS
%  ========================================================================
function test_evalp(tol)
    c = PhasorArray.cos();
    % cos(0) = 1
    assert(abs(evalp(c, 0) - 1) < tol, 'cos(0)');
    % cos(pi) = -1
    assert(abs(evalp(c, pi) - (-1)) < tol, 'cos(pi)');
    % cos(pi/2) = 0
    assert(abs(evalp(c, pi/2)) < tol, 'cos(pi/2)');
end

%% ========================================================================
%  10. DETERMINANT TESTS
%  ========================================================================
function test_det_runs()
    A = PhasorArray.eye(2) + 0.1 * PhasorArray.random(2, 2, 2);
    d = det(A);
    assert(isa(d, 'PhasorArray'), 'det should return PhasorArray');
end

function test_detLeibniz_2x2(tol)
    A = PhasorArray.random(2, 2, 3);
    dL = detLeibnizHmc(A);
    % Compare with formula: A11*A22 - A12*A21
    dManual = A{1,1} * A{2,2} - A{1,2} * A{2,1};
    err = energy(dL - dManual);
    assert(err < tol, sprintf('2x2 Leibniz error: %e', err));
end

function test_detLeibniz_3x3(tol)
    A = PhasorArray.random(3, 3, 2);
    dL = detLeibnizHmc(A);
    % Compare with time-domain evaluation
    t_test = [0, pi/3, pi/2, pi, 3*pi/2];
    for k = 1:numel(t_test)
        At = evalp(A, t_test(k));
        dLt = evalp(dL, t_test(k));
        dRef = det(At);
        err = abs(dLt - dRef);
        assert(err < 1e-6, sprintf('3x3 Leibniz mismatch at t=%g: err=%e', t_test(k), err));
    end
end

%% ========================================================================
%  11. FREQUENCY BASE TESTS
%  ========================================================================
function test_expandSquish(tol)
    A = PhasorArray.random(2, 2, 3);
    m = 3;
    B = expandBase(A, m);
    assert(B.h == m * A.h, 'expandBase h mismatch');
    
    C = squishBase(B, m);
    assert(C.h == A.h, 'squishBase h mismatch');
    
    % Round-trip should be lossless
    err = energy(A - C);
    assert(err < tol, sprintf('expandBase/squishBase round-trip error: %e', err));
end

function [A, Bmat, BvecLeft, BvecRight] = make_pd_case_3x3_h15()
% Build a 3x3 periodic SPD-like matrix with h=15 by rejection on det(A(t)).
    n = 3;
    h = 15;
    tGrid = linspace(0, 2*pi, 61);
    maxTries = 30;

    for k = 1:maxTries
        R = PhasorArray.random(n, n, h);
        S = 0.2 * (0.5 / n) * (R + R');          % Hermitian periodic perturbation
        A = PhasorArray.eye(n) + S;
        Bmat = PhasorArray.random(n, n, h);
        BvecLeft = PhasorArray.random(n, 1, h);
        BvecRight = PhasorArray.random(1, n, h);

        detA = det(A);
        detVals = zeros(size(tGrid));
        minSv = inf;
        for i = 1:numel(tGrid)
            detVals(i) = real(evalp(detA, tGrid(i)));
            s = svd(evalp(A, tGrid(i)));
            minSv = min(minSv, min(s));
        end

        if all(isfinite(detVals)) && all(detVals > 1e-2) && minSv > 0.2
            return
        end
    end

    error('Could not generate a valid 3x3 periodic matrix A with strictly positive det(A(t)).');
end

function [A, BmatLeft, BvecLeft, BmatRight, BvecRight] = make_rect_case_4x3_h15()
% Build a well-conditioned periodic rectangular matrix A (4x3), h=15.
    nrow = 4;
    ncol = 3;
    h = 15;
    tGrid = linspace(0, 2*pi, 61);
    maxTries = 120;

    A0 = [eye(ncol); zeros(nrow - ncol, ncol)];  % full column-rank baseline

    for k = 1:maxTries
        R = PhasorArray.random(nrow, ncol, h);
        A = PhasorArray(A0) + 0.02 * R;

        minSv = inf;
        for i = 1:numel(tGrid)
            s = svd(evalp(A, tGrid(i)));
            minSv = min(minSv, min(s));
        end

        if isfinite(minSv) && minSv > 0.2
            BmatLeft = PhasorArray.random(nrow, 5, h);     % for A\B (n x l)
            BvecLeft = PhasorArray.random(nrow, 1, h);     % for A\B
            BmatRight = PhasorArray.random(5, ncol, h);    % for B/A (l x m)
            BvecRight = PhasorArray.random(1, ncol, h);    % for B/A
            return
        end
    end

    error('Could not generate a well-conditioned rectangular periodic matrix A (4x3).');
end

function check_left_divide_equivalence(A, B, tolCons, tolRes, label, enforceSmallResidual)

    cleanupWarnings = suppress_expected_divide_warnings(~enforceSmallResidual);
    Xdiv = mldivide(A, B, "autoUpdateh", true, "verbose", 0);

    clear cleanupWarnings;
    Xinv = inv(A) * B;

    assert(size(Xdiv,1) == size(Xinv,1) && size(Xdiv,2) == size(Xinv,2), ...
        'Size mismatch in left-division comparison (%s).', label);

    eDiv = energy(A * Xdiv - B);
    eInv = energy(A * Xinv - B);
    relEnergy = abs(eDiv - eInv) / max([eDiv, eInv, 1]);
    assert(relEnergy < 1e-6, ...
        'Residual energy mismatch (%s): rel=%e, Eslash=%e, Einv=%e', ...
        label, relEnergy, eDiv, eInv);
end

function check_right_divide_equivalence(B, A, tolCons, tolRes, label, enforceSmallResidual)
    cleanupWarnings = suppress_expected_divide_warnings(~enforceSmallResidual);
    Xdiv = mrdivide(B, A, "autoUpdateh", true, "verbose", 0);
    clear cleanupWarnings;
    Xinv = B * inv(A);

    assert(size(Xdiv,1) == size(Xinv,1) && size(Xdiv,2) == size(Xinv,2), ...
        'Size mismatch in right-division comparison (%s).', label);

    eDiv = energy(Xdiv * A - B);
    eInv = energy(Xinv * A - B);
    relEnergy = abs(eDiv - eInv) / max([eDiv, eInv, 1]);
    assert(relEnergy < 1e-6, ...
        'Residual energy mismatch (%s): rel=%e, Eslash=%e, Einv=%e', ...
        label, relEnergy, eDiv, eInv);
end

function e = time_domain_consistency(X1, X2)
% Compare two PhasorArray objects in time domain when harmonic bases differ.
    tGrid = linspace(0, 2*pi, 41);
    e = 0;
    for i = 1:numel(tGrid)
        D = evalp(X1, tGrid(i)) - evalp(X2, tGrid(i));
        e = max(e, norm(D, 'fro'));
    end
end

function cleanupObj = suppress_expected_divide_warnings(shouldSuppress)
    if ~shouldSuppress
        cleanupObj = onCleanup(@()[]);
        return;
    end

    warningIds = {
        'phasorArray:mlHmcDivide:stagnation'
        'phasorArray:mlHmcDivide:residual'
        'phasorArray:mlHmcDivide:maxhReached'
        'phasorArray:mlHmcDivide:slowConvergence'
        'phasorArray:mlHmcDivide:unreachable'
    };

    prevStates = cell(size(warningIds));
    for k = 1:numel(warningIds)
        prevStates{k} = warning('query', warningIds{k});
        warning('off', warningIds{k});
    end

    cleanupObj = onCleanup(@() restore_warning_states(prevStates));
end

function restore_warning_states(prevStates)
    for k = 1:numel(prevStates)
        warning(prevStates{k}.state, prevStates{k}.identifier);
    end
end

function e = time_domain_residual_left(A, X, B)
    tGrid = linspace(0, 2*pi, 41);
    e = 0;
    for i = 1:numel(tGrid)
        D = evalp(A, tGrid(i)) * evalp(X, tGrid(i)) - evalp(B, tGrid(i));
        e = max(e, norm(D, 'fro'));
    end
end

function e = time_domain_residual_right(X, A, B)
    tGrid = linspace(0, 2*pi, 41);
    e = 0;
    for i = 1:numel(tGrid)
        D = evalp(X, tGrid(i)) * evalp(A, tGrid(i)) - evalp(B, tGrid(i));
        e = max(e, norm(D, 'fro'));
    end
end
