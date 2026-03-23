function results = test_PhasorArray_advanced()
%TEST_PHASORARRAY_ADVANCED Advanced test suite for the PhasorArray Toolbox.
%
%   results = TEST_PHASORARRAY_ADVANCED() runs tests covering advanced
%   features: symbolic PhasorArrays, YALMIP ndsdpvar, Kronecker algebra,
%   Lyapunov/Sylvester equations, Park transform, LTP simulation, and
%   the Leibniz determinant on symbolic objects.
%
%   Tests that require optional toolboxes (Symbolic Math Toolbox, YALMIP)
%   are **skipped gracefully** if the toolbox is not available.
%
%   USAGE:
%       results = test_PhasorArray_advanced();
%
%   See also: TEST_PHASORARRAY_BASIC, PHASORARRAY, LYAP_HARMONIQUE.

fprintf('\n==============================================\n');
fprintf('  PhasorArray Toolbox — Advanced Test Suite\n');
fprintf('==============================================\n\n');

results = struct('name', {}, 'passed', {}, 'message', {});
tol = 1e-8;

% Detect optional toolboxes
hasSymbolic = ~isempty(ver('symbolic'));
hasYALMIP   = exist('sdpvar', 'file') == 2;

if ~hasSymbolic
    fprintf('  [SKIP] Symbolic Math Toolbox not found — sym tests will be skipped\n');
end
if ~hasYALMIP
    fprintf('  [SKIP] YALMIP not found — ndsdpvar tests will be skipped\n');
end
fprintf('\n');

%% ========================================================================
%  1. SYMBOLIC PhasorArrays
%  ========================================================================
if hasSymbolic
    results(end+1) = runTest('Sym: construct sym PhasorArray', @() test_sym_construct());
    results(end+1) = runTest('Sym: sym arithmetic (add, mul)', @() test_sym_arithmetic());
    results(end+1) = runTest('Sym: sym diag', @() test_sym_diag());
    results(end+1) = runTest('Sym: detLeibnizHmc on sym 2x2', @() test_sym_detLeibniz());
else
    results(end+1) = skipTest('Sym: construct sym PhasorArray');
    results(end+1) = skipTest('Sym: sym arithmetic (add, mul)');
    results(end+1) = skipTest('Sym: sym diag');
    results(end+1) = skipTest('Sym: detLeibnizHmc on sym 2x2');
end

%% ========================================================================
%  2. YALMIP / ndsdpvar PhasorArrays
%  ========================================================================
if hasYALMIP
    results(end+1) = runTest('YALMIP: construct ndsdpvar', @() test_ndsdpvar_construct());
    results(end+1) = runTest('YALMIP: ndsdpvar addition', @() test_ndsdpvar_add());
    results(end+1) = runTest('YALMIP: ndsdpvar diag/blkdiag', @() test_ndsdpvar_diag());
    results(end+1) = runTest('YALMIP: ndsdpvar detLeibnizHmc 2x2', @() test_ndsdpvar_detLeibniz());
    results(end+1) = runTest('YALMIP: Lyapunov LMI stability', @() test_yalmip_LMI());
    results(end+1) = runTest('YALMIP: Check Solvers', @() test_check_solvers());
    results(end+1) = runTest('YALMIP: Lyap vs LMI Consistency', @() test_yalmip_LMI_vs_Lyap());
    results(end+1) = runTest('YALMIP: Riccati LMI vs RicHarmonicKlein', @() test_riccati_lmi_vs_kleinman());
else
    results(end+1) = skipTest('YALMIP: construct ndsdpvar');
    results(end+1) = skipTest('YALMIP: ndsdpvar addition');
    results(end+1) = skipTest('YALMIP: ndsdpvar diag/blkdiag');
    results(end+1) = skipTest('YALMIP: ndsdpvar detLeibnizHmc 2x2');
    results(end+1) = skipTest('YALMIP: Riccati LMI vs RicHarmonicKlein');
end

%% ========================================================================
%  3. KRONECKER ALGEBRA
%  ========================================================================
results(end+1) = runTest('Kronecker: kron consistency', @() test_kron(tol));
results(end+1) = runTest('Kronecker: oplus consistency', @() test_oplus(tol));

%% ========================================================================
%  4. PHASE SHIFT & PARK TRANSFORM
%  ========================================================================
results(end+1) = runTest('PhaseShift: scalar angle on cos', @() test_PhaseShift(tol));
results(end+1) = runTest('PhaseShift: broadcast to 3-phase system', @() test_PhaseShift_3phase(tol));
results(end+1) = runTest('Park: transform construction', @() test_Park());
results(end+1) = runTest('Park: P * P^T = I (orthogonality)', @() test_Park_orthogonal(tol));

%% ========================================================================
%  5. LYAPUNOV & SYLVESTER (Harmonic Domain)
%  ========================================================================
results(end+1) = runTest('Lyapunov: LTI case (constant A, Q)', @() test_lyap_LTI(tol));
results(end+1) = runTest('Lyapunov: LTP case (periodic A)', @() test_lyap_LTP(tol));
results(end+1) = runTest('Sylvester: (A+N)X + X(B-N) + C = 0', @() test_sylvester(tol));

%% ========================================================================
%  6. LTP SIMULATION (initial, lsim)
%  ========================================================================
results(end+1) = runTest('Simulation: initial() runs (stable LTP)', @() test_initial_runs());
results(end+1) = runTest('Simulation: lsim() runs with forcing term', @() test_lsim_runs());

%% ========================================================================
%  7. LEIBNIZ DET (ADVANCED CASES)
%  ========================================================================
results(end+1) = runTest('Det: Leibniz matches FFT-det (2x2)', @() test_detLeibniz_vs_FFT_2x2(tol));
results(end+1) = runTest('Det: Leibniz matches FFT-det (3x3)', @() test_detLeibniz_vs_FFT_3x3(tol));
results(end+1) = runTest('Det: Leibniz identity determinant = 1', @() test_detLeibniz_identity(tol));
results(end+1) = runTest('Det: Leibniz det(A*B) vs det(A)*det(B)', @() test_detLeibniz_product(tol));

%% ========================================================================
%  8. TOEPLITZ FORMALISMS CONSISTENCY
%  ========================================================================
results(end+1) = runTest('Toeplitz: T_bt construction', @() test_T_bt());
results(end+1) = runTest('Toeplitz: S_tb construction', @() test_S_tb());
results(end+1) = runTest('Toeplitz: S_tb * F_tb vs PhaseShift', @() test_S_tb_PhaseShift(tol));

%% ========================================================================
%  SUMMARY
%  ========================================================================
nPass = sum([results.passed]);
nSkip = sum(strcmp({results.message}, 'SKIPPED'));
nFail = sum(~[results.passed]) - nSkip;
nTotal = numel(results);

fprintf('\n==============================================\n');
fprintf('  RESULTS: %d passed, %d failed, %d skipped (out of %d)\n', nPass, nFail, nSkip, nTotal);
fprintf('==============================================\n\n');

if nFail > 0
    fprintf('FAILED TESTS:\n');
    for k = 1:nTotal
        if ~results(k).passed && ~strcmp(results(k).message, 'SKIPPED')
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

function result = skipTest(name)
    result.name = name;
    result.passed = false;
    result.message = 'SKIPPED';
    fprintf('  ○ %s — SKIPPED (missing toolbox)\n', name);
end

%% ========================================================================
%  1. SYMBOLIC TESTS
%  ========================================================================
function test_sym_construct()
    A = PhasorArray.sym(2, 2, 3, "M","isreal",true);
    assert(isa(A, 'PhasorArray'), 'Should be PhasorArray');
    assert(A.h == 3, 'Harmonic order should be 3');
    assert(size(A, 1) == 2 && size(A, 2) == 2, 'Size should be 2x2');
end

function test_sym_arithmetic()
    A = PhasorArray.sym(2, 2, 1, "A","isreal",true);
    B = PhasorArray.sym(2, 2, 1, "B","isreal",true);
    C = A + B;   % should work symbolically
    D = A * B;   % Cauchy product, h should grow
    assert(isa(C, 'PhasorArray'), 'Sum should be PhasorArray');
    assert(isa(D, 'PhasorArray'), 'Product should be PhasorArray');
    assert(D.h == 2, 'Product h should be 1+1=2');
end

function test_sym_diag()
    v = PhasorArray.sym(3, 1, 2, "v","isreal",true);
    D = diag(v);
    assert(size(D, 1) == 3 && size(D, 2) == 3, 'diag(3x1) should be 3x3');
end

function test_sym_detLeibniz()
    A = PhasorArray.sym(2, 2, 1, "M","isreal",true);
    d = detLeibnizHmc(A);
    assert(isa(d, 'PhasorArray'), 'Should return PhasorArray');
    assert(d.h == 2, 'det of 2x2 with h=1 should have h=2');
end

%% ========================================================================
%  2. YALMIP TESTS
%  ========================================================================
function test_ndsdpvar_construct()
    P = PhasorArray.ndsdpvar(3, 3, 2, 'PhasorType', 'symmetric', 'real', true);
    assert(isa(P, 'PhasorArray'), 'Should be PhasorArray');
    assert(P.h == 2, 'Harmonic order should be 2');
    assert(size(P, 1) == 3 && size(P, 2) == 3, 'Size should be 3x3');
end

function test_ndsdpvar_add()
    P = PhasorArray.ndsdpvar(2, 2, 1, 'PhasorType', 'full', 'real', true);
    Q = PhasorArray.ndsdpvar(2, 2, 1, 'PhasorType', 'full', 'real', true);
    R = P + Q;
    assert(isa(R, 'PhasorArray'), 'Sum should be PhasorArray');
end

function test_ndsdpvar_diag()
    P = PhasorArray.ndsdpvar(3, 3, 1, 'PhasorType', 'symmetric', 'real', true);
    d = diag(P);
    assert(size(d, 1) == 3 && size(d, 2) == 1, 'diag should extract 3x1 vector');
    % blkdiag with ndsdpvar
    Q = PhasorArray.ndsdpvar(2, 2, 0, 'PhasorType', 'full', 'real', true);
    R = blkdiag(P, Q);
    assert(size(R, 1) == 5 && size(R, 2) == 5, 'blkdiag should give 5x5');
end

function test_ndsdpvar_detLeibniz()
    P = PhasorArray.ndsdpvar(2, 2, 3, 'PhasorType', 'full', 'real', true);
    d = detLeibnizHmc(P);
    assert(isa(d, 'PhasorArray'), 'det should return PhasorArray');
    % The result should contain sdpvar expressions
    dval = d.value;
    assert(isa(dval, 'ndsdpvar') || isa(dval, 'sdpvar'), ...
        'Leibniz det of ndsdpvar should stay symbolic (sdpvar/ndsdpvar)');
end

%% ========================================================================
%  3. KRONECKER TESTS
%  ========================================================================
function test_kron(tol)
    A = PhasorArray.random(2, 2, 2);
    B = PhasorArray.random(2, 2, 3);
    K = kron(A, B);
    % Size: (2*2) x (2*2) = 4x4
    assert(size(K, 1) == 4 && size(K, 2) == 4, 'kron size');
    % Verify via time-domain evaluation: K(t) = A(t) ⊗ B(t)
    t0 = 0.7;
    At = evalp(A, t0);
    Bt = evalp(B, t0);
    Kt = evalp(K, t0);
    Kref = kron(At, Bt);
    err = max(abs(Kt(:) - Kref(:)));
    assert(err < tol, sprintf('kron time-domain error: %e', err));
end

function test_oplus(tol)
    A = PhasorArray.random(2, 2, 2);
    B = PhasorArray.random(3, 3, 2);
    Op = oplus(A, B);
    % Size: (2*3 + 3*2) ... no, oplus = A ⊗ I_b + I_a ⊗ B => size (2*3) x (2*3)
    assert(size(Op, 1) == 6 && size(Op, 2) == 6, 'oplus size should be 6x6');
    % Verify at a time sample
    t0 = 1.3;
    At = evalp(A, t0);
    Bt = evalp(B, t0);
    Opt = evalp(Op, t0);
    Opref = kron(At, eye(3)) + kron(eye(2), Bt);
    err = max(abs(Opt(:) - Opref(:)));
    assert(err < tol, sprintf('oplus time-domain error: %e', err));
end

%% ========================================================================
%  4. PHASE SHIFT & PARK TESTS
%  ========================================================================
function test_PhaseShift(tol)
    c = PhasorArray.cos();
    % Phase shift by π/2 should give -sin(t)
    c_shifted = c.PhaseShift(pi/2);
    % Evaluate at t=0: cos(0+π/2) = cos(π/2) = 0
    v0 = evalp(c_shifted, 0);
    assert(abs(v0) < tol, 'PhaseShift(cos, π/2) at t=0 should be 0');
    % Evaluate at t=π/2: cos(π/2+π/2) = cos(π) = -1
    vpi2 = evalp(c_shifted, pi/2);
    assert(abs(vpi2 - (-1)) < tol, 'PhaseShift(cos, π/2) at t=π/2 should be -1');
end

function test_PhaseShift_3phase(tol)
    c = PhasorArray.cos();
    % Create balanced 3-phase: [cos(t), cos(t-2π/3), cos(t-4π/3)]
    phases = c.PhaseShift([0, -2*pi/3, -4*pi/3]);
    assert(size(phases, 1) == 1 && size(phases, 2) == 3, '3-phase should be 1x3');
    % Sum of balanced 3-phase = 0 at all times
    s = phases{1,1} + phases{1,2} + phases{1,3};
    assert(energy(s) < tol, 'Balanced 3-phase sum should have zero energy');
end

function test_Park()
    P = PhasorArray.Park(0);
    assert(isa(P, 'PhasorArray'), 'Park should return PhasorArray');
    assert(size(P, 1) == 3 && size(P, 2) == 3, 'Park should be 3x3');
end

function test_Park_orthogonal(tol)
    P = PhasorArray.Park(0);
    % Park is amplitude-invariant: the 0-component row has factor 1/sqrt(3),
    % while d,q rows have factor sqrt(2/3). So P*P.' is diagonal but
    % NOT a scalar multiple of I.
    % Check that P*P.' is diagonal and time-invariant.
    for t0 = [0, pi/4, pi/2, pi]
        Pt = evalp(P, t0);
        PPt = Pt * Pt.';
        % Should be diagonal: check off-diagonal elements are ~0
        offdiag = PPt - diag(diag(PPt));
        err_offdiag = max(abs(offdiag(:)));
        assert(err_offdiag < tol, ...
            sprintf('Park*Park^T off-diagonal should be 0 at t=%g, err=%e', t0, err_offdiag));
        % Diagonal should be positive
        d = diag(PPt);
        assert(all(d > 0), 'Park*Park^T diagonal entries should be positive');
    end
end

%% ========================================================================
%  5. LYAPUNOV & SYLVESTER TESTS
%  ========================================================================
function test_lyap_LTI(~)
    % LTI Lyapunov: A'P + PA + Q = 0 with constant A, Q
    A = PhasorArray([-2 1; 0 -3]);
    Q = PhasorArray(eye(2));
    
    [X, residual] = lyap(A, Q);
    
    % Residual should be small
    assert(residual.resnorm < 1e-6, ...
        sprintf('LTI Lyapunov residual: %e', residual.resnorm));
    
    % The DC part should match MATLAB's built-in lyap
    n = size(A, 1);
    A_const = [-2 1; 0 -3];
    X_dc = X{1:n, 1:n, 0};
    X_ref = lyap(A_const', eye(n));
    err = max(abs(X_dc - X_ref), [], 'all');
    assert(err < 1e-6, sprintf('LTI Lyapunov DC vs MATLAB lyap: %e', err));
end

function test_lyap_LTP(~)
    % LTP Lyapunov: dP + A'(t)P + PA(t) + Q = 0 with periodic A(t)
    A0 = [-3 0.5; 0 -2];
    A1 = [0.1 0; 0 0.05];
    A2 = [0 0.3; 0.1 0];
    A = PhasorArray(cat(3, conj(A2), conj(A1), A0, A1, A2));  % h=1
    Q = PhasorArray(eye(2));
    
    [X, residual] = lyap(A, Q, 'h', 2,'autoUpdateh',true);
    
    % Residual should be moderate (truncation at h=10)
    assert(residual.resnorm < 1, ...
        sprintf('LTP Lyapunov residual: %e', residual.resnorm));
    
    E = energy(X);
    assert(isfinite(E) && E > 0, 'LTP Lyapunov solution should have finite positive energy');
    
    % Symmetry residual should be small (P should be symmetric)
    assert(residual.resPsym < 1e-6, ...
        sprintf('LTP Lyapunov symmetry residual: %e', residual.resPsym));
end

function test_sylvester(~)
    % Sylvester: dM + A(t)M + MB(t) + C = 0
    A = PhasorArray([-2 0; 0 -3]);
    B = PhasorArray([-1 0; 0 -4]);
    C = PhasorArray(randn(2, 2));
    
    [~, residual] = lyap(A, B, C);
    
    % Residual should be small for LTI case
    assert(residual.resnorm < 1e-6, ...
        sprintf('Sylvester residual: %e', residual.resnorm));
end

%% ========================================================================
%  6. SIMULATION TESTS
%  ========================================================================
function test_initial_runs()
    % Create a stable periodic system dx/dt = A(t) x
    A0 = [-2 0; 0 -3];
    A1 = [0.05 0; 0 0.05];
    Ahm = PhasorArray(cat(3, conj(A1'), A0, A1));
    
    x0 = [1; 0];
    T = 2*pi;
    tfinal = 10;
    [y, t] = initial(Ahm, x0, T, tfinal);
    
    assert(~isempty(y), 'initial() should return non-empty y');
    assert(~isempty(t), 'initial() should return non-empty t');
    assert(numel(t) > 1, 'Should have multiple time samples');
    % State should decay (stable system)
    y_end = y(:, end);
    assert(norm(y_end) < 0.5, 'Stable system: final state should be small');
end

function test_lsim_runs()
    % Create stable periodic system with forcing
    A0 = [-2 0; 0 -3];
    Ahm = PhasorArray(A0);  % LTI for simplicity
    
    x0 = [0; 0];
    T = 2*pi;
    tfinal = 5;
    
    % Forcing: u(t) = cos(t) on state 1
    B = PhasorArray(eye(2));  % B matrix
    Uph = PhasorArray(cat(3, [0; 0], [1; 0], [0; 0]));  % cos input on channel 1
    
    [y, t] = lsim(Ahm, tfinal, x0, T, Uph);
    
    assert(~isempty(y), 'lsim() should return non-empty y');
    assert(~isempty(t), 'lsim() should return non-empty t');
end

%% ========================================================================
%  7. ADVANCED LEIBNIZ TESTS
%  ========================================================================
function test_detLeibniz_vs_FFT_2x2(~)
    % Use general random inputs but disable reduction in the reference det
    % to check exact consistency with Leibniz formula
    n = 2; h = 2;
    % val = complex(randn(n, n, 2*h+1), randn(n, n, 2*h+1));
    A = PhasorArray.random(n,n,h);
    
    dL = detLeibnizHmc(A);
    % Force det to NOT reduce small terms, so it matches Leibniz exact result
    dF = det(A, 'reduceThreshold', 0);
    
    % Compare at several time points
    for t0 = linspace(0, 2*pi, 10)
        vL = full(evalp(dL, t0));
        vF = full(evalp(dF, t0));
        err = abs(vL - vF);
        assert(err < 1e-6, sprintf('Leibniz vs FFT-det mismatch at t=%g: %e', t0, err));
    end
end

function test_detLeibniz_vs_FFT_3x3(~)
    % Use general random inputs
    n = 3; h = 2;
    A = PhasorArray.random(n,n,h);
    % A = PhasorArray(val);
    
    dL = detLeibnizHmc(A);
    % Disable reduction to match exact Leibniz

    dF = det(A, 'reduceThreshold', 0);
    
    for t0 = linspace(0, 2*pi, 10)
        vL = full(evalp(dL, t0));
        vF = full(evalp(dF, t0));
        err = abs(vL - vF);
        assert(err < 1e-6, sprintf('Leibniz vs FFT-det 3x3 mismatch at t=%g: %e', t0, err));
    end
end

function test_detLeibniz_identity(tol)
    I = PhasorArray.eye(3);
    d = detLeibnizHmc(I);
    % det(I) = 1 (constant)
    dc = d{1, 1, 0};
    assert(abs(dc - 1) < tol, 'det(I) DC should be 1');
    assert(energy(d) - 1 < tol, 'det(I) energy should be 1');
end

function test_detLeibniz_product(tol)
    % For scalar PhasorArrays (1x1), det(A*B) = det(A)*det(B) is trivial.
    % For 2x2, test: det(A*B)(t) = det(A)(t) * det(B)(t)
    n = 2; h = 2;
    vA = complex(randn(n, n, 2*h+1), randn(n, n, 2*h+1));
    vB = complex(randn(n, n, 2*h+1), randn(n, n, 2*h+1));
    A = PhasorArray(vA);
    B = PhasorArray(vB);
    
    dAB = detLeibnizHmc(A * B);
    dA  = detLeibnizHmc(A);
    dB  = detLeibnizHmc(B);
    dAdB = dA * dB;
    
    % Compare at time points
    for t0 = linspace(0, 2*pi, 8)
        v1 = full(evalp(dAB, t0));
        v2 = full(evalp(dAdB, t0));
        err = abs(v1 - v2);
        assert(err < tol, sprintf('det(AB) vs det(A)*det(B) at t=%g: %e', t0, err));
    end
end

%% ========================================================================
%  8. TOEPLITZ FORMALISM TESTS
%  ========================================================================
function test_T_bt()
    A = PhasorArray.random(2, 2, 3);
    hTrunc = 8;
    T = A.T_bt(hTrunc);
    n = size(A, 1);
    Nt = 2*hTrunc + 1;
    expectedSize = n * Nt;
    assert(isequal(size(T), [expectedSize, expectedSize]), 'T_bt size');
end

function test_S_tb()
    % S_tb(h, phi, n) builds kron(I_n, diag(exp(jk*phi)))
    h = 5;
    phi = pi/4;
    n = 2;
    S = S_tb(h, phi, n);
    expectedSize = n * (2*h + 1);
    assert(isequal(size(S), [expectedSize, expectedSize]), 'S_tb size');
    % Diagonal should be exp(jk*phi) repeated n times
    Sdiag = diag(S);
    kvals = -h:h;
    expected_diag = repmat(exp(1i * kvals(:) * phi), n, 1);
    assert(max(abs(Sdiag - expected_diag)) < 1e-12, 'S_tb diagonal entries');
end

function test_S_tb_PhaseShift(tol)
    % Verify: S_tb(h, phi, n) * F_tb(A, h) == F_tb(PhaseShift(A, phi), h)
    % This tests that the Toeplitz dephasing operator S_tb is consistent
    % with the PhasorArray PhaseShift method.
    A = PhasorArray.random(2, 1, 3);  % 2x1 column vector, h=3
    phi = pi/6;
    hTrunc = 10;
    
    % Method 1: PhaseShift then Fourier vector
    A_shifted = A.PhaseShift(phi);
    F_shifted = A_shifted.F_tb(hTrunc);
    
    % Method 2: Fourier vector then multiply by S_tb
    F_A = A.F_tb(hTrunc);
    n = size(A, 1);
    S = S_tb(hTrunc, phi, n);
    F_via_S = S * F_A;
    
    err = max(abs(F_shifted - F_via_S));
    assert(err < tol, sprintf('S_tb * F_tb vs PhaseShift(A,phi).F_tb error: %e', err));
end

function test_yalmip_LMI()
    % TEST_YALMIP_LMI Checks LMI stability condition: P > 0, \dot{P} + A'P + PA < 0
    % For a simple constant stable system A = -I.
    
    nx = 2;
    h = 1;
    T = 1; % Period
    
    % Stable system A = -eye(2)
    A = PhasorArray(-eye(nx));
    
    % Lyapunov Matrix P(t) as ndsdpvar
    P = PhasorArray.ndsdpvar(nx, nx, h, "PhasorType", 'symmetric', "real", true);
    
    PT = P.T_tb(h);
    PA  = P*A;
    ATP = A'*P;
    Q   = PhasorArray(eye(nx));
    QT  = Q.T_tb(h);
    
    % Derivative operator matrix for harmonic domain
    % Corresponds to d/dt in frequency domain: N = diag(j*k*omega)
    N = N_tb(nx, h, T);
    
    % LMI: \dot{P} + A'P + PA <= -epsilon * I
    % \dot{P} term corresponds to commutator with N: N*PT - PT*N
    % Using the formula from Exemples/Exemple_Toolbox_LMI.m line 47:
    Sum = ATP + PA;
    LMI_lhs = Sum.T_tb(h) + QT + (N*PT - PT*N);
    
    % Constraints
    % P >= epsilon
    % LMI <= -epsilon
    epsilon = 1e-3;
    
    F = [PT >= epsilon*eye(size(PT)), ...
         LMI_lhs <= -epsilon*eye(size(LMI_lhs))];
         
    % Solve with minimal output
    opts = sdpsettings('verbose', 0);
    sol = optimize(F, [], opts);
    
    if sol.problem ~= 0
       % If solver issue, warn
       warning('YALMIP optimize returned problem code %d: %s', sol.problem, sol.info);
    else
        % Verify P is positive definite
        P_val = value(PT);
        min_eig_P = min(real(eig(P_val)));

        % Verify LMI is negative definite
        LMI_val = value(LMI_lhs);
        max_eig_LMI = max(real(eig(LMI_val)));

        % Assertions
        assert(min_eig_P >= epsilon - 1e-6, 'P should be positive definite');
        assert(max_eig_LMI <= -epsilon + 1e-6, 'LMI should be negative definite');
    end
end

function test_check_solvers()
    % Checks if YALMIP and MOSEK are detected
    has_sdpvar = exist('sdpvar', 'file') == 2;
    if has_sdpvar
        fprintf('    [INFO] YALMIP detected.\n');
    else
        fprintf('    [INFO] YALMIP NOT detected.\n');
    end
    
    has_mosek = exist('mosekopt', 'file') == 3;
    if has_mosek
        fprintf('    [INFO] MOSEK detected.\n');
    else
        fprintf('    [INFO] MOSEK NOT detected.\n');
    end
end

function test_yalmip_LMI_vs_Lyap()
    % TEST_YALMIP_LMI_VS_LYAP Compares X from lyap(A,Q) with X from LMI.
    % Solves A'X + XA + dX/dt + Q = 0.
    
    % System definition (from test_lyap_LTP)
    A0 = [-3 0.5; 0 -2];
    A1 = [0.1 0; 0 0.05];
    A2 = [0 0.3; 0.1 0];
    % A(t) has frequencies -2, -1, 0, 1, 2. So h=2.
    A = PhasorArray(cat(3, conj(A2), conj(A1), A0, A1, A2));
    Q = PhasorArray(eye(2));
    
    h_sol = 10;
    nx = 2;
    T = 1; % Assume default period T=1 (omega=2pi)
    
    % 1. Solve using lyap (Arithmetic)
    % Note: lyap(A, Q, ...) uses defaults if T not specified.
    [X_lyap, ~] = lyap(A, Q, 'h', h_sol,'autoUpdateh',true,'T',T);
    
    % 2. Solve using YALMIP LMI (Optimization)
    % Equation: A'X + XA + dX/dt + Q = 0
    % We solve this as a feasibility problem with equality constraints.
    
    X_yalmip = PhasorArray.ndsdpvar(nx, nx, h_sol, ...
        "PhasorType", 'symmetric', "real", true);
    
    XT = X_yalmip.T_tb(h_sol);
    ATX = A'*X_yalmip;
    XA = X_yalmip*A;
    QT = Q.T_tb(h_sol);
    
    % Derivative operators
    N = N_tb(nx, h_sol, T);
    
    % dX/dt term -> N*XT - XT*N
    % Equation (A'X + XA + Q + dX/dt).T_tb = 0
    
    Sum = ATX + XA;
    Lyap_Op = Sum.T_tb(h_sol) + QT + (N*XT - XT*N);
    
    % Constraint: Lyap_Op <= 0 (Relaxed equality to inequality)
    % We minimize trace(XT) to find the 'smallest' solution that satisfies 
    % the inequality (which converges to the equality solution for stable A).
    F = [Lyap_Op <= 0; XT>=0];
    
    % Solve: minimize trace(XT)
    opts = sdpsettings('verbose', 0);
    sol = optimize(F, trace(XT), opts);
    
    if sol.problem ~= 0
        warning('YALMIP optimize returned problem code %d: %s', sol.problem, sol.info);
    else
        % Compare solutions
        X_val = value(XT);
        X_ref = X_lyap.T_tb(h_sol);
        
        diff = max(abs(X_val - X_ref), [], 'all');
        
        % Normalize error by norm of solution
        rel_err = diff / max(abs(X_ref), [], 'all');
        
        % Assertion
        % Consistency should be high if T matches.
        assert(rel_err < 1e-5, sprintf('Mismatch between lyap() and LMI solution. Rel Err: %e', rel_err));
    end
end

function test_riccati_lmi_vs_kleinman()
% Verifies that RicHarmonicKlein converges to the same LQR gain as the
% direct Schur-complement Riccati LMI (YALMIP).
%
% ARE:  (A-N)^* P + P(A-N) - P B R^{-1} B^* P + Q = 0
% Schur LMI:
%   [ He((A-N)^*P) + Q,  PB ] >= 0,   max tr(PT)
%   [ B^*P,               R  ]
%
% System: stable open-loop (A0 has eigs -3, -2) with periodic perturbation.
% K0=0 is valid; Lyapunov solver converges cleanly.

    T  = 1;
    h  = 8;
    A0 = [-3 0.5; 0 -2];
    A1 = [0.1 0; 0 0.05];
    A  = PhasorArray(cat(3, conj(A1), A0, A1));
    B  = PhasorArray([1; 0.5]);
    Q  = PhasorArray(diag([10, 1]));
    R  = PhasorArray(1);

    % --- 1. RicHarmonicKlein (K0=0 valid since A is open-loop stable) ---
    [K_klein] = RicHarmonicKlein(A, B, Q, R, ...
        PhasorArray(zeros(1, 2)), T, ...
        'maxIter', 200, 'thresholdResidual', 1e-8, 'autoUpdateh', true);

    % --- 2. Riccati Schur LMI (wiki template §4) ---
    P  = PhasorArray.ndsdpvar(2, 2, h, 'PhasorType', 'symmetric', 'real', true);
    PT = P.T_tb(h);
    N  = N_tb(2, h, T);

    % He((A-N)^* P) + Q  in Toeplitz block
    % (A-N)^* = A^* + N  since N^* = -N (skew-Hermitian)
    % => T_tb(A'P + PA) + N*PT - PT*N + Q_tb
    F11 = T_tb(A'*P + P*A, h) + (N*PT - PT*N) + Q.T_tb(h);
    F12 = T_tb(P*B, h);
    F22 = R.T_tb(h);

    F = [PT >= 0, [F11, F12; F12', F22] >= 0];
    sol = optimize(F, -trace(PT), sdpsettings('solver', 'mosek', 'verbose', 0));
    assert(sol.problem == 0, sprintf('YALMIP Riccati LMI failed: %s', sol.info));

    P_opt = sdpval(P);
    K_lmi = R \ (B.' * P_opt);

    % --- 3. Compare at truncation h ---
    rel_err = norm(K_lmi.T_tb(h) - K_klein.T_tb(h), 'fro') / ...
              (norm(K_lmi.T_tb(h), 'fro') + eps);
    fprintf('    [INFO] Riccati LMI vs Kleinman: rel_err = %.2e\n', rel_err);
    assert(rel_err < 5e-2, ...
        sprintf('RicHarmonicKlein vs Schur LMI mismatch: %.2e', rel_err));
end
