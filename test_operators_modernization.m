% Test script for modernized mldivide and mrdivide
A = PhasorArray.random(3,3,5) + 5*eye(3);
B = PhasorArray.random(3,1,5);

fprintf('\n--- Testing mldivide (\\) ---\n');
try
    X1 = A \ B;
    res1 = norm(A*X1 - B);
    fprintf('mldivide residual: %e\n', res1);
    if res1 < 1e-3
        disp('SUCCESS: mldivide converged to expected accuracy.');
    else
        disp('WARNING: mldivide residual higher than expected.');
    end
catch e
    fprintf('ERROR in mldivide: %s\n', e.message);
end

fprintf('\n--- Testing mrdivide (/) ---\n');
try
    % For mrdivide, we test X * A = B  => X = B / A
    % Let's use B as a row vector for variety
    Brow = PhasorArray.random(1,3,5);
    X2 = Brow / A;
    res2 = norm(X2 * A - Brow);
    fprintf('mrdivide residual: %e\n', res2);
    if res2 < 1e-3
        disp('SUCCESS: mrdivide converged to expected accuracy.');
    else
        disp('WARNING: mrdivide residual higher than expected.');
    end
catch e
    fprintf('ERROR in mrdivide: %s\n', e.message);
end

fprintf('\n--- Testing with Name-Value arguments ---\n');
try
    % Verify that we can still pass options like threshold
    X3 = A.mlHmcDivide(B, 'thresholdResidual', 1e-8, 'autoUpdateh', true);
    res3 = norm(A*X3 - B);
    fprintf('Custom threshold (1e-8) residual: %e\n', res3);
catch e
    fprintf('ERROR in custom options: %s\n', e.message);
end
