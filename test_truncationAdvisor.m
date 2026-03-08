function ok = test_truncationAdvisor()
%TEST_TRUNCATIONADVISOR Smoke tests for truncationAdvisor utility.

fprintf('\n========================================\n');
fprintf('  truncationAdvisor - Test\n');
fprintf('========================================\n');

ok = true;

% Test 1: default call
adv = truncationAdvisor("lyap", "verbose", false);
assert(isfield(adv, 'recommended'));
assert(adv.recommended.h0 > 0);
assert(adv.recommended.hmax >= adv.recommended.h0);

% Test 2: LMI split orders
advLmi = truncationAdvisor("lmi", "profile", "balanced", "verbose", false);
assert(~isnan(advLmi.recommended.hLMI));
assert(advLmi.recommended.hLMI >= advLmi.recommended.h0);
assert(advLmi.recommended.hP > 0);
assert(advLmi.recommended.hA > 0);

% Test 3: override support
advCustom = truncationAdvisor("riccati", "h0", 12, "hmax", 80, "targetResidual", 1e-8, "verbose", false);
assert(advCustom.recommended.h0 == 12);
assert(advCustom.recommended.hmax == 80);
assert(abs(advCustom.recommended.targetResidual - 1e-8) < 1e-14);

fprintf('All truncationAdvisor tests passed.\n');
end
