# Changelog

Notable changes to the PhasorArray Toolbox. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/); versions follow `CITATION.cff`.

Entries state what changed and who is affected; where a design choice is not obvious,
the reasoning is kept inline rather than in a separate decision record.

## [Unreleased]

### Breaking

- **Random generators and `ndsdpvar` take a single `symmetry` argument.** It replaces
  `time_structure` on `rand_phasor` / `PhasorArray.random` and the `PhasorType`/`real`
  pair on `PhasorArray.ndsdpvar`.

  ```matlab
  % Before
  A = PhasorArray.random(3,3,5,"time_structure","symetric");
  P = PhasorArray.ndsdpvar(4,4,5,'PhasorType','symmetric','real',true);

  % Now  ("real" is the default, so the second line needs nothing)
  A = PhasorArray.random(3,3,5, symmetry = ["real" "symmetric"]);
  P = PhasorArray.ndsdpvar(4,4,5);
  ```

  `symmetry` accepts any combination of 14 names, listed in `phasorSymmetry`. Transposition,
  conjugation and time reversal are three commuting involutions, so they generate
  $(\mathbb{Z}/2)^3$: eight elements, each giving a `+` and a `-` property. The vocabulary is
  therefore exhaustive, and the 14 names combine into 51 classes against 8 reachable before.

  **Requesting two properties usually yields three.** `real` and `symmetric` force
  `hermitian`; the satisfied set is a subgroup with a character, so it always has 0, 1, 3 or 7
  entries, never 2, 4, 5 or 6. Ask for the closure with `symmetryClosure`, or take it as the
  second output of `phasorSymmetry` / `ndsdpvar`. Contradictory requests now error
  (`PhasorArray:symmetry:contradiction`) instead of silently returning zeros — the check is a
  sign collision on eight integers, so it needs no tolerance.

  **`retroHermitian` is now `paraHermitian`**, and `PhasorType='hermitian'` was also
  para-Hermitian despite its name: both mean $A(-t)' = A(t)$, i.e. every coefficient
  hermitian. Time-domain hermitian, $A(t)' = A(t)$, is the separate name `hermitian`
  and was not reachable before.

  **Who is affected:** `time_structure`, `hurwitzeig` and `Q` are **gone**, along with
  `PhasorType`/`real` on `ndsdpvar`. All examples, templates and docs are migrated.

  The three spectral routes constrained the spectrum, not the symmetry group, so they never
  belonged on the same axis. Each now has a dedicated generator, and both are more accurate
  than what they replace:

  ```matlab
  % Before                                          % Now
  rand_phasor(2,2,7,"hurwitzeig",[-1 -0.3], ...)     PhasorArray.randomWithNPole(diag([-1 -0.3]), 7)
  PhasorArray.random(n,n,h,"time_structure","sdp")   PhasorArray.randomSPD(n, h)
  ```

  `randomWithNPole` returns the prescribed exponents exactly (`[-0.3 -1]` measured against
  the same request) and is already real in time, so the `mreal(...)` wrapper the examples
  carried is now a verified no-op and was dropped. `randomSPD` produces an exactly symmetric
  positive-definite array 15x faster than the rejection sampling it replaces. `rand_phasor`
  goes from 128 lines to 73 and has a single branch.

- **Comparison operators now act on `A(t)`, not on the coefficient array.** `<`, `>`, `<=`,
  `>=`, `==` and `~=` each delegated one line to `pvalue`, comparing Fourier coefficients.
  That was wrong in three separate ways, and the first was silent:

  - MATLAB's ordering on complex numbers tests **real parts only**. Measured: `z < z+100i`
    is `false`. Half of every coefficient was discarded without a warning, and the verdict
    bore no relation to any ordering — on two SPD arrays with
    `min_t lambda_min(P(t)-Q(t)) = -0.577`, so `P` is *not* above `Q` in the Loewner sense,
    the old `P>Q` still returned 22 of 28 coefficients "true".
  - `A == B` raised `MATLAB:sizeDimensionsMustMatch` whenever the harmonic orders differed,
    while `A + B` padded them correctly. The comparison ignored the class invariant.
  - `&`, `|` and `~` threw `Argument must be real` on any complex-coefficient array, which
    is nearly all of them: a real-in-time signal has conjugate-symmetric complex coefficients.

  Each operator now returns an `[n x m]` logical whose entry `(i,j)` judges the scalar
  periodic function `A_ij(t)` against `B_ij(t)` over a whole period — the matrix shape, not
  the coefficient array. `MINUS` supplies the harmonic padding and scalar broadcasting, so
  differing orders work. Verified identical to a direct 5000-point sweep.

  `==` and `~=` need no sampling: the map between a periodic signal and its coefficients is
  a bijection, so equality of functions is equality of the padded coefficients. The ordering
  operators sample, and refuse a complex difference
  (`PhasorArray:compare:complexOrder`) instead of quietly dropping its imaginary half.

  Called by name, the four ordering operators grade the answer instead of only passing or
  failing it:

  ```matlab
  [holds, frac, crenel, Cph] = lt(A, B);
  %  holds  [n x m]        true when it holds for every t
  %  frac   [n x m]        fraction of the period over which it holds
  %  crenel [n x m x Nt]   the indicator, on theta = (0:Nt-1)/Nt*2*pi
  %  Cph                   that indicator as a PhasorArray, order Nt/2-1
  ```

  `A < B` takes the first output only, as MATLAB requires. `holds` is `frac == 1`, and `frac`
  is both `mean(crenel,3)` and the DC coefficient of `Cph` — measured identical, and matching
  a 20001-point reference to 2.7e-03. Outputs beyond the first are computed only if asked for.

  **`Cph` is a square wave**, so its coefficients decay only as `1/k` — measured `c_k*k`
  constant at 0.25 — and 46% of its energy sits outside DC at every grid size from 16 to 1024.
  Refining the grid resolves more of an infinite tail instead of converging. The one
  coefficient that is not Gibbs is DC, which is `frac`.

  `and`, `or` and `not` were removed: no time-domain reading of them exists, and MATLAB's own
  "undefined operator" message is clearer than the failure they produced. Nothing in the
  toolbox, the examples or the suites called any of these.

- **`ctrretro` returned `A(t)'` instead of `A(-t)'`**, duplicating `mctranspose`
  character for character. Fixed by dropping a `flip`: conjugating the pages already
  reverses the harmonic index. It had no callers.

- **`expm` now returns a `PhasorArray`**, as its help always claimed. It returned a raw
  double array, so the result could not be fed back into any `PhasorArray` operation.
  The phasor values are unchanged (6.7e-16 against a page-wise `expm` on the time grid).

  Its `plot`, `autoTrunc`, `reduceThreshold` and `reduceMethod` arguments were declared
  and never read; they now work. `autoTrunc` matters in practice — the time grid is far
  finer than the bandwidth of `exp(A)`, so the raw result carries 127 harmonics whatever
  the input order, against 11 once truncated. It stays `false` by default, so numerical
  output is untouched unless asked for.

- **Default period is now `T = 2*pi` everywhere.** The library previously mixed
  two conventions: some functions defaulted to `T = 1` (implying `ω = 2π`),
  others to `T = 2*pi` (implying `ω = 1`). The second is now the single rule: it makes the
  default pulsation `ω = 1`, which removes stray `2π` factors from the complex exponentials
  and matches the convention used in the publications.

  **Who is affected:** any call that relied on the default instead of passing
  `T` explicitly. Results are now computed over `[0, 2π]` instead of `[0, 1]`.

  ```matlab
  % Before: period 1 in some functions, 2*pi in others
  A = PhasorArray.random(2,2,3);
  At = PhasorArray2time(A);          % which period?

  % Now: 2*pi everywhere, and explicit is still best
  At = PhasorArray2time(A, 2*pi);
  ```

  Functions whose default changed include `N_bt`, `spN_tb`, `spN_bt`,
  `FloquetDec`, `TransitionMatrixOverT`, `rand_phasor`, and the `@PhasorArray`
  / `@PhasorSS` constructors and simulation entry points. **Passing `T`
  explicitly is unaffected**, and is the recommended habit.

- **`iszero` now raises an error on `sdpvar`/`ndsdpvar`-backed arrays**
  (`PhasorArray:iszero:decisionVariable`) instead of failing with an opaque
  `mustBeFloat` message. The predicate is undecidable before the optimisation
  problem is solved — YALMIP turns `x == 0` into a *constraint*, not a test.
  Apply it to the solved array instead.

- **`sparray2TBlocks` rejects `sym`/`sdpvar` coefficients explicitly**
  (`PhasorArray:sparray2TBlocks:specialType`). It never supported them; it used
  to fail inside `cell2mat` with an unrelated message. Use `array2TBlocks`.

### Added

- **`phasorSymmetry` and `symmetryClosure`.** `phasorSymmetry(A, req)` projects a
  `PhasorArray` onto a symmetry class; `symmetryClosure(req)` returns the properties that
  request actually implies, without building anything. The projection is a single group
  average $\frac{1}{|H|}\sum_{g\in H}\varepsilon(g)\,g\!\cdot\!A$, which is idempotent by
  construction, so no alternating iteration is needed.

  For `ndsdpvar`, the request is mapped onto a native YALMIP declaration whenever it fits,
  so no decision variable is created only to be projected away: `real`+`symmetric` on a
  3×3 with `h=2` declares 30 free variables, the same as the previous API, against 90 for
  an unconstrained complex variable.

- **The R2021b floor the README advertises is now actually held.** `tensorprod` is R2022a, and
  the toolbox met it with three different strategies: `PhasorArrayTimes` gated on the release
  and kept a matrix fallback, `PhasorArray2time` forced itself onto a slower sparse method, and
  `@PhasorSS` had no guard at all — one of its ten calls sat inside a `try` whose `catch` set
  the result to **zero**, so an old release would have got a silent wrong answer.

  Twelve of the thirteen calls turn out to be the *same* contraction, `(M, E, 3, 1)`: combining
  the harmonic dimension against a sampled basis. They now share one kernel, `harmonicCombine`,
  which needs no release gate at all — a mode-3 contraction is a matrix product on the
  unfolding, bit-identical to `tensorprod` and within 4% of it. Old releases therefore gain the
  *fast* path rather than a compatibility one, and `PhasorArray2time`'s forced downgrade is
  gone: its three methods now agree to 2.2e-16.

  The two genuine exceptions keep their own treatment. `PhasorArrayTimes` performs a double
  contraction `([2,3],[1,3])` — the Cauchy product — which is why `PhasorArrayTimes2` exists as
  a separate, genuinely slower fallback; its gate stays. `PhasorArray2time`'s second method
  contracted two 2-D operands, which is a plain product and is now written as one.

  `eig` gated on `verLessThan('matlab','9.12')` for `pageeig`, naming R2022a. Release numbers
  are the wrong instrument — being off by one calls a function that is not there — so it now
  asks `which('pageeig')`. Both branches verified to agree exactly.

  `contract3`, added earlier in this changelog's cycle, is this kernel under a name that leaked
  its dimension index. Renamed.

- **`PhasorArray.logm` removed; it was the editor's blank template.** Present since the initial
  commit, `function [outputArg1,outputArg2] = logm(inputArg1,inputArg2)` with a body that
  returns its inputs and a help block still reading "Summary of this function goes here". It
  overrode MATLAB's `logm` for the class, so `logm(A)` raised `Not enough input arguments`
  pointing inside the stub, and `logm(A,B)` returned `A` unchanged — a wrong answer, silently.
  The classdef also declared it as `out = logm(A,nvp)`, one output where the file defines two.

  `FloquetDec` is unaffected: its `logm(Q)` takes the eigenvalue matrix of `Phi(T)`, a plain
  double, so it always reached the builtin. Verified after removal — exponents `-0.3, -1` as
  prescribed.

  Removing it beats implementing it on the spot: `logm(A)` now raises
  `MATLAB:logm:inputType`, which says what is wrong. A periodic matrix logarithm is a Floquet
  question, not a pointwise one, and deserves its own work.

- **The ACDC template solved the wrong Lyapunov equation.** `LyapHarmonic` and `SylvHarmonic`
  take a **pulsation** as their last positional argument; the template passed `T = 1/f`, its
  **period**. At 50 Hz that is 0.02 against 314, a factor of 15700. Its `P` therefore left a
  residual of **1.4e-03** where the correct pulsation gives **4.2e-17**.

  The audit had passed it throughout, because it checks that `P` is positive definite — which
  it is, `min eig 2.31e-05` either way — not that it solves anything. A certificate that
  satisfies no equation certifies nothing.

  Migrating to `lyap`, which takes the period by name, fixes it by construction: residual
  `1.1e-15`.

- **The tutorials call `lyap` rather than the kernels beneath it.** `BasicToolbox` and
  `GettingStarted` reached for `SylvHarmonic` directly, twice each, in the shape
  `SylvHarmonic(A', A, Q, ...)` — which *is* the Lyapunov equation, written as a Sylvester.
  `lyap(A, Q, ...)` says so, and gives the reader the tool they will actually reuse. Verified
  identical on all three forms, `0` to `1e-15`.

- **`functionSignatures.json` rebuilt; it had never worked.** Present since the initial commit
  and never touched, the editor completion file was inert in four independent ways:

  - it described `Lyap_Harmonique`, a name renamed to `LyapHarmonic` long ago — the same lost
    rename as `Sylv_harmonique` in `place`;
  - that one key appeared **three times** in a single JSON object, so two of the three entries
    were discarded at parse whatever else was true;
  - the three entries described three *different* solvers — Lyapunov, Sylvester, Riccati —
    under that one name;
  - and it sat in `Fonctions/`, which holds only the two test scripts. MATLAB applies the file
    to functions in its own folder, and the solvers live in `pArrayBasicOperations/`.

  Its Riccati description also carried a sign error, `+ X B R^-1 B' X` where the equation has
  `- S B R^-1 B' S`.

  Rewritten as three correctly named entries in `Fonctions/pArrayBasicOperations/`, covering
  the positional arguments and the name-value options a caller actually reaches for.

- **The examples call `hare`, the documented entry point, instead of the engine underneath.**
  `hare` mirrors `care`/`icare`, unifies the plain and descriptor equations behind one `'E'`
  argument, and defaults `autoUpdateh` to true — yet every example, template and doc snippet
  reached past it for `RicHarmonicKlein` directly. It had no callers at all.

  All seven sites now use it. `hare` gained `maxh` and `skipValidate`, without which two of
  them could not migrate; deeper tuning (`updateMethod`, `stagnationWindow`,
  `warmStartFraction`) stays on `RicHarmonicKlein` and is documented as such. On identical
  options the two agree bit for bit.

  Two call sites pinned `warmStartFraction` — 0.7 in `periodic_lqr_riccati`, 0.95 in
  `ECC_ex` — so they would have kept the slow path whatever the default became. Removed. The
  LaTeX documentation offered `"max_iter"`, `"residualThreshold"`, `"htrunc"` and `"hmax"`,
  none of which are real option names, so the published snippet could only error. Corrected.

- **`MathieuEquation_Tutorial` §15 no longer builds `K0` by lifting.** The trick assembled
  `A_hss = T(A) - N` at `h = 7` and called `lqr` on the 75x75 result. It is ill-posed: the
  Mathieu plant is undamped, so every lifted exponent `mu + jk*omega` lies on the imaginary
  axis, and the modulator states are pure integrators reachable only through `x`. Six modes
  come out neither stable nor controllable — smallest PBH singular value **4.6e-16** — and
  `lqr` refuses the problem, correctly.

  The same matrices, saved once and loaded into both releases, settle where the change lies:
  R2024b Update 6 returns a finite gain of norm 10.45, R2026a Update 4 raises
  `Control:design:lqr1`. The lenient one was the old release, and the example depended on it.

  `K0 = []` hands the initial gain to the solver's own fallback cascade, which stabilises the
  system in **0.3 s** with a closed-loop exponent of `-0.506`. The passage explaining the
  lifting trick is gone with it, replaced by a note on why it cannot work here.

- **`RicHarmonicKlein` warm-starts at the previous order instead of 70% of it**
  (`warmStartFraction` default 0.7 to 1.0). Stepping back looks for a smaller `h` as the gain
  converges, and it does find one — but the search costs far more than the saving. Measured
  across five problems, `f = 1.0` wins every time by 1.5x to 4x:

  | problem | f = 0.70 | f = 0.95 | f = 1.00 |
  |---|---|---|---|
  | stable h=6 | 3.6 s (h=25) | 2.4 s (h=25) | **1.5 s** (h=26) |
  | stable h=10 | 5.1 s (h=30) | 2.7 s (h=30) | **1.4 s** (h=30) |
  | unstable h=6 | 7.1 s (h=34) | 4.2 s (h=34) | **2.5 s** (h=35) |
  | unstable h=10 | 27.4 s (h=40) | 9.9 s (h=40) | **6.9 s** (h=46) |
  | unstable h=10, large | 125.0 s (h=62) | 41.6 s (h=71) | **27.7 s** (h=88) |

  Residuals stay at or below 2.7e-08 throughout.

  The design error is that the backstep optimises the **size of the returned object** during
  the search, when that size is adjustable afterwards for nothing. Truncating `K` from `h = 94`
  down to `h = 40` costs **0.4 ms**, leaves the closed-loop exponent at `-0.72071` unchanged,
  and discards `1.1e-06` of its energy — `K` simply does not live above 40. Solving fast and
  truncating dominates searching slowly for the same order.

  An extrapolated warm start was built and measured before being dropped: reading the decay
  slope to predict the minimal order estimates the need of the iteration just finished, not of
  the next one, whose operands have drifted up. It tied with `f = 1.0` at best across the same
  five problems, for 35 lines and two options.

- **`adaptiveHSolve` no longer crawls when the operator is wide.** Its `1.1*hOp` gate holds
  the stepper at `+1` until the solution order passes the operator's spectral width, on the
  reasoning that the residual is not yet decaying regularly. That assumes convergence lands
  well above `hOp` — which fails on a Kleinman iterate, where `A - B*K` accumulates harmonics
  faster than the solution needs them. Traced at `verbose=3`: `hOp = 89`, gate at 97.9,
  convergence at 88, so the gate never opened and the loop took **34 unit steps of 2–5 s each,
  114.7 s for one outer iteration**, every row labelled `initial`.

  The gate now also opens after `cfg.maxUnitSteps` unit steps, default 5. That bounds the
  waste whatever `hOp` does, and is cheap to get wrong: the extrapolated jump was already
  damped to 80% and clamped to `h/2`, so a premature slope costs one solve, not thirty-four.
  A second lock needed the same treatment — the slope window starts at the first sample past
  `1.1*hOp`, which does not exist in this case, and now falls back to the last `maxUnitSteps`
  orders.

  On the Riccati that exposed it: **324.9 s to 120.8 s, 2.7x**, with `h_history` unchanged at
  `[38 79 88 84 70 62 62]`, the same 7 iterations, residual `2.411e-07` and closed loop
  `-0.7207` to the digit. Same trajectory, same answer. On a case where the gate already
  opened (`hOp = 45`) the order path is identical step for step and the solution matches to
  `0.00e+00`.

  `maxUnitSteps` deliberately avoids the word stagnation: `stagnationWindow` watches the
  **residual** and aborts, this watches the **stepper** and pushes on.

- **`adaptiveHSolve`**, the shared harmonic-order refinement engine. Each solver passes a
  closure `solveAtH(h)` and a configuration struct; the engine owns the loop. It estimates
  two convergence slopes from the residual history,

  ```
  s_exp = [log(e2) - log(e1)] / (h2 - h1)          exponential regime
  s_alg = [log(e2) - log(e1)] / [log(h2) - log(h1)]  algebraic regime
  ```

  and declares an algebraic regime for `-1.5 < s_alg < -0.1`, otherwise an exponential one
  for `s_exp < -1e-4`. It can then extrapolate whether the residual target is reachable
  before `maxh` and stop early instead of letting precision collapse silently. `lyap`,
  `lyapG` and `mlHmcDivide` now share it, and publish the same `info` diagnostics.

- **`randomWithNPole` takes a period.** Its rotation methods parameterise the angle
  `theta` by the normalised angle, so the time derivative carries `omega = 2*pi/T`.
  That factor was missing, and the prescribed exponents drifted as soon as the caller
  evaluated at another period: `-1 -> -1.0726` at `T = 0.1`. With `T` supplied the
  error is flat, `6.3e-05` for `'truncated'` and `1.6e-13` for `'generic'`.

  `'structured'`, the default, needed no correction and takes `T` only for symmetry of
  the interface: its mismatch term is strictly lower triangular, so it leaves the
  diagonal of `J` — hence the exponents — untouched at any period. Verified exact
  (`4e-13`) at `T = 2*pi`, `1` and `0.1`.

  `randomSPD` deliberately has **no** `T`: no derivative enters `P = L L'`, so
  positivity is a property of the coefficients and holds at any period (measured
  identical at the three). Documented rather than given a parameter that does nothing.

  Default `T = 2*pi` reproduces the previous output bit for bit.

- **`hmq_sim` accepts `solver = "ode78"` and `"ode89"`.** It hard-coded `ode15s`, a stiff
  solver of order at most 5. Periodicity is oscillation, not stiffness: an oscillatory
  right-hand side forces small steps for *accuracy*, which implicit machinery does not help
  with, while the low order costs dearly. Measured on a 4x4 with prescribed Floquet exponents,
  `RelTol = AbsTol = 1e-12`:

  | solver | exponent error | time | steps |
  |---|---|---|---|
  | ode15s | 1.29e-08 | 2.40 s | 4072 |
  | ode78  | 6.84e-12 | 1.44 s |  281 |

  So ode78 is both more accurate and faster here — it is not a trade-off. `'adaptative'`
  (ode15s) remains the default, because genuinely stiff systems, such as the switching
  converter template, still want it. Both solvers are R2021b, the documented floor.

- **`TransitionMatrixOverT` warns when the monodromy matrix cannot carry the answer.**
  `Phi(T)` holds the multipliers `exp(mu_i*T)`; once the smallest sits below the roundoff of
  the largest, the fast exponents are unrecoverable — at any tolerance, with any solver.
  Measured on a 2x2 with prescribed exponents, ode15s and ode78 agree in being wrong:

  | span of \|E\| | exponent error |
  |---|---|
  | 1.9e-03 | 1.5e-08 |
  | 3.5e-06 | 8.3e-06 |
  | 1.1e-11 | 1.5e-02 |

  The function used to return those numbers silently. It now raises
  `TransitionMatrixOverT:multiplierUnderflow` and points at `HmqNEig`, which diagonalises
  `T_h(A)-N` and never exponentiates: exact to **1.4e-13** across the same range, including
  `mu = -200` where the monodromy route is off by 196, and 20-80x faster. The warning quotes
  the span only — the internal estimate that triggers it tracks the magnitude too loosely to
  be published as an accuracy.

- **`setHarmonic`, and a warning when a brace assignment breaks realness.** Writing one
  harmonic left its opposite untouched, so `A{1,1,2} = 0.3+0.1i` silently cost `A(t)` its
  realness — measured symmetry defect `3.3e-01`, `isreal` dropping to false, no diagnostic.
  The invariant was the user's to remember, at every edit.

  Brace assignment now raises `PhasorArray:braceAssign:brokenConjugatePair` when the pair it
  writes stops being conjugate. The test is **local**, on the written pair only: a global
  `isreal` is `O(n*m*h)` and would have added 60% to 205% to the cost of an assignment,
  measured across `h = 8` to `120`.

  `setHarmonic(A, IDX, VAL)` writes both halves in one call, `IDX` taking a bare `k`, a triple
  `[i j k]`, or a cell `{i, j, k}` whose entries may each be `':'`, a scalar or a vector:

  ```matlab
  A = setHarmonic(A, 2, 0.3+0.1i);           % harmonic +-2, A(t) stays real
  A = setHarmonic(A, {1,':',3}, [1 2]);      % row 1 of harmonic +-3
  A = setHarmonic(A, 3, 1i, isReal=false);   % +3 only, deliberately complex
  ```

  `isReal` defaults to `isreal(A)`, so an array that is real stays real and one that is not is
  left alone. Degrading realness on purpose remains possible and explicit — nothing is
  repaired behind the caller's back. Harmonics beyond the current order are padded in.

  Two ambiguities are resolved rather than guessed: a three-element numeric index is always
  `[i j k]`, so several harmonics need the cell form; and naming both `+k` and `-k` while
  mirroring errors out (`setHarmonic:mirroredTwice`) instead of letting the second write
  silently overwrite the first.

- **`regularize`**, the mollification of the erratum to *"A TBLMI Framework for Harmonic
  Robust Control"* (Vernerey, Riedinger, Daafouz). Banded truncation converges in operator
  norm only when the Fourier series converges **uniformly**, which `A` in `L^inf` alone does
  not give. Convolving with the `C^inf` bump `exp(-1/(1-t^2))` produces a `C^inf` `A_eps` for
  which it does, and the trajectories converge back uniformly as `eps -> 0`.

  Convolution in time is multiplication of the coefficients, so it costs one windowing:
  `A_k` becomes `phihat(eps*k*omega)*A_k`. `phi` is real and even, hence so is the window —
  no phase is introduced, and `eps = 0` returns the input unchanged.

  ```matlab
  Aeps = regularize(A, 0.05);        % T defaults to 2*pi
  ```

  Convergence is second order in `eps`, as a symmetric mollifier gives: measured
  `1.3e-01, 1.7e-02, 2.0e-03, 1.8e-04, 2.0e-05` for `eps = 1, 0.3, 0.1, 0.03, 0.01`.

  **It does suppress Gibbs**, which is what makes it useful beyond the theory. On the square
  wave returned by the comparison operators, the overshoot above the plateau falls from
  **14.1% to 0.2%** at `eps = 0.1`, while the plateaus themselves are untouched (median error
  away from the transitions reaching 0.0000).

  `phihat` decreases from 1 to 0 over `[0, 5]` and first vanishes at 4.9965, its envelope then
  falling sub-exponentially — 4.6e-02 near 10, 6.8e-04 near 40. Keeping `eps*h*omega` inside
  that first lobe leaves the window positive and monotone; past it the window has zeros, which
  is faithful mollification but no longer reads as a low-pass.

  **It is deliberately opt-in, not applied after every transform.** On band-limited data it
  removes nothing and costs `O(eps²)` of bias — measured 2.0e-05 at `eps = 0.01` rising to
  1.7e-02 at 0.3. And it does not undo aliasing, which is the artefact an evaluate-then-
  transform round trip actually risks: content at `k = ±35` sampled at `Nt = 32` folds onto
  `k = 3` with amplitude 0.5, and regularization leaves it at 0.4859 even at `eps = 0.2`,
  because a fold deposits energy at low `k` where the window is 1 by design. Oversampling is
  the defence there, and `abs` and the comparison operators already apply it.

  The erratum takes the same position: practical implementations may use the truncated matrix
  directly, the regularization being what guarantees convergence in the infinite-dimensional
  setting.

- **`abs` on a `PhasorArray`**, absent until now, returning `|A(t)|` taken pointwise in time:
  sample over one period, take the modulus, transform back. Complex `A(t)` is fine, the
  modulus being real either way. Checked against the closed form for `|cos|` — DC `0.636588`
  against `2/pi = 0.636620`, second harmonic `0.212239` against `2/(3*pi) = 0.212207`.

  `|A(t)|` has a corner wherever `A` crosses zero, so it is not band-limited and the result is
  a truncation — but a well-behaved one. Unlike the comparison indicator, it is continuous, so
  it does not ring; the coefficients decay as `1/k²` and the retained energy is 99.93% already
  at `h = A.h`. The pointwise error nevertheless converges only as `1/h`: measured 41%
  relative at `h = A.h`, 8% at `4*A.h`, 1.2% at the full sampling order.

  `ABS(A)` therefore returns everything the sampling supports rather than silently truncating
  back to `A.h`, which would be a 41% approximation. `ABS(A, h=H)` truncates on request.

- **`randomSPD` takes `alphaMin`**, a prescribed lower bound: `P(t) ⪰ alphaMin·I` for all `t`.
  `P = L(t)L(t)'` is already positive semi-definite pointwise, so adding `alphaMin·I` gives
  the bound exactly, with no rejection sampling and no eigenvalue check. It is a pure DC term,
  so the bandwidth and the lossless truncation are unaffected. Default `0` reproduces the
  previous output bit for bit.

  The guarantee is deliberately **pointwise**, not on the truncated operator. `T_h(P)` is the
  compression of the infinite Toeplitz operator, so `lambda_min(T_h(P))` only decreases
  towards `min_t lambda_min(P(t))` as `h` grows — measured `0.0678, 0.0629, 0.0615, 0.0611`
  at `h = 6, 12, 24, 48` against a pointwise `0.060995`, approaching from above and never
  crossing. `T_h(P) ⪰ alpha·I` is therefore necessary but not sufficient; the equivalence with
  `P(t) ⪰ alpha·I` holds at `h = inf` alone.

- **`pageEnergy` profiling and a `matlab.unittest` suite** — see *Changed* and *Tests*.

### Fixed

- **`isimag` was wrong for every non-zero purely imaginary array.** It tested
  `A == mimag(A)`, but for a purely imaginary `A` the relation is
  `A == 1i*mimag(A)` — the factor `i` was missing. Since `A == mimag(A)` only
  holds for `A == 0`, the predicate was in practice equivalent to `iszero`. It
  now answers correctly for `A = i`, `A = i·cos(t)`, and so on. Nothing inside
  the toolbox called it, which is why it went unnoticed.
- `mreal`, `mimag` and `mconj` now work on `ndsdpvar`-backed PhasorArrays.
  Two YALMIP gaps on the N-D type were being hit: `mrdivide` is not overloaded,
  and `real`/`imag` work on a base variable but not on a derived expression.
  The formulas are unchanged (verified bit-identical on doubles).
- `iszero` works on symbolic arrays (it went through `pagenorm`, which requires
  floating-point input).
- **The R2021b floor is now actually held on the `@PhasorArray` path.** The class header
  promises R2021b, but three entry points reached for newer functions with no guard:

  - `iszero` called `pagenorm` (**R2022b**). Replaced by a local `pageInfNorm`
    reproducing `norm(.,inf)` page by page, including its switch to vector semantics on
    a single-row or single-column page — the subtlety that makes a naive `max(abs(...))`
    wrong as soon as the page is a matrix. Matches `pagenorm` to 4.4e-16 over 500 random
    shapes, exactly on the edge cases.
  - `hmq_sim` called `tensorprod` (**R2022a**) three times. It is the engine behind
    `PhasorArray.lsim`, `FloquetDec` and `TransitionMatrixOverT`, so simulation and the
    whole Floquet machinery were unusable on R2021b.
  - `freqVarying_hmqSim` called it twice.

  Both now use the new `contract3(M, E)`: a mode-3 contraction is a plain matrix product
  on the unfolding, so the replacement is bit-identical over 300 random shapes and within
  4% of `tensorprod` up to 40×40×301. `PhasorArray2time` and `PhasorArrayTimes` already
  had release guards and keep them.

  `@PhasorSS` still calls `tensorprod` unguarded, which is consistent: `toLPVss` and
  `toLTVss` need the Control System Toolbox R2022b, so the class is out of reach on
  R2021b whatever happens. `@sparsePhasorArray` is unmaintained.
- `isreal` / `isimag` / `iscomplex` no longer surface raw YALMIP internals on
  `ndsdpvar`. They return the documented "cannot determine" result and now say
  so with a warning instead of answering `false` silently.
- `TimeArray2plot` replaced: the previous implementation had no callers and its
  `explosed = false` branch referenced an undefined variable.
- Stale documentation corrected in `RicHarmonicKlein`: `relChangeThreshold` was
  described as `‖Sk-Sk-1‖/‖Sk‖` while the code computes the **K**-based ratio.
- Removed two leftover console dumps of Floquet exponents in the Riccati
  solvers' `K0` failure path.

### Performance

- **`isreal` and `isimag` no longer allocate a projection.** They used to build the
  `mreal`/`mimag` result (three full allocations) and then call `ismembertol` twice — a
  set-membership test that sorts its inputs. They now test the equivalent symmetry of the
  harmonic coefficients directly, one `flip` and one subtraction:

  ```
  A(t) real           <=>  A_{-k} =  conj(A_k)
  A(t) purely imag.   <=>  A_{-k} = -conj(A_k)
  ```

  The gain grows with the array; it is not worth quoting as a single figure:

  | array | before | after | speed-up |
  |---|---|---|---|
  | 3×3×7 | 0.00009 s | 0.00008 s | 1.1× |
  | 10×10×51 | 0.00134 s | 0.00008 s | 15.7× |
  | 20×20×151 | 0.01201 s | 0.00049 s | 24.5× |
  | 40×40×301 | 0.11098 s | 0.00390 s | 28.5× |

  Since `mtimes` calls `isreal` to decide whether to re-symmetrise the product, the saving
  compounds: the realness check drops from 73% to 16% of a matrix product, and the
  `tolReal` bisection (~60 calls) goes from 7.9 s to 0.29 s at the largest size.

  The tolerance stays **relative** (`1e-12 · max|A|`, matching `ismembertol`'s default).
  This matters: round-off breaks the symmetry bit-exactly as soon as a computation runs
  long, but the drift stays around `1e-16` *relative* to `max|A|` while the absolute drift
  grows with magnitude. An absolute test would wrongly reject long chains of genuinely real
  arithmetic.

  `sym`/`sdpvar` keep the previous projection route. `issymmetric` and `ishermitian` are
  unchanged — they are per-page tests and were never on this path.

### Changed

- **`array2BToepliz` is now spelled `array2BToeplitz`** — the second `t` was missing.
  The old spelling is gone; it is a kernel function with no call site outside
  `Fonctions/`, so no script is affected.

- **`N_tb`, `N_bt`, `spN_tb` and `spN_bt` call their harmonic order `h`, not `nh`**,
  matching `S_tb`, `S_bt`, `pr`, `PR_In` and `SylvHarmonic`. The argument is positional,
  so calls are unaffected — and the help already said "Scalar h" next to an argument
  named `nh`.

- **Help headers name their own function again**, which is what makes `help` and
  `lookfor` resolve them: `array2TFTB` announced `%F_tb`, `PhasorArrayOplus` announced
  `%OPLUS`, and both `PhasorArrayTimes` and `PhasorArrayTimes2` opened on `%Takes`.
  All 66 files in `pArrayBasicOperations` were checked.

- **`help` works again on 15 kernel functions.** MATLAB only reads the comment block
  that immediately follows the function signature, so an `arguments` block placed above
  the documentation makes `help` and `lookfor` return nothing. Adding `arguments` blocks
  had done exactly that to `BT2array`, `BT_2_TB`, `TB_2_BT`, `PR_In`, `spPR_In`, `pr`,
  `S_tb`, `S_bt`, `PhasorArrayPad`, `pagekron`, `dephase`, `isFunny`, `isToeplitz`,
  `nHarm`, `shuffle_matrix` and `BToeplitz` — `S_tb` and `S_bt` among them, which the
  examples call directly. The blocks were moved back above `arguments`; no code changed.

  `Phasor2CosSin` and `PhasorArrayKron`, plus `disph`, `epsHAPV`, `manageTiledLayout`
  and `plotAngularSFT`, have no documentation to restore and still return nothing.

- **`help T_tb` and `help T_bt` answer for the first time.** The sweep was extended from
  one function per file to all 487 functions in the tree, which reached the 174 methods
  inside the `classdef` bodies. `T_tb` and `T_bt` — the most-called functions in the
  toolbox, 105 call sites across examples, templates and docs — carried no help at all:
  they are thin aliases and the documentation lived on `TB` and `BT`, the original names.
  They now document themselves and point to the aliases, along with `spBT` and the
  accessors `dc`, `ac`, `value`, `Value` and `isspecial`. `fromF_tb` had a full help block
  sitting below its `arguments` block; it was lifted.

  What remains undocumented is the framework surface MATLAB documents itself —
  `parenReference`, `braceAssign`, `getPropertyGroups` and the operator overloads.

- **Descriptor and base Riccati solvers merged.** `RicHarmonicKlein` now takes
  the mass matrix as name-value `'E'`; `RicHarmonicKleinGen` is a thin
  compatibility shim over it. Existing calls to either name keep working, and
  results are bit-identical. `lyapG` handles the generalised equation
  $\frac{d}{dt}(E'PE) + \dots = 0$ without ever inverting `E(t)`, by assembling the
  augmented Toeplitz-Block operator directly. The `direction` (forward/backward) and
  `derivativeForm` (product/sandwich) options make the dual formulations explicit.
- `RicHarmonicKleinGen` gained the `K0` fallback cascade the base solver already
  had (DC LQR, then truncated-HSS LQR). It previously gave up with an error
  where the base engine recovers.
- `hare` keeps the specialised fast path when no mass matrix is supplied
  (~1.6× cheaper than solving with `E = I`).
- **The five energy methods share one kernel, `energyOf`.** `energy`, `realEnergy`,
  `imagEnergy`, `ACenergy` and `DCenergy` differ only in which part of the array they
  hand it. All rest on Parseval: the $L^2$ energy $\int|X(t)|^2\,dt$ is the Frobenius
  norm of the Fourier coefficients, $\sum_k|X_k|^2$. `pageEnergy` (unchanged) profiles
  that sum harmonic by harmonic.
- `mlHmcDivide`'s `info` struct gains a `time_history` field, so the three
  adaptive solvers now publish the same diagnostics.
- **YALMIP decision graphs are preserved through the low-level operators.** Divisions are
  written as fractional multiplications (`* (1/2) * (-1i)` rather than `/ (2i)`) and
  matrix assembly uses `cat` rather than `cell2mat`, because both would otherwise collapse
  an `sdpvar`/`ndsdpvar` to `double` and silently break LMI synthesis. Verified
  bit-identical on doubles.
- `installToolbox` no longer puts `scratch/` on the MATLAB path. A file there
  whose name collides with a toolbox or MATLAB function used to shadow it — and
  since the folder is untracked, on one machine only.
- `installToolbox` reports when `savepath` fails. It cannot write `pathdef.m`
  without elevation on a standard Windows install, and used to fail silently.

### Tests

- **Migrated to `matlab.unittest`.** The previous suites were monolithic `assert()`
  scripts that stopped at the first failure. 102 tests now run as
  `matlab.unittest.TestCase` classes, so failures are reported individually, and
  `run_all_tests.m` provides a single entry point for CI. The legacy
  `test_PhasorArray_basic` / `test_PhasorArray_advanced` scripts still pass (49/49 and
  34/34) and remain the reference until the migration is signed off.
- `PhasorArrayRegressionTest` pins the atomic operators `T_tb`, `T_bt`, `N_tb`, `N_bt`
  against hardcoded matrices, fixing the index semantics rather than only relative
  identities such as $T_{tb}\,T_{tb}^{-1} = I$. Note that these expectations were
  themselves wrong on arrival — `T_tb`/`T_bt` were swapped and the `j` was missing from
  `N_tb`/`N_bt`, which are $\mathrm{kron}(I_n,\ \mathrm{diag}(jk\omega))$ and
  $\mathrm{kron}(\mathrm{diag}(jk\omega),\ I_n)$ respectively. They now match the code.
- `PhasorArrayAdvancedTest` adds an LMI stress case: $A'P + PA + BK + K'B' + Q$ assembled
  with deliberately mismatched orders ($h_P=2$, $h_A=1$, $h_B=3$, $h_K=1$, $h_Q=0$) and a
  mix of `double` and `ndsdpvar`, so that padding and Cauchy products cannot quietly drop
  the decision graph to `double`.
- `PhasorArrayCoreApiTest` locks the class semantics: paren/brace indexing, asymmetric
  concatenation, the power-electronics transforms (Park, Concordia, dq0), base algebra,
  and Parseval consistency.

### Examples and templates

- **`Exemples/Exemple_Toolbox_LMI.m` produced an unusable gain.** The synthesis LMI
  reported success while returning `||K|| ~ 1e8` and an unstable closed loop, on every
  draw tested. Two causes, both in the formulation:

  - The objective minimised `trace(Q_0)`. The constraint
    `A Q + Q A' - dQ/dt + B Y + Y' B' <= -2*alpha*Q` is homogeneous of degree 1 in
    `(Q, Y)`, so nothing bounds the problem from below except the floor
    `T(Q) >= 1e-6*I`; the solver parks `Q` on it, `Q` turns singular and `K = Y/Q`
    diverges. The objective now maximises the trace, pushing `Q` away from that boundary.
  - There was no quadratic coupling between `Q` and `Y`. A bare decay-rate inequality
    does not bound the gain. The LQ Schur complement of
    `templates/Optimal-Control/periodic_lqr_schur.m` does, and is now used; the `alpha`
    shift moved into its (1,1) block.

  Also `T_h(B)*T_h(Y)` became `T_h(B*Y)`: truncation does not commute with the product.

  | | before | after |
  |---|---|---|
  | `min_t lambda_min(Q(t))` | -3.2e-07 | +6.6e-04 |
  | `norm(K, inf)` | 1.2e+08 | 1.0e+01 |
  | closed-loop `max Re(mu)` | +2.4e+06 | -4.41 (open loop -2.17) |

  The example now **certifies** the result instead of trusting the solver status: it
  rebuilds `R = A Q + Q A' + B Y + Y' B' - dQ/dt` with untruncated products and checks
  `Q(t) > 0`, `R(t) < 0`. Since `B Y = B K Q`, that is the Lyapunov inequality for
  `A + B*K`, so stability is proved rather than sampled. A feasible LMI on truncated
  operators proves nothing on its own.

- **`PhasorArray.place` always threw.** It called `Sylv_harmonique`, a name that no longer
  exists — the function was renamed `SylvHarmonic` and this call site never followed, so every
  invocation raised `MATLAB:UndefinedFunction`. The signatures are identical, five positional
  arguments in the same order, so the fix is the rename. With a square `B` the method now
  places the requested Floquet exponents to `7.6e-15` at `hLyap = 8`.

  Found while checking that `place` survived the `randomSPD` migration; nothing else calls it.

  A sweep of every identifier called across the 154 source files found **one** other name
  resolving nowhere: `validateInputs_v2`, called twice by
  `Fonctions/Display and data manipulation/angularsft_v2.m`. Unlike `Sylv_harmonique` this is
  not a lost rename — the function is absent from the repository, from the personal library,
  and from the archived `angularsft_v2_old.m`, which carries the same three subfunctions. It
  was never written, so `angularsft_v2` has never run, despite a 350-line README documenting
  it. Nothing calls it. Left alone: the angular/dSPACE cluster is frozen.

- **`TransitionMatrixOverT`'s `enforce_real` was documented backwards.** The help read
  "Force complex-valued ODE integration even if A(t) appears real", while the code sets it to
  `true` precisely *when* `isreal(A)` — it selects the real path, as its name says. The
  warning compounded it by telling the user to set `isRealValued`, an option this function
  does not expose; its name here is `enforce_real`, and it is forwarded to `hmq_sim` as
  `isRealValued`. Both corrected.

- **The spectral-decay figure never showed decay.** `rand(1,1,h+1) .* exp(-0.3*(0:h))`
  expands a row vector against dimension 2, giving a `1 x 16 x 16` array whose
  first column carries no decay -- and that is the column the script read back.
  Wrong since the figure was first generated. The regenerated figure falls over
  five decades. Its `sgtitle` also landed on top of the title `pageEnergy` draws
  itself, and is dropped.

- **Two `pageEnergy` calls in `docs/generateWikiFigures.m` were positional** where the
  signature is name-value, so they raised `MATLAB:TooManyInputs`. They had never been
  reached: the script dies earlier on `exportgraphics(..., '.svg')`, which MATLAB only
  accepts from R2025a and which raises `MATLAB:print:InvalidFileFormatForExport` before
  that. Run the figure generators on R2025a or newer.

- **`Exemples/MathieuEquation_Tutorial.m` gave the same exponent two opposite verdicts.**
  §2 tested `< 0` and §15 tested `< 1e-6`, so an exponent of `+5e-14` printed as
  `+0.000000 → UNSTABLE` in one section and `+0.0000 → STABLE` in the other. Both
  systems are genuinely on the imaginary axis — the undamped Mathieu sits on a tongue
  boundary, and the periodic integral action places a mode there by construction — so
  `stability_str` now takes the exponent and returns STABLE, MARGINAL or UNSTABLE.
  The eight call sites and the summary line report through it.

- **`templates/Periodic-Observer/h2_observer_lmi.m` printed the covariance negated.**
  `obj = -trace(P_0)` because the solver minimises; the report showed `value(obj)`. It
  now prints the trace itself. The label said "variance" while the dual template
  `periodic_kalman_lmi.m` gives that name to the *information* matrix; the two are now
  distinguished.

### Internal

No behaviour change; listed for reviewers.

- The phasor operands are named `pA1`..`pA5` rather than `o1`..`o5`, in the 19
  files that used the old names: 1423 occurrences, none outside `Fonctions/`.
  `pA` alone was already taken, so the numbering is kept even where a method
  takes a single operand.

- Adaptive harmonic-order refinement extracted from `lyap`, `lyapG` and
  `mlHmcDivide` into `adaptiveHSolve` (~500 duplicated lines).
- Toeplitz assembly prologue shared via `prepToeplitzInput`
  (`array2TBlocks`, `array2BToeplitz`, `sparray2TBlocks`).
- `isreal`/`isimag` share `projectionMatch`.
- Plotting moved out of `PhasorArray2time` (749 → 267 lines) into
  `TimeArray2plot`; most callers only want the evaluation. The visual behaviour is
  unchanged (dynamic grids, `ZeroCentered` symmetry, real/imaginary tiling through
  `manageTiledLayout`, `sgtitle`, global symmetric y-limits and `linkaxes`, the
  `linkprop` camera tie in 3-D, and the `angle (rad)` / `time (sec)` abscissa label).

  One deliberate difference: in the real-and-imaginary layout, the abscissa label
  now appears on both columns of the bottom row. It used to sit only on the
  imaginary one, leaving half the bottom row unlabelled.
- Name-value arguments normalised to the native `arguments` block with the `nvp.` prefix
  across the toolbox.
- ~300 lines of commented-out dead code removed, along with narrative comment blocks
  (ASCII flowcharts, "clean version" notes). What remains states the mathematics or a
  non-obvious algorithmic choice.
