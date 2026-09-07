# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> AI agent context for `phasorArray_Toolbox`. Read this before touching any file.

## Project Overview

MATLAB OOP toolbox for **harmonic modeling, analysis, and control of Linear Time-Periodic (LTP) and Bilinear systems**. Implements the bijection $L^2_{loc} \leftrightarrow \mathcal{H}$ between time-domain signals and their harmonic (Fourier/Toeplitz) representation, enabling standard LTI tools on periodic systems.

**Authors:** Maxime Grosso, Pierre Riedinger, Jamal Daafouz (CRAN, Université de Lorraine)
**Requirement:** MATLAB R2021b+ (no extra toolboxes for core features)
**License:** MIT — DOI: 10.5281/zenodo.17560958

---

## Repository Structure

```
phasorArray_Toolbox/
├── Fonctions/
│   ├── @PhasorArray/             ← Core class (38 files, 8500 lines, 223 methods;
│   │                                PhasorArray.m alone is 4900 lines)
│   ├── @PhasorSS/                ← Periodic state-space (LTP/LPV/LTV)
│   ├── @sparsePhasorArray/       ← Sparse variant
│   ├── pArrayBasicOperations/    ← Computational kernels (68 files)
│   ├── Display and data manipulation/  ← Visualization (18 files)
│   └── SimulationTools/          ← Floquet/simulation utilities (3 files)
├── Exemples/                     ← 7 application examples + GettingStarted.m
├── docs/                         ← Unified documentation (LaTeX, Wiki assets)
├── templates/                    ← Control design templates (ACDC, LQR, SPMSM…)
├── installToolbox.m              ← Path setup entry point
├── checkDependencies.m           ← Dependency checker
├── REQUIREMENTS.md               ← Toolbox dependency matrix
└── CITATION.cff                  ← Software citation metadata
```

---

## Installation & Setup

```matlab
% In MATLAB — run once per machine
run('installToolbox.m')

% Verify installation
checkDependencies("verbose", true)
```

---

## Running Tests

The suite is `matlab.unittest`: seven `TestCase` classes in `tests/`, run through
`run_all_tests.m`.

```matlab
results = run_all_tests();            % full regression suite
results = run_all_tests("install");   % Install-tagged smoke set only
```

The two answer different questions. The full suite is the regression net and
exercises paths a user never touches directly — fallback kernels, symbolic
payloads, solver residuals. The install set is one check per layer, no optional
toolbox, a couple of seconds: a red install run means the *installation* is
wrong, a red full run means the *code* is.

A test joins the smoke set by moving into a tagged block:

```matlab
methods (Test, TestTags = {'Install'})
```

`tests/` holds `PhasorArrayCoreTest`, `PhasorArrayCalculusTest`,
`PhasorArrayHarmonicOperatorsTest`, `PhasorArraySolversTest`,
`PhasorArraySimulationTest`, `PhasorArrayTimeDomainTest` and
`PhasorArrayCompatibilityTest`.

`Fonctions/test_PhasorArray_basic.m` and `Fonctions/test_PhasorArray_advanced.m`
(the earlier struct-returning runners) were removed on 2026-09-07: nothing
called them, `run_all_tests` did not execute them, and their coverage was
subsumed by `tests/`, verified theme by theme.

---

## Core Architecture

### Class Hierarchy

| Class | Role |
|---|---|
| `@PhasorArray` | 3D array `[n × m × (2h+1)]` — periodic matrix in harmonic domain |
| `@PhasorSS` | State-space `{A, B, C, D}` of `PhasorArray` + LPV parameter `p` |
| `@sparsePhasorArray` | Memory-efficient sparse variant |

### PhasorArray Storage Convention

A `PhasorArray` of dimension `[n × m]` truncated at harmonic order `h` is stored as a `[n × m × (2h+1)]` double array, where slice `k` = harmonic coefficient `k - h - 1` (i.e., center slice = DC component).

### Key Computational Kernels (`pArrayBasicOperations/`)

| File | Role |
|---|---|
| `PhasorArrayTimes.m` | Convolution multiplication via `tensorprod` (R2022a+) |
| `PhasorArrayTimes2.m` | Fallback multiplication (R2021b, matrix-based) |
| `SylvHarmonic.m` | Harmonic Sylvester solver — **also the Lyapunov path**, as a special case |
| `RicHarmonicKlein.m` | Iterative Riccati solver (adaptive h, LQR fallback) |
| `adaptiveHSolve.m` | Shared adaptive-h driver — `lyap`, `lyapG`, `mlHmcDivide` and `place` all route through it, each supplying its own `solveOnce` / `buildPreamble` / `computeResidual` / `packInfo` |
| `array2TBlocks.m` / `array2BToeplitz.m` | Harmonic array → Toeplitz-Block (TB) / Block-Toeplitz (BT) operators, harmonics ascending (−h..+h) in both |

### Version Compatibility

- **R2021b**: uses `PhasorArrayTimes2.m` (slower, Kronecker-based)
- **R2022a+**: uses `PhasorArrayTimes.m` (`tensorprod`, accelerated)
- Detection is automatic in `PhasorArray.m`

---

## Branches

A three-stage promotion pipeline. The name says the stage, not a claim about
stability.

| Branch | Stage |
|---|---|
| `dev` | Work in progress. No stability or compatibility guarantee. |
| `stable` | Passed review, waiting out a soak period before being versioned. Created on the first promotion; does not exist yet. |
| `main` | Tested, released, tagged, and the only branch carrying the Zenodo DOI. |

Promotion gates:

- `dev` -> `stable`: `run_all_tests()` green, the examples tour green, CHANGELOG up to date.
- `stable` -> `main`: a soak period with no issue, version bump, tag, GitHub release, DOI, `CITATION.cff`.

Do **not** merge into `main` without explicit instruction.

`stable_experimental` was deleted in August 2026: it had been an exact
duplicate of `main` throughout, and the pipeline it was meant to embody had
never actually run.

---

## Code Conventions

### Method Return Patterns

Solvers (`lyap`, `lyapG`, `mlHmcDivide`, `place`) return a structured `info`
output. Fields produced by every one of them:

```matlab
[X, info] = lyap(PA, Q)
% info.status      0=CONVERGED, 1=STAGNATED, 2=MAXH_REACHED, 3=FIXED_H,
%                  4=UNREACHABLE (extrapolation says the threshold is out of reach)
% info.statusMsg   human-readable form of the status
% info.h           harmonic order the solution was accepted at
% info.resnorm     absolute residual norm
% info.resrelnorm  relative residual norm
% info.residualPhasor
% info.h_history info.res_history info.resrel_history
% info.time_history info.regime_history
% info.s_alg_history info.s_exp_history
```

Plus two fields the adaptive driver computes for every solver:

```matlab
% info.hForTargetResidual  order a near-zero residual would need, extrapolated
% info.targetResidual      the residual that extrapolation aimed at
%                          (both NaN when the solver ran at a fixed order)
% info.solskewnorm         ||(X - X')/2||_F — how far a solution that should be
%                          Hermitian actually is. A norm, so it compares
%                          directly with resnorm. NaN when X is not square.
```

**The contract is defined once**, in
`pArrayBasicOperations/packSolverInfo.m`. Every solver calls it — including the
fixed-order branches, which pass the values they have and let the helper fill
the rest with the documented defaults, so the field set never varies. Add a
guaranteed field there and every solver gains it; a field only one solver
reports goes through the `extra` argument, which refuses to overwrite a
guaranteed one.

`testInfoContractIsIdenticalAcrossSolvers` (in `tests/PhasorArraySolversTest.m`)
asserts the exact field set for all six solver entry points.

### Error / Warning Policy

Use `warning('PhasorArray:...')` with toolbox-specific IDs — **not** `fprintf` or `disp` for diagnostics.

### Harmonic Order `h`

- Always pass `h` explicitly when constructing fixed-order arrays.
- Adaptive solvers increase `h` until convergence or `maxH` is reached.

### MATLAB Style

- `arguments` blocks for all public functions
- `mustBe*` validators for input validation
- No global variables
- `applyStyle` for plot formatting

---

## Optional Dependencies

| Dependency | Unlocks |
|---|---|
| Control System Toolbox (R2022b+) | `toLPVss`, `toLTVss` export |
| Symbolic Math Toolbox | Symbolic `PhasorArray` |
| YALMIP | `ndsdpvar`, LMI synthesis |
| MOSEK / SDP solver | Fast SDP resolution |
| Signal Processing Toolbox | **Removed** as of March 2026 |

---

## Agent Instructions

- **Before refactoring `@PhasorArray`**: map the dependencies first — 223 methods with non-obvious interdependencies. If the repo carries a `.codegraph*/` index, `codegraph impact <symbol>` answers this in one call.
- **Tests are `matlab.unittest`** — add them to `tests/`, as methods of the relevant `TestCase` class. Do not add to the legacy `Fonctions/test_PhasorArray_*.m` runners.
- **`LyapHarmonic.m` is superseded and unused** — the Lyapunov path is `@PhasorArray/lyap.m` → `SylvHarmonic`, which handles truncation properly (`'rectangle'` vs `'square'`); `LyapHarmonic` has no such notion and its Riccati mode is an `error(...)` stub. Do not build on it. Slated for removal.
- **Untracked directories are the developer's own** — if a working sandbox exists locally, do not clean it up autonomously.
- **Do not add `Signal Processing Toolbox` calls** — it was intentionally removed.
- **CI/CD workflows are absent** — `.github/` is tracked and reserved for public GitHub config (PR template today, `workflows/` when CI lands).

- **This file is committed.** It must stand on its own for anyone who clones the
  repository: never cite an audit report, a working note, or any path that is
  gitignored or untracked. If a fact only lives in a local document, state the
  fact here instead of pointing at the document.
