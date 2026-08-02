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
│   ├── @PhasorArray/             ← Core class (6500+ lines, 150+ methods)
│   ├── @PhasorSS/                ← Periodic state-space (LTP/LPV/LTV)
│   ├── @sparsePhasorArray/       ← Sparse variant
│   ├── pArrayBasicOperations/    ← Computational kernels (56 files)
│   ├── Display and data manipulation/  ← Visualization (18 files)
│   └── SimulationTools/          ← Floquet/simulation utilities (3 files)
├── Exemples/                     ← 8 application examples + GettingStarted.m
├── docs/                         ← Unified documentation (LaTeX, Wiki assets, audits)
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

Custom lightweight test runner — **not** `matlab.unittest.TestCase`.

```matlab
% Basic test suite (~40 tests: constructors, arithmetic, indexing, Toeplitz, reduction)
results = test_PhasorArray_basic();

% Advanced test suite (~20 tests: symbolic, YALMIP, Floquet, LMI — auto-skips missing toolboxes)
results = test_PhasorArray_advanced();

% Inspect failures
failed = results(~[results.passed]);
disp({failed.name; failed.message}')
```

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
| `SylvHarmonic.m` | Harmonic Sylvester equation solver |
| `LyapHarmonic.m` | Harmonic Lyapunov solver |
| `RicHarmonicKlein.m` | Iterative Riccati solver (adaptive h, LQR fallback) |
| `array2TBlocks.m` / `array2BToeplitz.m` | Harmonic array → Toeplitz-Block (TB) / Block-Toeplitz (BT) operators, harmonics ascending (−h..+h) in both |

### Version Compatibility

- **R2021b**: uses `PhasorArrayTimes2.m` (slower, Kronecker-based)
- **R2022a+**: uses `PhasorArrayTimes.m` (`tensorprod`, accelerated)
- Detection is automatic in `PhasorArray.m`

---

## Active Side-Project Branches

| Branch | Status |
|---|---|
| `stable_experimental` | Pre-main consolidation branch. Contains V3 Floquet solvers and unified docs. |

Do **not** merge these into `main` without explicit instruction.

---

## Code Conventions

### Method Return Patterns

Solvers (`lyap`, `mlHmcDivide`, `RicHarmonicKlein`) return a structured `info` output:
```matlab
[X, info] = lyap(PA, Q)
% info.status: 0=CONVERGED, 1=STAGNATED, 2=MAXH_REACHED, 3=FIXED_H
% info.h_final, info.iterations, info.residual
```

### Error / Warning Policy

Use `warning('PhasorArray:...')` with toolbox-specific IDs — **not** `fprintf` or `disp` for diagnostics. See `docs/audit/disp_warning_audit.md` for migration status.

### Harmonic Order `h`

- Always pass `h` explicitly when constructing fixed-order arrays.
- Adaptive solvers increase `h` until convergence or `maxH` is reached.

### MATLAB Style

Follow `.agents/matlab-guidelines.md` for:
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

- **Before refactoring `@PhasorArray`**: run `audit-matlab` to map dependencies — the class has 150+ methods with non-obvious interdependencies.
- **Tests are custom structs**, not `matlab.unittest` — do not migrate without discussion.
- **`scratch/` is an untracked dev sandbox** — if the user creates one locally, do not clean it up autonomously.
- **Do not add `Signal Processing Toolbox` calls** — it was intentionally removed.
- **CI/CD workflows are absent** — `.github/` is tracked and reserved for public GitHub config (PR template today, `workflows/` when CI lands); agent instruction modules live in the untracked `.agents/`.
