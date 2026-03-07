# PhasorArray Toolbox — Requirements

## Summary

The PhasorArray Toolbox is designed to work on a **wide range of MATLAB versions** with graceful degradation when optional features or toolboxes are not available.

| Requirement | Status | Version/Details | Notes |
|:---|:---:|:---|:---|
| **MATLAB** | **Required** | **R2021b+** | Core engine; older versions may work with reduced functionality |
| **Control System Toolbox** | Optional | Any recent | Only for `PhasorSS.toLPVss()` / `toLTVss()` methods |
| **Symbolic Math Toolbox** | Optional | Any recent | Only for symbolic PhasorArrays (`PhasorArray.sym()`) |
| **YALMIP** | Optional | Latest | For LMI-based control (`ndsdpvar`, Toeplitz-Block LMIs) |
| **MOSEK** (or other SDP solver) | Optional | Compatible with YALMIP | For fast SDP solving; YALMIP can use others |

---

## Core Functionality (R2021b+, MATLAB Only)

The **majority** of the toolbox requires **only base MATLAB** (R2021b or later) and **no additional toolboxes**.

### What Works Without Add-Ons

- **PhasorArray Class**: Create, manipulate, and operate on periodic matrices.
- **Arithmetic Operations**: Add, multiply, invert PhasorArrays (`+`, `*`, `\`, `/`, etc.).
- **Time-Domain Evaluation**: `evalTime`, `evalp`, `initial`, `lsim`.
- **Harmonic Domain Analysis**: Toeplitz formalism (`BT`, `TB`), Fourier operators (`F_tb`, `FvTB`).
- **Eigenvalue Analysis**: Floquet exponents via `HmqNEig`.
- **Solvers**: Harmonic Lyapunov (`lyap`), Riccati (`RicHarmonicKlein`), Sylvester (`sylvester`).
- **Visualization**: `plot`, `stem`, `barsurf`, `spy` for Toeplitz structure.
- **Examples**: Most scripts in `Exemples/` folder use only base MATLAB.

### MATLAB Version-Specific Optimizations

The toolbox automatically adapts to the MATLAB version:

#### R2022a+ (Recommended)

- Uses `tensorprod()` for **fast harmonic convolution** in `PhasorArrayTimes.m` and `PhasorArray2time.m`.
- ~10× speedup for large harmonic orders.

#### R2021b (Minimum Supported)

- Falls back to **matrix-based convolution** (`PhasorArrayTimes2.m`).
- Slower but fully functional.
- **Automatic detection** via `isMATLABReleaseOlderThan("R2022a")`.

---

## Optional Toolboxes

### 1. Control System Toolbox

**Required For:**
- `PhasorSS.toLPVss()` — Convert PhasorSS to `lpvss` object.
- `PhasorSS.toLTVss()` — Convert PhasorSS to `ltvss` object.

**Introduced In:** R2022b (for `lpvss` and `ltvss` classes).

**Note:** If you don't need LPV/LTV conversion to MathWorks' native format, you can use:
- Direct harmonic simulation with `PhasorSS.lsim()`.
- Manual state-space construction and time-domain evaluation.

**Without It:**
- Calling `toLPVss()` or `toLTVss()` will throw an error:
  ```matlab
  Undefined function 'lpvss' for input arguments of type 'function_handle'.
  ```

---

### 2. Symbolic Math Toolbox

**Required For:**
- `PhasorArray.sym(n, m, h, name)` — Create symbolic PhasorArrays.
- Symbolic derivation and formal manipulation of periodic systems.

**Test Coverage:**
- `test_PhasorArray_advanced.m` includes symbolic tests (gracefully skipped if unavailable).

**Without It:**
- Symbolic operations are disabled.
- Numerical workflows (99% of use cases) work perfectly.

---

### 3. ~~Signal Processing Toolbox~~ ✅ **No Longer Required**

**Status:** ✅ **Eliminated as of March 2026**

Previously, `padarray()` was used in a few internal functions. This has been replaced by a native `phasorPad()` implementation using only `cat()` and `zeros()`, which is:
- Compatible with **all MATLAB versions** (R2021b+)
- Compatible with **symbolic arrays** (`sym`)
- Compatible with **YALMIP decision variables** (`sdpvar`, `ndsdpvar`)

No action needed — this dependency has been removed.

---

### 4. YALMIP

**Required For:**
- `PhasorArray.ndsdpvar(n, m, h)` — Create periodic decision variables.
- Toeplitz-Block LMI formulation for robust control.
- Examples: `Exemple_Toolbox_LMI.m`, `SPMSM_template.m`.

**Installation:**
```matlab
% Download from GitHub
% https://github.com/yalmip/YALMIP
addpath(genpath('path/to/YALMIP'));
savepath;
```

**Test Coverage:**
- `test_PhasorArray_advanced.m` checks for YALMIP availability:
  ```matlab
  hasYALMIP = exist('sdpvar', 'file') == 2;
  ```

**Without It:**
- LMI-based synthesis is unavailable.
- Direct algorithms (Riccati, Lyapunov) still work.

---

### 5. MOSEK (or other SDP Solver)

**Required For:**
- **Fast** and **numerically robust** solving of Toeplitz-Block LMIs.

**Installation:**
- [MOSEK Academic License](https://www.mosek.com/products/academic-licenses/) (free for research).
- YALMIP automatically detects installed solvers.

**Alternatives:**
- **SeDuMi**, **SDPT3**, **SDPA** — All free, but slower for large problems.
- YALMIP will use whatever is available.

**Test Coverage:**
- `test_check_solvers()` reports which solvers YALMIP detects.

**Without It:**
- You can still formulate LMIs, but solving them will fail if no SDP solver is installed.

---

## Verification Script

Run the built-in dependency checker:

```matlab
% Basic functionality
test_PhasorArray_basic

% Advanced features (symbolic, YALMIP, LMI)
test_PhasorArray_advanced
```

The test suite will:
- ✅ **Pass** tests that are supported.
- ⏭️ **Skip** tests for missing optional toolboxes (with clear messages).
- ❌ **Fail** only if a required dependency (R2021b+) is missing.

**Example Output:**
```
  PhasorArray Toolbox — Advanced Test Suite
════════════════════════════════════════════
  [SKIP] Symbolic Math Toolbox not found — sym tests will be skipped
  [SKIP] YALMIP not found — ndsdpvar tests will be skipped

  ✓ Symbolic: Matrix construction — SKIPPED (missing toolbox)
  ✓ YALMIP: construct ndsdpvar — SKIPPED (missing toolbox)
  ✓ YALMIP: Lyapunov LMI stability — SKIPPED (missing toolbox)

════════════════════════════════════════════
  12 passed, 8 skipped, 0 failed
```

---

## Minimum Installation (No Toolboxes)

If you only have **MATLAB R2021b+** (no toolboxes), you can still:

1. Create and manipulate `PhasorArray` objects.
2. Solve Lyapunov and Riccati equations.
3. Compute Floquet exponents and stability margins.
4. Simulate LTP systems in time-domain.
5. Visualize harmonic content.
6. Run **all examples** in `Exemples/BasicToolbox.m` and `Exemples/MathieuPendulum.m`.

**You cannot:**
- Use LMI-based control (needs YALMIP + solver).
- Create symbolic PhasorArrays (needs Symbolic Math Toolbox).
- Convert to `lpvss`/`ltvss` (needs Control System Toolbox).

---

## Recommended Setup for Research Use

For the **full experience** as described in the ECC 2026 paper:

```matlab
% 1. MATLAB R2022a+ (for tensorprod speedup)
% 2. Install YALMIP
addpath(genpath('~/YALMIP'));
% 3. Install MOSEK (academic license)
% 4. Optional: Control System Toolbox (for lpvss/ltvss)
```

**One-Time Check:**
```matlab
test_PhasorArray_advanced;  % See what's available
```

---

## Summary Table: What Needs What?

| Feature | MATLAB | Control Toolbox | Symbolic | Signal Proc. | YALMIP | Solver |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| **PhasorArray basics** | ✅ | — | — | — | — | — |
| **Arithmetic (+, *, /, \)** | ✅ | — | — | — | — | — |
| **Time-domain sim** | ✅ | — | — | — | — | — |
| **Floquet analysis** | ✅ | — | — | — | — | — |
| **Lyapunov / Riccati** | ✅ | — | — | — | — | — |
| **Kronecker products** | ✅ | — | — | — | — | — |
| **Symbolic derivation** | ✅ | — | ✅ | — | — | — |
| **LPV/LTV export** | ✅ | ✅ | — | — | — | — |
| **LMI synthesis** | ✅ | — | — | — | ✅ | ✅ |

Legend:
- ✅ Required
- ⚠️ Minor dependency (can be worked around)
- — Not needed

---

## Citation & Support

If you use this toolbox in your research, please cite:

```bibtex
@software{phasorarray2025,
  author = {Grosso, M. and Riedinger, P. and Daafouz, J.},
  title = {The PhasorArray Toolbox for Harmonic Modelling and Control of Periodic systems},
  year = {2025},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.17560958}
}
```

For issues or questions:
- 📧 maxime.grosso@protonmail.com
- 🐛 [GitHub Issues](https://github.com/mxmGrss/phasorArray_Toolbox/issues)
- 📖 [Full Documentation (Wiki)](https://github.com/mxmGrss/phasorArray_Toolbox/wiki)

---

**Last Updated:** March 2026 (v1.0)
