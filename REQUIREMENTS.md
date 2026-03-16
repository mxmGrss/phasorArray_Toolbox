# PhasorArray Toolbox — Requirements

## Overview

The PhasorArray Toolbox is designed to operate across a broad range of MATLAB versions, providing core functionality natively and utilizing optional toolboxes for extended capabilities.

| Dependency | Status | Version Requirements | Purpose |
| :--- | :---: | :--- | :--- |
| **MATLAB** | **Required** | **R2021b or later** | Core engine; foundational operations. |
| **Control System Toolbox** | Optional | R2022b or later | Export to `lpvss` / `ltvss` native objects. |
| **Symbolic Math Toolbox** | Optional | Any recent | Symbolic instantiation and derivation (`PhasorArray.sym()`). |
| **YALMIP** | Optional | Latest | LMI-based synthesis (`ndsdpvar`, Toeplitz-Block LMIs). |
| **MOSEK** / SDP Solver | Optional | Compatible with YALMIP | Fast numerical resolution of SDP problems. |

---

## Core Capabilities (R2021b+, MATLAB Only)

The foundations of the toolbox require base MATLAB (R2021b or later) with no additional dependencies.

### Native Functionality

- **PhasorArray Class**: Instantiation and manipulation of periodic matrices.
- **Arithmetic Operations**: Matrix operations (`+`, `*`, `\`, `/`).
- **Time-Domain Evaluation**: Methods such as `evalTime`, `evalp`, `initial`, `lsim`.
- **Harmonic Domain Analysis**: Toeplitz formalisms (`BT`, `TB`), Fourier operators (`F_tb`, `FvTB`).
- **Eigenvalue Analysis**: Computation of Floquet exponents via `HmqNEig`.
- **Solvers**: Harmonic Lyapunov (`lyap`), Riccati (`RicHarmonicKlein`), and Sylvester (`sylvester`) resolutions.
- **Visualization**: `plot`, `stem`, `barsurf`, and `spy` for structural analysis.

### Version-Specific Optimizations

The computational backend adapts based on the available MATLAB release:

#### R2022a+ (Recommended)
- Utilizes `tensorprod()` for accelerated harmonic convolution in `PhasorArrayTimes.m` and `PhasorArray2time.m`.
- Achieves significant performance improvements for large harmonic truncations.

#### R2021b (Minimum Supported)
- Defaults to matrix-based convolution (`PhasorArrayTimes2.m`).
- Numerically equivalent but computationally slower. Automatically detected via `isMATLABReleaseOlderThan("R2022a")`.

---

## Optional Dependencies

### 1. Control System Toolbox

**Required For:**
- Object conversion to MathWorks native formats (`PhasorSS.toLPVss()`, `PhasorSS.toLTVss()`).

**Notes:** Direct harmonic simulation via `PhasorSS.lsim()` or manual state-space evaluations remain functional without this toolbox. Invoking conversion methods without the toolbox will result in standard undefined function errors.

### 2. Symbolic Math Toolbox

**Required For:**
- The creation and manipulation of symbolic elements (`PhasorArray.sym(n, m, h, name)`), facilitating formal derivation of periodic systems.

**Notes:** Numerical workflows, encompassing the vast majority of applications, remain entirely unaffected by the absence of this toolbox.

### 3. YALMIP

**Required For:**
- Instantiating periodic decision variables (`PhasorArray.ndsdpvar(n, m, h)`).
- Formulating Toeplitz-Block Linear Matrix Inequalities (LMIs) for robust control.

**Installation:** Available via the [YALMIP GitHub repository](https://github.com/yalmip/YALMIP). Ensure the path is integrated into the MATLAB environment.

**Notes:** Direct algebraic solutions (Lyapunov, Riccati) remain functional in the absence of YALMIP.

### 4. SDP Solver (e.g., MOSEK)

**Required For:**
- Numerical resolution of Toeplitz-Block LMIs formulated via YALMIP.

**Alternatives:** While [MOSEK](https://www.mosek.com/) is recommended for performance, YALMIP natively supports free alternatives such as SeDuMi, SDPT3, and SDPA.

---

## Verification Suite

A built-in test suite is provided to verify the integrity of the installation and detect available capabilities.

```matlab
% Base functionality verification
test_PhasorArray_basic

% Advanced functionality verification (Symbolic, YALMIP integration)
test_PhasorArray_advanced
```

The test framework categorizes results as follows:
- **Passed**: Functionality is verified.
- **Skipped**: Optional toolboxes are absent; dependencies are unmet but non-critical.
- **Failed**: Core regression or unmet base requirements.

---

## Summary of Feature Dependencies

| Feature | MATLAB Core | Control Toolbox | Symbolic Math | YALMIP | SDP Solver |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **PhasorArray Instantiation** | Required | - | - | - | - |
| **Arithmetic Operators** | Required | - | - | - | - |
| **Time-Domain Simulation** | Required | - | - | - | - |
| **Floquet Analysis** | Required | - | - | - | - |
| **Lyapunov / Riccati Solvers** | Required | - | - | - | - |
| **Symbolic Derivation** | Required | - | Required | - | - |
| **LPV/LTV Export** | Required | Required | - | - | - |
| **LMI Synthesis** | Required | - | - | Required | Required |

---

## Citation & Support

For academic or research utilization, please cite:

```bibtex
@software{phasorarray2025,
  author = {Grosso, M. and Riedinger, P. and Daafouz, J.},
  title = {The PhasorArray Toolbox for Harmonic Modelling and Control of Periodic systems},
  year = {2025},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.17560958}
}
```

**Resources:**
- [Issue Tracker (GitHub)](https://github.com/mxmGrss/phasorArray_Toolbox/issues)
- [Documentation Wiki](https://github.com/mxmGrss/phasorArray_Toolbox/wiki)

---

**Last Updated:** March 2026 (v1.0)
