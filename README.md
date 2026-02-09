# PhasorArray Toolbox

**Harmonic Modeling, Analysis, and Control of Periodic Systems made simple.**

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/mxmGrss/phasorArray_Toolbox/graphs/commit-activity)

The **PhasorArray Toolbox** is an object-oriented library that takes the pain out of analyzing and controlling **Linear Time-Periodic (LTP)** and **Bilinear Systems**. 

By leveraging **Harmonic Modeling** and **Toeplitz Matrix** algebra, it effectively turns complex periodic differential equations into standard algebraic problems. This means you can apply robust LTI control techniques—like Eigenvalue analysis or $H_\infty$ synthesis—to systems that were previously too difficult to model (e.g., Power Converters, Rotating Machinery).
**Built on Solid Science**: This toolbox enables the numerical application of advanced theoretical results presented in recent control literature [1-8], ensuring **consistency** and **guarantees** where ad-hoc methods often fail.

---

## 🚀 Why This Toolbox?

If you deal with periodic systems, you know the struggle: time-varying matrices, infinite harmonics, and unstable numerical methods. 

**PhasorArray solves this by:**

*   **Algebraicizing the Periodic**: Converts $\dot{x}(t) = A(t)x(t)$ into $\dot{X} = (\mathcal{A} - \mathcal{N})X$. Suddenly, your periodic problem looks just like an LTI one.
*   **Guaranteeing Consistency**: We implement distinct, rigorous algorithms (Riedinger et al.) for truncation [6]. This ensures that your harmonic operations actually converge to the true infinite-dimensional solution ($\mathcal{T}(AB)_h$ to $\mathcal{T}(A)\mathcal{T}(B)$).
*   **Targeting Real Applications**:
    *   **LTP Systems**: Floquet analysis and stabilization made easy.
    *   **Bilinear Systems**: Design **Forwarding Controllers** for global stability on Power Converters (see [2], [4]).
    *   **Variable Frequency**: Use **Dynamic Phasors** to control **PMSM** drives even during speed transients (see [3]).

### Based On & Related To
This framework encompasses and generalizes methods known in the literature as:

*   **Harmonic State Space (HSS)** Modeling
*   **Dynamic Phasors** & **Generalized State Space Averaging (GSSA)**
*   **Harmonic Balance** Method
*   **Linear Time-Periodic (LTP)** & **Polynomial Parameter-Periodic (PPP)** Systems

---

## ✨ Key Features

### 1. The `PhasorArray` Class

The core of the toolbox. It hides the complexity of Fourier series and Toeplitz matrices behind a clean, intuitive object.

*   **Natural Syntax**: Add, multiply, invert periodic matrices as if they were constants: `C = A * B + inv(D)`.
*   **Smart Constructors**:
    *   `PhasorArray(A)`: Build from raw data.
    *   `funcToPhasorArray(fh, T, N)`: Pass a function handle `@(t) ...`, get a PhasorArray (auto-FFT).
    *   `time2Phasor(t, x)`: Convert time-series data directly.
*   **Flexible Forms**: Switch between `SinCos`, `AngleAmp`, or `RealImag` representations instantly.
*   **Reduction Tools**: Use `neglect(tol)` to automatically prune negligible harmonics and keep your models light.

### 2. Advanced Solvers

Don't just model—solve.

*   **Eigenvalues**: Compute **Floquet exponents** with `HmqNEig` to determine stability.
*   **Lyapunov**: Assess stability and compute periodic Gramians with `lyap`.
*   **Riccati (HARE)**: Synthesize **LQR** and **$H_{\infty}$** controllers using efficient iterative algorithms (`RicHarmonicKlein`).
*   **Sylvester**: Solve observer design problems with `sylvester`.

### 3. Robust Control (LMI)

Seamless integration with **YALMIP** allows you to harness the power of convex optimization.

*   **`ndsdpvar`**: Create periodic decision variables effortlessly.
*   **Toeplitz-Block LMIs**: Formulate and solve $H_2$/$H_{\infty}$ problems to design controllers that are **robust** across entire operating ranges.

---

## 🛠️ Installation

1.  **Clone** the repository:
    ```bash
    git clone https://github.com/mxmGrss/phasorArray_Toolbox.git
    ```

2.  **Open MATLAB** and navigate to the folder.

3.  **Run the installer**:
    ```matlab
    installToolbox
    ```
    This script will set up your path and automatically open the documentation.

---

## ⚡ Quick Start Example

Here is a minimal example showing how to create a periodic matrix, analyze its stability (Floquet), and solve a Lyapunov equation.

```matlab
% 1. Create a Periodic Matrix A(t)
% Define via function handle or harmonic coefficients
T = 1; % Period
w = 2*pi/T;
At = @(t) [0, 1; -1 - 0.5*cos(w*t), -0.2]; % Harmonic oscillator with varying stiffness
A = PhasorArray.funcToPhasorArray(At, T, 10); % 10 harmonics

% 2. Analyze Stability (Floquet Exponents)
% Compute eigenvalues of the harmonic operator (A - N)
eigenvalues = A.HmqNEig(10, T, "fundamental");
disp('Floquet Exponents:');
disp(eigenvalues);

% 3. Solve Harmonic Lyapunov Equation: (A-N)'P + P(A-N) + Q = 0
Q = PhasorArray.eye(2); 
P = lyap(A, Q, "T", T);

% 4. Visualize the Result
figure;
plot(P); 
title('Periodic Solution P(t) of Lyapunov Equation');
```

---

## 🔬 Applications & Scientific Context

### Grid-Tied Power Converters [2], [4]

For **Bilinear Systems** like AC/DC converters, the toolbox enables the design of **Forwarding Controllers**. 

*   **Problem**: Standard controllers often fail to guarantee stability under large disturbances.
*   **Solution**: Our approach guarantees **Global Stability** and allows for precise **Harmonic Mitigation** (Active Filtering) by directly shaping the closed-loop harmonic spectrum.

### Variable Speed Drives & PMSM [3]

The toolbox extends to **LPV (Linear Parameter-Varying)** models using **Dynamic Phasors**.

*   **Application**: Permanent Magnet Synchronous Motors (PMSM).
*   **Goal**: Control torque and speed while actively rejecting current harmonics caused by distortions, even as the speed varies.
*   **Method**: The speed $\omega(t)$ becomes a scheduling parameter in the harmonic model: 
    $$\dot{X} = (\mathcal{A} - \omega(t)\mathcal{N})X + \mathcal{B}U$$

---

## 📚 Documentation

Detailed documentation for all classes and functions is available directly within MATLAB:

```matlab
doc PhasorArray
doc PhasorSS
```

See the `Exemples/` folder for complete scripts covering:
*   Periodic LQR design.
*   LMI-based reference tracking.
*   Modeling of DC/DC converters.

---

## 🔗 References

### Toolbox & Applications
1.  **[Toolbox]** M. Grosso, P. Riedinger, and J. Daafouz, *"The PhasorArray Toolbox for Harmonic Analysis and Control Design,"* arXiv preprint arXiv:2510.21294, 2025.
2.  **[TCST]** M. Grosso, P. Riedinger, J. Daafouz, S. Pierfederici, H. J. Idrissi, and B. Lapôtre, *"Harmonic Control of Three-Phase AC/DC Converter With Time-Domain Guarantees,"* *IEEE Transactions on Control Systems Technology*, 2025.
3.  **[ECCE]** M. Grosso, P. Riedinger, J. Daafouz, S. Pierfederici, H. J. Idrissi, and B. Lapôtre, *"Frequency-Varying Harmonic Domain Control for PMSMs with Current Harmonic Mitigation,"* in *2025 IEEE Energy Conversion Conference Congress and Exposition (ECCE)*, pp. 1-8, 2025.
4.  **[IECON]** M. Grosso, P. Riedinger, J. Daafouz, S. Pierfederici, H. J. Idrissi, and B. Lapôtre, *"Control of three-phase PWM rectifier using multiple frame dq transform and harmonic modeling,"* in *IECON 2024-50th Annual Conference of the IEEE Industrial Electronics Society*, 2024.

### Theoretical Foundations
5.  **[TAC 24]** F. Vernerey, P. Riedinger, and J. Daafouz, *"A TBLMI Framework for Harmonic Robust Control,"* *IEEE Transactions on Automatic Control*, 2025 (Erratum 2025).
6.  **[TAC 22]** P. Riedinger and J. Daafouz, *"Solving infinite-dimensional harmonic Lyapunov and Riccati equations,"* *IEEE Transactions on Automatic Control*, vol. 68, no. 10, pp. 5938-5953, 2022.
7.  **[TAC 21]** N. Blin, P. Riedinger, J. Daafouz, L. Grimaud, and P. Feyel, *"Necessary and sufficient conditions for harmonic control in continuous time,"* *IEEE Transactions on Automatic Control*, vol. 67, no. 8, pp. 4013-4028, 2021.
8.  **[CDC 23]** F. Vernerey, P. Riedinger, and J. Daafouz, *"On solving infinite-dimensional Toeplitz block LMIs,"* in *2023 62nd IEEE Conference on Decision and Control (CDC)*, pp. 4717-4722, 2023.

---

*Developed by Maxime Grosso, Pierre Riedinger, and Jamal Daafouz at CRAN (Centre de Recherche en Automatique de Nancy).*
