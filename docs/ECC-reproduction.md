# ECC harmonic-control example

Run the example from the toolbox root:

```matlab
run('Exemples/ECCpaperAsOfV2_0.m')
```

The script covers Fourier construction, harmonic arithmetic, state-space modeling, Lyapunov and Riccati equations, feedback design and simulation on a two-state periodic system.

## Model and solver settings

| Parameter | Value |
|---|---|
| Period | `T=1` |
| Fourier sampling | `N=6` |
| Control input | `B=[1;sin(2*pi*t/T)]` |
| LQR weights | `Q=10*eye(2)`, `R=1` |
| Initial gain | `K0=[10,10]` |
| Riccati solve | Initial order 6, maximum order 500, maximum 50 iterations, residual tolerance `1e-6` |
| Finite LMI orders | 20/10/10 |

The triangle and square waves are expressed with base MATLAB formulas. The state-space model has two outputs, one input and a 2-by-1 feedthrough matrix.

## Harmonic arithmetic

Inversion and division use the regular matrix

$$A_{\mathrm{alg}}(\theta)=\begin{bmatrix}2&0.2\\0.2&1.5\end{bmatrix}+0.25\cos(\theta)I.$$

Its smallest eigenvalue is at least 1.05 at every phase. A seeded random right-hand side (`rng(0)`) makes the division example reproducible. Both the harmonic-division residual and the inverse residual must be below `1e-8`.

## Control results

The Riccati solve uses `hare` with adaptive harmonic order. The script checks a relative residual below `1e-6` and closed-loop exponents within `1e-3` of `[-3.4466,-2.3234]` at analysis order 20.

The LMI calculation checks successful solution of the finite section with YALMIP. For pointwise periodic positivity, follow the [Gram certificate example](https://github.com/mxmGrss/phasorArray_Toolbox/wiki/Interfacing-with-YALMIP).

## Correspondence with the article

| File | Contents |
|---|---|
| [ECCpaperAsOfV2_0.m](../Exemples/ECCpaperAsOfV2_0.m) | Executable example for toolbox v2, including residual and convergence checks |
| [ECC_ex.m](../Exemples/ECC_ex.m) | The article's 18 MATLAB listings in their published order |

The executable example retains the time-domain plant, control weights, initial condition and simulation grid. It uses a corrected third Fourier coefficient `A_3(1,1)=-4/(3*pi)^2`, dimensionally consistent feedthrough, and the regular matrix above for arithmetic. The control plant is singular at `t=1/8`, so its pointwise inverse is undefined there.

The published Riccati call uses `RicHarmonicKlein` with `warmStartFraction=0.95`; the executable example uses `hare` with its default `1.0`. Adaptive orders depend on solver settings and are reported with the residuals.

The listing archive corresponds to `ECC_PhasorArray.tex`, submitted version V2, March 2026 (SHA256 `050539dd0f072abefeb632bb18cc03dbeecb5879fcac9f0d2a5d3aded4f84ef0`).
