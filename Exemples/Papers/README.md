# Periodic analysis and control examples

Run a script from the toolbox root after installation:

```matlab
run('Exemples/Papers/riedinger_lyapunov.m')
```

The [worked examples](https://github.com/mxmGrss/phasorArray_Toolbox/wiki/Paper-Examples) explain the models, equations and results.

| Problem | Script | Main result |
|---|---|---|
| Solve a periodic Lyapunov equation | [riedinger_lyapunov.m](riedinger_lyapunov.m) | Full equation residual versus harmonic order. |
| Compute and certify an LQR gain | [riedinger_riccati.m](riedinger_riccati.m) | Feedback K, candidate P and an independent error bound. |
| Assign Floquet exponents | [riedinger_pole_placement.m](riedinger_pole_placement.m) | Feedback K and error relative to the target exponents. |
| Impose a Riccati inequality with Gram | [vernerey_lqr.m](vernerey_lqr.m) | Polynomial candidate and numerical inequality certificates. |
| Bound disturbance amplification | [vernerey_hinf.m](vernerey_hinf.m) | Rational feedback gain and H-infinity certificates. |
| Simulate a rectifier startup | [blin_rectifier.m](blin_rectifier.m) | Periodic operating point and startup trajectory. |
| Build a variable-speed motor model | [automatica_pmsm_model.m](automatica_pmsm_model.m) | Plant, output maps and augmented models. |
| Explore a full harmonic-control workflow | [ECCpaperAsOfV2_0.m](../ECCpaperAsOfV2_0.m) | Executable ECC example; [reproduction notes](../../docs/ECC-reproduction.md). |

## Dependencies

| Script | Additional software |
|---|---|
| Scalar Lyapunov, PMSM model | MATLAB and the toolbox |
| Riccati | Control System Toolbox may be used for gain initialization |
| Pole placement | Symbolic Math Toolbox (`FloquetDec`) |
| LQR and H-infinity Gram | YALMIP and MOSEK |
| Rectifier | Optimization Toolbox |

Keep the scripts in this folder: MATLAB resolves the plant coefficients in [riedingerVernereyPlant.m](private/riedingerVernereyPlant.m) and the rectifier objective in [blinRectifierResidual.m](private/blinRectifierResidual.m) through the `private` subfolder.

## Parameters and outputs

Riccati and Gram examples use a two-state, single-input Fourier plant with period `T=1` and degree `p=15`. Pole placement uses `p=20` and a triangle factor of `0.5`. The solution degree is set separately by `h` or `q`.

Riccati, placement and LQR use feedback `u=-K*x`. H-infinity uses `u=K*x`; its function handle `K(theta)` evaluates the rational gain at `theta=2*pi*t/T`. Returned gains act on the original state coordinates.

The Gram bounds apply to the finite Fourier plant and the displayed polynomial inequalities. `checkGramCertificate` includes coefficient residuals and the numerical Gram eigenvalue in a floating-point lower bound. Check `accepted` after every solve, including after changing the degree or margin.

## Numerical checks

Reference results obtained with MATLAB R2025b on 27 October 2026:

| Calculation | Checked result |
|---|---|
| Scalar Lyapunov, h=24 | Full residual 1.1767e-11. |
| Riccati, h=120 | Mean trace 1595.6055015; certificate bound 8.5013e-6 in scaled coordinates. |
| Placement | Poles approximately -1.99993 ± 1.64821i; target error 1.38e-13 on the polynomial model. |
| LQR Gram, q=40 | Mean trace 1590.9613773; validated lower margins 9.46e-5 and 2.11e-2. |
| H-infinity Gram, q=20 | gamma=24.922; validated lower margins 9.96e-4 and 1.00e-3. |
| Rectifier | Equilibrium mean 199.9998 V; full equilibrium residual 0.0479; Lyapunov relative residual 1.01e-13. At 0.2 s the startup voltage is 198.026 V. |
| Automatica model | Both output conventions execute; dimensions and physical energy balance checked.  |

## References

- Riedinger and Daafouz, [Solving Infinite-Dimensional Harmonic Lyapunov and Riccati equations](https://arxiv.org/abs/2203.09774): Lyapunov and Riccati examples.
- Riedinger and Daafouz, [Harmonic Pole Placement](https://arxiv.org/abs/2203.09140): Floquet exponent assignment.
- Vernerey, Riedinger and Daafouz, [A TBLMI Framework for Harmonic Robust Control](https://arxiv.org/abs/2311.05934): LQR and H-infinity examples. The scripts use solution degrees 40 and 20, respectively.
- Blin et al., *Necessary and Sufficient Conditions for Harmonic Control in Continuous Time*, IEEE TAC 67(8), 4013–4028, 2022: rectifier example.
