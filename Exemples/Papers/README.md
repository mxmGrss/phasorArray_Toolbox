# Examples aligned with the corrected manuscripts

Install the toolbox on `dev`, then run a script from the toolbox root. The examples specify their model coefficients, truncation orders and conventions below. Equation numbering refers to corrected formulations and may differ between article versions.

| Manuscript and equation | Script | Scope |
|---|---|---|
| Riedinger–Daafouz, Lyapunov/Riccati, corrected formulation, §5.3, Table 1 | [riedinger_lyapunov.m](riedinger_lyapunov.m) | Scalar formula and full least-squares residual versus order. |
| Same manuscript, §7, Eq.35, Table 2 | [riedinger_riccati.m](riedinger_riccati.m) | p=15, h=120, Q=diag(1,10000), R=1; independent certificate in z=diag(1,64)x. |
| Riedinger–Daafouz, Pole Placement, corrected formulation, §IV-A, Eq.21 | [riedinger_pole_placement.m](riedinger_pole_placement.m) | Floquet target and adaptive Sylvester solve on p=20. |
| Vernerey–Riedinger–Daafouz, TBLMI, corrected formulation, §VI-A, Eq.44 | [vernerey_lqr.m](vernerey_lqr.m) | Direct-P LQR, real Gram, p=15, q=40; not the q=200 row of Table I. |
| Same manuscript, §VI-B, Eq.30 | [vernerey_hinf.m](vernerey_hinf.m) | Bw=I, p=15, q=20, fixed suboptimal gamma=24.922; real Gram. The paper uses q=60. |
| Blin et al., Necessary and Sufficient Conditions, corrected formulation, §6.1–6.2 | [blin_rectifier.m](blin_rectifier.m) | Single-phase rectifier: equilibrium optimization and saturated-feedback startup. |
| Variable-frequency PMSM model | [automatica_pmsm_model.m](automatica_pmsm_model.m) | Plant and augmented models with two Park conventions and weighting choices; model construction only. |
| ECC toolbox article | [ECCpaperAsOfV2_0.m](../ECCpaperAsOfV2_0.m) | Existing executable API adaptation; [reproduction notes](../../docs/ECC-reproduction.md). |

## Run

```matlab
run('Exemples/Papers/riedinger_lyapunov.m')
run('Exemples/Papers/riedinger_riccati.m')
run('Exemples/Papers/riedinger_pole_placement.m')
run('Exemples/Papers/vernerey_lqr.m')
run('Exemples/Papers/vernerey_hinf.m')
run('Exemples/Papers/blin_rectifier.m')
run('Exemples/Papers/automatica_pmsm_model.m')
```

The main scripts contain the actual toolbox calls. [riedingerVernereyPlant.m](private/riedingerVernereyPlant.m) holds only analytic plant data; [blinRectifierResidual.m](private/blinRectifierResidual.m) expresses the weighted equilibrium objective and its harmonic solve. No private research toolbox or Signal Processing Toolbox is used.

## Data and conventions that matter

- All three corrected Riedinger/Vernerey models use input harmonic **3**, with **no phase offset in a21**. Pole Placement has triangular amplitude **8/pi²**; Lyapunov/Riccati and TBLMI use **16/pi²**. The shared helper exposes this difference.
- The corrected scalar Lyapunov formula is `-1-cos(t)+2*sin(t)+cos(2*t)`.
- Riccati and LQR return K for **u=-Kx**. H-infinity uses **u=Kx**.
- The infinite Fourier series in A are explicitly truncated to p. Certificates apply to this polynomial model. They are not automatically transferred to the original discontinuous coefficients.
- The LQR Gram checks validate the displayed Riccati inequality and positivity of P. They do not certify equality with the Riccati solution or global optimality of the finite-degree gain.
- Gram solver status and certificate acceptance are checked separately. Bounds use floating point, not interval arithmetic. LQR normalizes the scalar objective by 10000 for numerical conditioning without changing its optimizer.
- The H-infinity certificate concerns the rational gain Y(t)/P(t). K is a function handle evaluated at phase theta, using numeric matrix division after reconstruction. This avoids an additional Fourier approximation of the inverse.
- The rectifier uses the paper’s order-5 square equilibrium problem and reports its full residual. A refined Lyapunov solve and startup simulation are included; disturbance-rejection and forwarding figures are not reproduced by this short script.
- The placement example refines the Sylvester order on the p=20 polynomial model; it does not reproduce all original-data trajectories or the entire m=4,6,8,10 comparison.

## Variable-frequency PMSM model

The model uses mechanical omega in [10,200] rad/s and four pole pairs. It exposes two output conventions: `Cnegative` uses the negative-sine q row of Park, and `Cpositive` uses the positive-sine q row. The cost matrices `Qfull` and `Qintegrators` respectively weight all augmented states or only the integrators. Choose one consistent convention when designing a controller.

The vertex symbols require the harmonic term -omega*N, with N defined for a 2*pi phase period. This example checks model dimensions and physical energy balance; it does not synthesize a controller or reproduce closed-loop performance figures. The generic SPMSM template uses different model parameters.

## Dependencies

- MATLAB core for scalar Lyapunov and PMSM model.
- Control System Toolbox may be used by the Riccati stabilizing-gain initialization.
- Symbolic Math Toolbox for the existing FloquetDec used by placement.
- YALMIP and MOSEK for the Gram examples; the q=20 H-infinity default is a short case. The paper’s q=60 case is substantially heavier and requires its own numerical certificate validation.
- Optimization Toolbox for the rectifier equilibrium.

## Publication references

- Riedinger–Daafouz: [Lyapunov/Riccati](https://arxiv.org/abs/2203.09774), [Pole Placement](https://arxiv.org/abs/2203.09140).
- Vernerey–Riedinger–Daafouz: [TBLMI](https://arxiv.org/abs/2311.05934).
- Blin et al.: Necessary and Sufficient Conditions for Harmonic Control in Continuous Time, IEEE TAC 67(8), 4013–4028, 2022.

These compact scripts reproduce selected calculations, not every figure, tracking scenario or alternative LQR/H2 formulation.

## Validation on MATLAB R2025b (27 October 2026)

| Calculation | Checked result |
|---|---|
| Scalar Lyapunov, h=24 | Full residual 1.1767e-11; the Table 1 residual sequence is reproduced. |
| Riccati, h=120 | Mean trace 1595.6055015; certificate bound 8.5013e-6 in scaled coordinates, matching the reported Table 2 row. |
| Placement | Poles approximately -1.99993 ± 1.64821i; target error 1.38e-13 on the polynomial model. |
| LQR Gram, q=40 | Mean trace 1590.9613773; validated lower margins 9.46e-5 and 2.11e-2. |
| H-infinity Gram, q=20 | gamma=24.922; validated lower margins 9.96e-4 and 1.00e-3. |
| Rectifier | Equilibrium mean 199.9998 V; full equilibrium residual 0.0479; Lyapunov relative residual 1.01e-13. At 0.2 s the startup voltage is 198.026 V: this finite startup is not exact convergence to the periodic reference. |
| Automatica model | Both output conventions execute; dimensions and physical energy balance checked. Model construction only. |

The H-infinity default uses q=20 and a margin of 0.001 at gamma=24.922. Changing the order or margin requires rerunning the numerical certificate checks: solver success alone does not establish positivity. Different settings can yield different gains; the example does not regenerate the paper’s q=60 gain plot.
