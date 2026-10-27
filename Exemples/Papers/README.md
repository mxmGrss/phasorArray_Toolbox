# Worked examples: periodic analysis and control

Each script solves a specified problem with the current toolbox API. The models, conventions and expected results below allow the examples to be used without reading the source articles. The references at the end give the scientific background and the correspondence with published calculations.

Install the toolbox on `dev`, then run a script from the toolbox root. Start with the scalar Lyapunov example; use Riccati or pole placement for feedback design, and the Gram examples for optimization-based inequalities.

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

## Shared two-state model

The Riccati, pole-placement and Gram examples use `dx/dt = A(t)x + B(t)u`, with two states, one input and period `T = 1`. The helper `riedingerVernereyPlant(p,3,tau)` constructs the following finite Fourier model; all its data are included in the toolbox:

$$
A(t)=A_0+\sum_{k=1}^{p}\left(A_k e^{jk\omega t}+\overline{A_k}e^{-jk\omega t}\right),
\qquad A_0=\begin{bmatrix}1&2\\-1&1\end{bmatrix},\qquad \omega=2\pi/T.
$$

For positive `k`, with `o_k = 1` for odd `k` and zero otherwise, the entries of `A_k` are:

| Entry | Fourier coefficient |
|---|---|
| (1,1) | `-2i*o_k/(pi*k)` |
| (1,2) | `8*tau*o_k/(pi*k)^2` |
| (2,1) | `(-1)^k/(1i*pi*k)` |
| (2,2) | `1i` at k=1, `1+1i` at k=3, `1` at k=5, zero otherwise |

The input matrix is `B(t) = [1 + 2*cos(2*omega*t) + 4*sin(3*omega*t); 0]`. Riccati and Gram use `p=15, tau=1`; placement uses `p=20, tau=0.5`. No phase offset is added to entry (2,1).

`p` fixes the plant degree; `h` or `q` fixes the degree of the unknown solution. Increasing the solution degree leaves the plant unchanged. `T` is the time period; phase evaluation uses `theta = 2*pi*t/T`. The certificates below concern this finite polynomial plant.

## What each script computes

### Solve a periodic Lyapunov equation

For the scalar system `dx/dt = a(t)x`, take `a(t) = -1-cos(t)+2*sin(t)+cos(2*t)` and `T=2*pi`. Find periodic P satisfying

$$\dot P(t)+2a(t)P(t)+1=0.$$

`lyap` computes a least-squares approximation at each harmonic order. The script reports the residual of the full equation, including coefficients outside the solve window. The last residual should be about `1.18e-11` at order 24. This example needs only MATLAB and the toolbox. Scientific background: [1](#references-and-correspondence-with-the-articles).

### Compute and certify a periodic LQR gain

Use the shared plant with `Q=diag(1,10000)` and `R=1` in the cost integral of `x'*Q*x + u'*R*u`. The Riccati equation is

$$\dot P+A^\top P+PA-PBR^{-1}B^\top P+Q=0,\qquad u=-Kx,\quad K=R^{-1}B^\top P.$$

`hare` computes the candidate at degree 120. `certifyHare` independently checks the closed-loop stability and an error bound in scaled coordinates `z=diag(1,64)*x`; the returned P and K are in the original coordinates. Expect mean trace `1595.6055`, `certificate.accepted=true` and a bound near `8.50e-6`. Control System Toolbox may be used for gain initialization. Scientific background: [1](#references-and-correspondence-with-the-articles).

### Assign closed-loop Floquet exponents

Use the shared plant with `p=20, tau=0.5`. The Floquet decomposition expresses `x=V(t)z` with constant dynamics `dz/dt=J*z`. The script chooses the target matrix `Lambda=-J'-I` and `G=B'/V'`, then solves

$$\dot P-AP+P\Lambda+BG=0,\qquad K=GP^{-1},\qquad u=-Kx.$$

`place` refines the solve order and checks invertibility of the coordinate change. Compare the exponents of `A-B*K` with `eig(Lambda)`: they should be approximately `-1.99993 ± 1.64821i`, with target error below `1e-3`. Symbolic Math Toolbox is required by `FloquetDec`. Scientific background: [2](#references-and-correspondence-with-the-articles).

### Formulate a Riccati inequality with real Gram constraints

Use the shared plant and the same LQR weights as above. Choose a real symmetric polynomial P of degree 40, maximize its average trace, and impose positivity of P and of

$$M(t)=\begin{bmatrix}\dot P+A^\top P+PA+Q&PB\\B^\top P&R\end{bmatrix}.$$

`gramConstraint` turns these continuous-symbol inequalities into a finite SDP. After optimization, `checkGramCertificate` validates the numerical result for each inequality. Expect two accepted certificates and mean trace about `1590.9614`. The checks certify the displayed inequalities; they do not assert that the finite-degree candidate equals the optimal Riccati solution. The state scaling `diag(1,100)` improves conditioning, and K is returned for `u=-Kx` in the original coordinates. Requires YALMIP and MOSEK. Scientific background: [3](#references-and-correspondence-with-the-articles).

### Bound disturbance amplification with H-infinity control

Extend the shared plant to `dx/dt=A(t)x+B(t)u+w`, with two disturbance inputs and performance output `z=[x(1);100*x(2);u]`. Seek feedback `u=K(t)x` whose induced L2 gain from w to z is bounded by `gamma=24.922`.

The SDP uses a positive polynomial P and the variable Y=KP, both of degree 20, with state scaling `diag(1,10)`. The block inequality is written explicitly in the listing. Both numerical Gram checks must be accepted; their lower bounds should be near `9.96e-4` and `1.00e-3`. The returned `K(theta)` evaluates the rational gain at a scalar phase `theta=2*pi*t/T` in the original state coordinates. Requires YALMIP and MOSEK. Scientific background: [3](#references-and-correspondence-with-the-articles).

### Find a rectifier operating point and simulate startup

The states are the AC inductor current i (A) and DC voltage v (V). The modulation s is limited to [-1,1]. With `R=1e-3 ohm`, `L=1e-3 H`, `C=4e-3 F`, `RL=10 ohm` and a 100 V, 50 Hz source,

$$L\dot i=-Ri-sv+100\sin(2\pi50t),\qquad C\dot v=si-v/R_L.$$

The supplied helper fits a degree-5 periodic operating point with mean voltage 200 V, then the script solves a Lyapunov equation and simulates saturated feedback from zero current and voltage. Expect mean operating voltage `199.9998 V`, full operating-point residual about `0.0479`, and startup voltage about `198.026 V` at 0.2 s. The full residual exposes the operating-point truncation error. Requires Optimization Toolbox. Scientific background: [4](#references-and-correspondence-with-the-articles).

### Build a variable-speed PMSM model

The plant states are three phase currents and mechanical angular speed. The script specifies resistance, inductance, flux, inertia, friction and four pole pairs, then builds A(theta), Bu and the output maps. Mechanical phase satisfies `dtheta/dt=omega`, with `omega` between 10 and 200 rad/s.

Ten output integrators augment the four plant states to 14 states. `Cnegative` and `Cpositive` expose the two q-axis sign conventions; `Qfull` weights all augmented states and `Qintegrators` only the integrators. At each speed vertex, the harmonic state operator is the Toeplitz matrix of the augmented symbol minus `omega*N`, where N uses a `2*pi` phase period. Inspect the matrices and the energy-balance assertion. Only MATLAB and the toolbox are required. This example stops at model construction; controller design is a separate step.

## Interpretation and limits

- The infinite Fourier series in A are explicitly truncated to p. Certificates apply to this polynomial model. They are not automatically transferred to the original discontinuous coefficients.
- The LQR Gram checks validate the displayed Riccati inequality and positivity of P. They do not certify equality with the Riccati solution or global optimality of the finite-degree gain.
- Gram solver status and certificate acceptance are checked separately. Bounds use floating point, not interval arithmetic. LQR normalizes the scalar objective by 10000 for numerical conditioning without changing its optimizer.
- The H-infinity certificate concerns the rational gain Y(t)/P(t). K is a function handle evaluated at phase theta, using numeric matrix division after reconstruction. This avoids an additional Fourier approximation of the inverse.
- The rectifier uses the paper’s order-5 square equilibrium problem and reports its full residual. A refined Lyapunov solve and startup simulation are included; disturbance-rejection and forwarding figures are not reproduced by this short script.
- The placement example refines the Sylvester order on the p=20 polynomial model; it does not reproduce all original-data trajectories or the entire m=4,6,8,10 comparison.

## Dependencies

- MATLAB core for scalar Lyapunov and PMSM model.
- Control System Toolbox may be used by the Riccati stabilizing-gain initialization.
- Symbolic Math Toolbox for the existing FloquetDec used by placement.
- YALMIP and MOSEK for the Gram examples; the q=20 H-infinity default is a short case. The paper’s q=60 case is substantially heavier and requires its own numerical certificate validation.
- Optimization Toolbox for the rectifier equilibrium.

## References and correspondence with the articles

The models, equations and expected outputs are specified above. These references give the derivations and broader numerical studies. Section and equation numbers below identify the corrected formulations used for these examples; numbering can differ between versions.

1. P. Riedinger and J. Daafouz, [Solving Infinite-Dimensional Harmonic Lyapunov and Riccati equations](https://arxiv.org/abs/2203.09774). Scalar Lyapunov: §5.3, Table 1; Riccati: §7, Eq.35, Table 2.
2. P. Riedinger and J. Daafouz, [Harmonic Pole Placement](https://arxiv.org/abs/2203.09140), §IV-A, Eq.21.
3. F. Vernerey, P. Riedinger and J. Daafouz, [A TBLMI Framework for Harmonic Robust Control](https://arxiv.org/abs/2311.05934). LQR: §VI-A, Eq.44; H-infinity: §VI-B, Eq.30.
4. Blin et al., *Necessary and Sufficient Conditions for Harmonic Control in Continuous Time*, IEEE TAC 67(8), 4013–4028, 2022; rectifier: §6.1–6.2.

These compact scripts cover selected calculations; the scope of each example is stated above.

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
