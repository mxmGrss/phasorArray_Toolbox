# Reproducing the ECC paper

The published-paper reference is `ECC_PhasorArray.tex` at the root of the
Overleaf archive confirmed by the author (submitted V2, late March 2026).
The presentation and later unpublished local changes are separate works.
Reference TeX SHA256: `050539dd0f072abefeb632bb18cc03dbeecb5879fcac9f0d2a5d3aded4f84ef0`.

## Published listings: ECC_ex.m

`Exemples/ECC_ex.m` concatenates the 18 MATLAB listing bodies verbatim and
in order, including the commented listing. Only a provenance header and
section separators are added. It retains the printed third coefficient,
`A\B`, the unseeded random matrix and the direct `RicHarmonicKlein` call with
`warmStartFraction=0.95`. No correction or acceptance assertion is inserted.

This preserves the published source, not the historical toolbox runtime:
random draws, unspecified defaults and adaptive solver implementations can
change across versions. In particular, execution without exceptions is not
proof that every numerical operation converges.

## Adaptation to toolbox v2.0.0: ECCpaperAsOfV2_0.m

`Exemples/ECCpaperAsOfV2_0.m` is the separate, instrumented adaptation.
The underscore makes the version suffix a valid MATLAB script name.
It preserves the time-domain A(t), T=1, FFT N=6, control B=[1;sin], Q=10I,
R=1, K0=[10,10], Riccati h=6/maxh=500/maxIter=50/threshold=1e-6,
LMI orders 20/10/10, and simulation initial condition and grid.
Intentional differences from the published source are:

- Correct A_3(1,1) to -4/(3*pi)^2, consistent with the specified triangle.
  The next construction overwrites A from the unchanged time function, so
  this correction does not alter subsequent control calculations.
- Set rng(0) for the random algebra illustration; the paper specifies no seed.
- Replace A\B with its underlying mlHmcDivide call to retain diagnostics and
  explicitly report nonconvergence.
- Use the public hare entry point. Its warmStartFraction is the current
  solver default 1.0, whereas the published call explicitly uses 0.95.
  This is an algorithmic adaptation, not identical solver parametrization.
- Add Lyapunov/Riccati convergence checks, a Riccati residual bound of 1e-6,
  a closed-loop exponent comparison within 1e-3 of [-3.4466,-2.3234] at h=20,
  and a YALMIP problem=0 check for the finite LMI.

## Known numerical limitations

A(t) is singular at t=1/8. The algebra illustration therefore does not
establish a regular inverse; inv(A) is sampled inversion and does not certify
invertibility between samples. The harmonic-division residual must not be
interpreted as successful convergence merely because the script completes.
The LMI status concerns the finite problem, not an additional certificate
for the infinite-dimensional operator.

The earlier execution of the published listings on R2025b yielded Riccati
relative residual 3.59e-8, S/K orders 73/74 and closed-loop exponents near
-3.4465 and -2.3237. The printed orders 60/61 are historical results, not
fixed acceptance criteria for a newer adaptive solver.