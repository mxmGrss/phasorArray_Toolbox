# Reproducing the ECC paper

The reference is the submitted V2 (late March 2026), named
`ECC_PhasorArray.tex` at the root of the Overleaf project. The presentation
and subsequent unpublished local changes are separate works.

`Exemples/ECC_ex.m` preserves the paper's time-domain A(t), period T=1,
FFT setting N=6, B=[1;sin], Q=10I, R=1 and LMI orders. It uses the current
public Riccati entry point `hare`. The random algebra example uses seed 0
for reproducibility; the paper does not fix a seed.

## Explicit differences and acceptance

- The initial coefficient listing has a sign erratum: A_3(1,1) must be
  -4/(3*pi)^2 for the triangular wave specified by A(t). This is corrected
  in the example. The next construction overwrites A from the unchanged
  time function, so subsequent control calculations are unaffected.
- A(t) is singular at t=1/8. The division example is retained, with the
  diagnostics of `mlHmcDivide`, the method underlying `A\B`. Its failure to
  reach the requested residual is reported, not hidden by changing A or
  loosening a threshold. `inv(A)` is sampled inversion, not a certificate
  of a regular inverse between samples.
- Lyapunov and Riccati must report convergence. The Riccati relative
  residual must be at most 1e-6.
- Closed-loop exponents evaluated at h=20 must remain in the left half-plane
  and within 1e-3 of the printed values -3.4466 and -2.3234. This is a
  reproduction tolerance, not a claim of identical historical rounding.
- YALMIP must report `problem=0` for the finite LMI. This is not by itself
  a positivity certificate for an infinite-dimensional operator.

On MATLAB R2025b, the unmodified paper listings produced a Riccati relative
residual 3.59e-8, S order 73, K order 74 and exponents approximately
-3.4465 and -2.3237. The printed orders 60 and 61 are not fixed acceptance
targets: adaptive solver implementations can take different paths.

Completion without exceptions is insufficient: the algebraic division is
an explicitly documented limitation. Do not describe this example as proving
convergent inversion for the paper's singular A(t).
