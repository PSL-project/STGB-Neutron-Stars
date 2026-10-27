# Audit plan and numerical cross-checks

## Reproducible sanity tests

1. Run `python -m pytest -q`.
2. Verify the `phi_c=0` limit of the full radial solver against independent GR TOV on SLy4 at `rho_c=8.1e14` g/cm3.
3. Check `alpha=0` light-SLy zero-mode threshold in `(-6,-4)` for the supplied EOS.
4. Compare one scalarized configuration with increasing `rmax`, decreasing `r0`, smaller maximum step, tighter tolerances. Record M, Q, R changes.
5. Check EOS table inversion and warnings. Review filtered G3/IOPB crust rows and APR superluminal points.
6. Inspect all seven figures for units, dashed GR references, no connections across missing scalarized branches, smooth center, continuous stellar surface fields and zero asymptotic scalar constant.

## Paper-level status

The solver follows Ventagli et al. (2021) action, GR-connected rr root, and radial reduction. Prior exploratory runs are consistent with *qualitative* trends in published charge/mass curves, but not a pointwise numerical reproduction. The original SLy fits/data table must be matched for stronger validation.

## Physics equations used by the linear threshold scan

On the GR background (kappa=8*pi):
- Ricci scalar `R=kappa*(epsilon-3*p)`.
- GB invariant `G=48*m**2/r**6 - 16*kappa*m*epsilon/r**3 - 4*kappa**2*epsilon*p`.
- Static scalar equation `phi''+[2/r+(Gamma'-Lambda')/2]*phi'=exp(Lambda)*(beta*R/2-alpha*G)*phi`.
- Exterior Schwarzschild solution `phi=constant+C*log(1-2*M/r)` is used to infer the asymptotic constant exactly in the GR-background exterior.

A zero of the asymptotic constant with no radial nodes is a fundamental zero-mode **threshold**; this calculation alone does not prove nonlinear branch existence or energy stability.
