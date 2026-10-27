# Figure interpretation guide (preliminary, not thesis-ready citations)

## 1. Effective mass and scalar profile
Plot the three terms `beta R/2`, `-alpha GB`, and their sum computed **on the GR background** in the left panel. Plot the fully coupled `phi_V` and a zero GR field in the right panel. The curvature is in km^-2; the GB contribution is multiplied by alpha [km²]. A locally negative effective mass squared does **not** automatically imply a growing normal mode. Compare with threshold figure 7.

## 2. M–R relation
At fixed EOS/beta/alpha and varying total central energy-equivalent density, black dashed curve is GR, solid curves are different alpha. The scalarized equilibrium sequence may bifurcate. Unverified configurations remain blank and are never joined. Turning points are **not** sufficient to establish radial stability in general modified-gravity models.

## 3. Scalar charge–mass relation
The exterior field is `phi = Q/r + M Q/r² + ...`, so `|Q/M|` is dimensionless. The GR branch has identically Q=0. Fitted charge uncertainties depend on outer radius and fit model. The code retains signed Q in CSV even though it plots |Q|.

## 4. Alpha dependence at fixed central energy
The top panel shows `100*(M_scalar - M_GR)/M_GR` at fixed EOS and total central energy density; the lower panel shows `|Q|/M_scalar`. GR is a dashed zero reference. Do not interpret changing alpha as keeping stellar baryon number fixed.

## 5. EOS dependence
Six small-multiple axes prevent 12 curves from overlapping. All panels use the same beta/alpha and the same central-energy-equivalent density grid. The shared axes allow direct visual comparison. No EOS table is inherently a physical truth; document uncertainty and data provenance.

## 6. Internal structure
At matching central total energy density, black dashed GR and blue solid scalarized pressure and energy-density radial profiles demonstrate hydrostatic changes. Near coincidence is a physically meaningful outcome, not necessarily a plotting error. Vertical dotted lines show (approximate) stellar surfaces. The surface is defined at the minimum **available positive tabulated pressure**.

## 7. Linear scalarization zero-mode threshold
Fundamental, no-node static scalar zero modes on GR stars yield threshold beta_critical(alpha), at fixed background EOS and central density. This is derived from a **linear eigenvalue equation**, not from visually locating m_eff²<0. A threshold curve is **not** a fully mapped nonlinear scalarized-solution existence region and is **not** a proof of dynamical stability.

## Caveats before thesis submission

- Repeat with at least two smaller tolerances, larger outer radius, smaller center cutoff, and sufficiently dense density/alpha grids; quote convergence errors on M, R and Q.
- Obtain the same analytical SLy/MPA1 EOS as the benchmark paper for a quantitative replication. The user's SLy4 table differs from the benchmark at about 1–2% in the published light-star mass check.
- Validate EOS continuity, causality, thermodynamic consistency and licensing, especially APR and small crust reversals in some tables.
- Add uncertainties to extracted ADM mass and scalar charge if comparing with observations; do not overlay observational constraints without a consistent observable calculation.

Scientific background: Ventagli, Antoniou, Lehébel and Sotiriou, *Phys. Rev. D* **104**, 124078 (2021), https://doi.org/10.1103/PhysRevD.104.124078.
