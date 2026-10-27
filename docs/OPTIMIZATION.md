# Solver performance update (optimization branch)

This update changes the **numerical implementation**, not the Ricci–Gauss–Bonnet action, radial field equations, shooting search interval, or selected EOS interpolators.

## What changed

1. `src/stgb/gravity.py`: reuses algebraic rr coefficients instead of calculating them twice at every radial ODE evaluation; solves the sparse 3-by-3 system by algebraic elimination instead of allocating an array and invoking a generic dense solver. Ill-conditioned/elimination-degenerate cases fall back to `numpy.linalg.solve`. Full matrix diagnostics are still available.
2. `src/stgb/gravity.py`: recognizes the exact zero-scalar GR central series and removes unnecessary nonlinear center solves, plus one redundant first continuation step for nonzero central scalars.
3. `src/stgb/stellar.py`: separates exterior and interior step limits. Exterior integration uses an adaptive RK45 integrator with a 4 km maximum step, rather than the old 0.8 km cap. The interior cap is relaxed from 0.16 km to 0.3 km. `rtol`, `atol`, event detection, background model, and surface definition are unchanged. Pass `--exterior-max-step 0.8` to restore the old exterior cap for convergence testing.
4. `src/stgb/analysis.py`: adds **optional** multiprocessing of independent mass-radius EOS/coupling sequences. All output is returned to the parent process. Identical `(EOS, alpha)` sequences required by different figure families are computed only once. `--workers 1` (default) preserves sequential behavior, while `--workers 4` works well on a six-core Windows machine.
5. Adds numerical-regression tests against the outputs of the earlier release, and `RUN_FAST_WINDOWS.bat` for launching with four workers.

## Usage

In Command Prompt / VS Code terminal:

```bat
py -3 run.py --quality quick --figures all --workers 4
py -3 run.py --quality thesis --figures all --workers 4
py -3 run.py --quality thesis --figures 2,3,5 --workers 4 --no-show --save
py -3 -m pytest -q
```

Double-click `RUN_FAST_WINDOWS.bat` for quick quality with parallel EOS/alpha sequences and the usual interactive preview.

On Windows the worker processes are safely started from `run.py` using a `__main__` guard. Multiprocessing is for the independent **sequence** tasks (figures 2/3/5), not every single ODE integration. Figure 7 threshold scans remain sequential. Higher worker counts consume more RAM, and short jobs may run **slower** owing to process-start overhead.

## Measurements (same Linux container, Python/SciPy installation)

At `SLy4`, `rho_c=8.1e14 g/cm^3`, `beta=-10`, `alpha=30`, default 400 km outer radius, 11 shoot grid points:

| Workload | Earlier solver | Optimized solver | Speedup |
|---|---:|---:|---:|
| Single star with `phi_c=0.02` | ~0.143 s | ~0.060 s | ~2.4x |
| Node-zero scalar shoot (11 grid points) | ~2.96 s | ~1.25 s | ~2.4x |
| Two sequences, 2 central densities each, short outer radius | ~4.15 s with one worker | ~2.41 s with two workers | ~1.7x |

Timing results will differ across Windows PCs, EOS, branch parameters, and figure sets. Measurements are small runs rather than a formal performance-characterization campaign; the multiprocess comparison uses **the optimized solver** with 1 vs 2 workers.

## Numerical regression checks

- Original and optimized solvers were run independently for regular configurations with `alpha=-200, 0, 30, 100`, including `phi_c=0` and nonzero values. Across successful points in this sample, relative differences in ADM mass remained below `2e-7`; the original and optimized solvers agreed on center failures at `alpha=200` for the tested high `phi_c` choices.
- The 11-point shooting root at `(alpha, beta)=(30,-10)` shifted by less than `1e-9` in central scalar amplitude; the derived mass and scalar charge changes were likewise tiny.
- The new `tests/test_optimization.py` verifies output values against the previous solver with explicit tolerances, tests sparse-system residuals, and tests the zero-field central coefficients.
- A two-worker vs one-worker mass-radius sequence comparison returned numerically **identical** row dictionaries for the tested configurations.

**This is a performance/numerical-regression test, not an independent physics validation.** No claim is made that the chosen branch is dynamically stable, that all nonlinear solution branches have been found, or that the solver has been point-by-point benchmarked against the journal literature. Retain or tighten `rtol` and change radial step limits when doing a final convergence study. Results produced on another operating system may differ slightly due to floating-point arithmetic and solver adaptation.

## Data permissions

The MIT software license does not cover third-party EOS tables automatically. Check the licenses and provenance of datasets inside `data/` before making a new public release; deleting data from the HEAD of a public branch does not erase copies in its Git history.
