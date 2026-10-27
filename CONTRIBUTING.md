# Contributing to STGB-Neutron-Stars

Thank you for considering a contribution. This project explores neutron-star solutions and spontaneous scalarization in scalar–tensor theories with Ricci-scalar and Gauss–Bonnet couplings. It is research software: reproducibility and scientific correctness matter at least as much as runtime.

## Before contributing

- Read [README.md](README.md) and any relevant documentation in [`docs/`](docs/).
- Search existing issues and pull requests before opening a duplicate.
- For substantial changes to the mathematical model, numerical algorithms, or published figures, open an issue first so assumptions and validation criteria can be discussed.
- Be aware that the code is exploratory. Passing automated tests does not by itself validate the underlying physical model.

## Report a problem or propose a change

Open a GitHub issue using the appropriate template. For numerical or scientific discrepancies, include as much of the following as possible:

- Equation of state (EOS), including its source and version.
- Coupling parameters, including `alpha` (in km^2, where applicable) and `beta`.
- Central density or pressure, scalar boundary conditions, and relevant solver settings.
- Exact command, Python version, operating system, and complete error message, if any.
- Expected behavior, observed behavior, and—where possible—a minimal reproducible example.
- A comparison with the GR limit, an independent calculation, or a literature reference when making a physics claim.

Please report possible security vulnerabilities privately according to [SECURITY.md](SECURITY.md), not as public issues.

## Contributing code

1. Create a branch from the relevant development branch; avoid committing directly to `main`.
2. Keep changes focused and explain *why* they are needed.
3. Add or update tests for changed behavior. Bug fixes should ideally include a regression test.
4. Install requirements and run the test suite:

   ```bash
   python -m pip install -r requirements.txt
   python -m pip install pytest
   python -m pytest -q
   ```

5. Open a pull request explaining the changes, their scientific or technical motivation, the tests performed, and any limitations.

If a change improves speed, report the benchmark setup and compare **both** runtimes and numerical outputs with the previous implementation. Do not trade accuracy for speed without explicitly documenting and justifying the tolerance.

## Scientific and numerical changes

Changes to field equations, scalar-field normalization, boundary conditions, shooting algorithms, EOS interpolation, or convergence criteria need special care. Whenever applicable:

- Recover the general-relativistic limit for vanishing scalar field.
- Check regularity at the stellar center and asymptotic behavior at large radius.
- Compare masses, radii, scalar charges, and profiles against the previous implementation or an independent reference.
- Document units and conventions, especially `alpha`, `beta`, the scalar-field definition, and geometric units.
- Keep linear scalarization thresholds distinct from nonlinear branch calculations and stability claims.
- Explain numerical tolerances, convergence tests, and known limitations.

## Data, figures, and attribution

Do **not** add proprietary, copyrighted, or otherwise redistribution-restricted EOS tables, article figures, or other third-party material without verifying permission. Record the origin, applicable license, and citation for every new dataset. When redistribution is not permitted or is uncertain, provide an acquisition or preprocessing script instead of uploading the data.

Generated plots should label their axes, units, EOS, and coupling parameters. Clearly separate illustrative calculations from independently validated or observationally constrained results.

## Conduct and licensing

Follow [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md). Unless separately agreed and documented, contributions submitted for inclusion in this project's MIT-licensed code are intended to be distributable under the project's existing license. Do not submit third-party code incompatible with that license.

The maintainer may request revisions or decline a contribution when correctness, licensing, scope, or maintenance requirements are not met.
