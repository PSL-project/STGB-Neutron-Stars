"""Numerical regression checks for the faster radial solver.

Reference values were measured using the unoptimized October 2026 package,
with SLy4, rho_c=8.1e14 g/cm^3 and Settings() defaults. They are not
independently validated astrophysical benchmark values.
"""
import numpy as np
import pytest
from pathlib import Path
from stgb.eos import load_eos
from stgb.gravity import Couplings, rhs, regular_center
from stgb.stellar import Settings, solve_scalar, observables

DATA = Path(__file__).resolve().parents[1] / 'data'

@pytest.fixture(scope='module')
def sly():
    return load_eos('SLy4', DATA)

@pytest.mark.parametrize('alpha, expected_radius, expected_mass, expected_q', [
    (-200., 10.818030615730107, 0.7317843100319652, 0.09504929160889714),
    (0., 11.648952745676981, 1.0798399173233897, 0.08153320384324786),
    (30., 11.752163337574263, 1.106475445817328, 0.08176698199088396),
    (100., 11.957569453815818, 1.1401293870562956, 0.08375154720469179),
])
def test_original_implementation_regression(sly, alpha, expected_radius, expected_mass, expected_q):
    # Absolute tolerances accommodate adaptive ODE step-size differences, while
    # remaining far tighter than the scientific accuracy of these models.
    star = solve_scalar(sly, 8.1e14, .02, Couplings(alpha, -10.), Settings())
    obs = observables(star)
    assert abs(star['radius'] - expected_radius) < 3e-4  # km
    assert abs(obs['M_solar'] - expected_mass) < 4e-6
    assert abs(obs['Q_over_M'] - expected_q) < 4e-6


def test_fast_radial_solver_residual(sly):
    couplings = Couplings(30., -10.)
    star = solve_scalar(sly, 8.1e14, .02, couplings, Settings())
    for r in (0.001, .01, .2, 1., 3., 6., 10., 11.5, 20., 90., 300.):
        if r < star['radius']:
            y = star['interior'].sol(r)
            ep = sly.epsilon
        else:
            y = star['exterior'].sol(r)
            ep = lambda p: 0.
        deriv, diagnostics = rhs(r, y, ep, couplings, diagnostics=True)
        assert np.all(np.isfinite(deriv))
        assert abs(diagnostics['rreq']) < 1e-10
        assert diagnostics['linear_resid'] < 1e-10


def test_exact_gr_central_series_without_nonlinear_root(sly):
    eps, p = sly.central(8.1e14)
    c = Couplings(100., -10.)
    _, (g, l, s) = regular_center(eps, p, 0., c)
    assert s == 0.
    assert np.isclose(l, c.kappa * eps / 3., rtol=1e-14)
    assert np.isclose(g, c.kappa * (eps + 3*p) / 6., rtol=1e-14)
