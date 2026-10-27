"""Run with python -m pytest -q; no GUI or network needed."""
import numpy as np
import pytest
from stgb.eos import load_eos,MEV_TO_GEO
from stgb.gravity import Couplings,algebraic_x,radial_rr_coeffs
from stgb.stellar import Settings,solve_gr,solve_scalar,observables,gr_curvature,shoot
from stgb.onset import find_critical_beta
from pathlib import Path

DATA=Path(__file__).resolve().parents[1]/'data'

@pytest.fixture(scope='module')
def sly():return load_eos('SLy4',DATA)
@pytest.fixture(scope='module')
def bg(sly):return solve_gr(sly,8.1e14,Settings())

def test_eos_units_and_monotonicity(sly):
    assert np.isclose(MEV_TO_GEO,1.323833e-6,rtol=1e-4)
    assert np.all(np.diff(sly.pressure)>0)
    assert sly.filtered_rows==0

def test_rr_branch_gr():
    r=9.;u=.03;p=.00012;c=Couplings(alpha=20,beta=-10)
    x=algebraic_x(r,u,0.,0.,p,c)
    assert np.isclose(x,(1+r*u)/(1+8*np.pi*r*r*p),rtol=2e-13)
    A,B,C=radial_rr_coeffs(r,u,0,0,p,c)
    assert abs(4*A*x*x+2*B*x+C)<2e-13

def test_gr_reference(bg):
    assert abs(bg['M_solar']-1.13779)<0.0002
    assert abs(bg['radius']-11.7920)<0.01

def test_gr_recovered_from_zero_scalar(bg,sly):
    star=solve_scalar(sly,8.1e14,0.,Couplings(30.,-10.),Settings())
    obs=observables(star)
    assert abs(obs['M_solar']-bg['M_solar'])<1e-5
    assert abs(star['radius']-bg['radius'])<1e-4
    assert obs['Q_km']==0

def test_schwarzschild_gauss_bonnet(bg):
    r=np.array([bg['radius']*3,bg['radius']*5]);R,G=gr_curvature(bg,r)
    assert np.allclose(R,0.,atol=1e-14)
    assert np.allclose(G,48*bg['mass']**2/r**6,rtol=1e-12)

def test_linear_zero_mode_near_literature_threshold(bg):
    b=find_critical_beta(bg,0.,-9.,-2.,samples=8)
    assert -6<b<-4

def test_nontrivial_nodeless_branch(sly):
    star,log=shoot(sly,8.1e14,Couplings(30,-10),Settings(shoot_points=13))
    assert star is not None
    assert star['nodes']==0
    assert abs(star['fit']['phi_inf_fit'])<3e-5
    assert 0<abs(star['fit']['Q_over_M'])<1
