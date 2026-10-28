"""Independent code-facing STGB radial system for Ventagli et al. (2021).

Uses their canonical scalar normalization, static spherical metric
 ds^2=-exp(Gamma) dt^2+exp(Lambda) dr^2+r^2 dOmega^2,
 kappa = 8*pi, gamma_derivative_coupling=0, scalar_mass=0.

STATE y=(Gamma, U=Gamma', phi, psi=phi', p).
Lambda is NOT integrated: x=exp(Lambda) is the GR-connected algebraic
root of the rr field equation. Lambda', U', psi' are obtained by solving
three linear relations: tt, scalar, and derivative of the rr constraint.

The derivation is independent of the typographical problem in the printed
Appendix A3 scalar equation; the scalar row uses Box(phi)=(beta R/2-alpha G)phi.

Requires numpy and scipy; EOS callback eps_of_p(p) must return geometrized
energy density (km^-2) given geometrized pressure (km^-2).
A stand-alone EOS+stellar structure+shooting+plotting script is included
below. Verified against its GR limit and several sample scalarized runs;
NOT YET quantitatively benchmarked against the published paper results.
"""
from dataclasses import dataclass
from typing import Callable
import numpy as np
from scipy.optimize import root

KAPPA = 8.0*np.pi

@dataclass(frozen=True)
class Couplings:
    alpha: float     # km^2
    beta: float      # dimensionless; paper's beta
    kappa: float = KAPPA


def radial_rr_coeffs(r, u, phi, psi, p, c: Couplings):
    a,b,k=c.alpha,c.beta,c.kappa
    A=1+k*r*r*p-0.5*b*k*phi**2
    B=(-2+b*k*phi**2-2*r*u+r*b*k*phi**2*u
       +4*r*b*k*phi*psi-8*a*k*phi*u*psi
       +r*r*b*k*phi*u*psi+k*r*r*psi**2)
    C=48*a*k*phi*u*psi
    return A,B,C


def _rr_root_from_coeffs(A, B, C, r):
    """Recover exp(Lambda) from already-computed rr coefficients."""
    D=B*B-4*A*C
    if not np.isfinite(D) or D < 0 or abs(A)<1e-14:
        raise ValueError(f'Nonreal/singular rr branch at r={r}, A={A}, D={D}')
    x=(-B+np.sqrt(D))/(4*A)
    if not np.isfinite(x) or x<=0:
        raise ValueError(f'Unphysical exp(Lambda)={x} at r={r}')
    return x


def algebraic_x(r,u,phi,psi,p,c):
    return _rr_root_from_coeffs(*radial_rr_coeffs(r,u,phi,psi,p,c),r)


def rhs(r,y,eps_of_p:Callable[[float],float],c:Couplings,diagnostics=False):
    Gam,u,phi,psi,p=map(float,y)
    a,b,k=c.alpha,c.beta,c.kappa
    if r<=0: raise ValueError('Use center series: r must be positive.')
    eps=float(eps_of_p(max(p,0.0))) if p>0 else 0.0
    # Matter exterior is vacuum; scalar and metric remain nontrivial.
    pdot=-0.5*(eps+p)*u if p>0 else 0.0
    A,B,C=radial_rr_coeffs(r,u,phi,psi,p,c)
    x=_rr_root_from_coeffs(A,B,C,r)
    # Row 1, tt: tL*Lambda' + tPsi*phi'' = -t0.
    tL=x*(-8*a*k*phi*psi+b*k*r*r*phi*psi+b*k*r*phi**2-2*r)+24*a*k*phi*psi
    tPsi=16*a*k*phi*(x-1)-2*b*k*r*r*phi*x
    t0=(x*x*(b*k*phi**2+2*k*r*r*eps-2)
        +x*((16*a*k-2*b*k*r*r+k*r*r)*psi**2-b*k*phi**2-4*b*k*r*phi*psi+2)
        -16*a*k*psi**2)
    # Row 2, scalar: sL*Lambda'+sU*Gamma''+phi''=-s0.
    R0=2*(1-1/x)/r**2-(u*u/2+2*u/r)/x
    GB0=-2*(x-1)*u*u/(x*x*r*r)
    sL=-psi/2-phi*(b/2*(u/2+2/r)-2*a*(x-3)*u/(x*r*r))
    sU=phi*(b/2-4*a*(x-1)/(x*r*r))
    s0=(2/r+u/2)*psi-x*phi*(b*R0/2-a*GB0)
    # Row 3, preserve the rr constraint E=4A*x*x+2B*x+C=0.
    Ar=2*k*p*r
    Ap=-b*k*phi
    B_r=(u*b*k*phi**2+2*u*b*k*phi*psi*r-2*u+4*b*k*phi*psi+2*k*psi**2*r)
    B_phi=k*(-8*u*a*psi+2*u*b*phi*r+u*b*psi*r*r+2*b*phi+4*b*psi*r)
    B_psi=k*(-8*u*a*phi+u*b*phi*r*r+4*b*phi*r+2*psi*r*r)
    B_u=(-8*a*k*phi*psi+b*k*phi**2*r+b*k*phi*psi*r*r-2*r)
    C_phi=48*u*a*k*psi
    C_psi=48*u*a*k*phi
    C_u=48*a*k*phi*psi
    E_r=4*x*x*Ar+2*x*B_r
    E_phi=4*x*x*Ap+2*x*B_phi+C_phi
    E_p=4*x*x*k*r*r
    cL=x*(8*A*x+2*B)
    cU=2*x*B_u+C_u
    cPsi=2*x*B_psi+C_psi
    c0=E_r+E_phi*psi+E_p*pdot
    # Eliminate Lambda' analytically and solve the remaining 2x2 system.
    # This matrix has the exact sparse pattern
    #   [tL,  0, tPsi]
    #   [sL, sU,    1]
    #   [cL, cU, cPsi]
    # and is evaluated thousands of times inside solve_ivp. Avoiding a generic
    # numpy.linalg.solve call and two temporary arrays significantly reduces
    # Python overhead without modifying the field equations.
    use_fallback = abs(tL) <= 1e-12 * max(1., abs(tPsi))
    if not use_fallback:
        L0=-t0/tL
        L1=-tPsi/tL
        h1=-s0-sL*L0
        h2=-c0-cL*L0
        k1=1.+sL*L1
        k2=cPsi+cL*L1
        det=sU*k2-cU*k1
        scale=abs(sU*k2)+abs(cU*k1)
        use_fallback=(not np.isfinite(det) or
                      abs(det)<=2e-13*max(scale,1e-300))
        if not use_fallback:
            Uprime=(h1*k2-k1*h2)/det
            psiprime=(sU*h2-cU*h1)/det
            Lprime=L0+L1*psiprime
            if not (np.isfinite(Lprime) and np.isfinite(Uprime) and
                    np.isfinite(psiprime)):
                use_fallback=True
    if use_fallback:
        mat=np.array([[tL,0.,tPsi],[sL,sU,1.],[cL,cU,cPsi]])
        vec=np.array([-t0,-s0,-c0])
        try:
            Lprime, Uprime, psiprime=np.linalg.solve(mat,vec)
        except np.linalg.LinAlgError as exc:
            raise ValueError(f'Singular radial system at r={r}') from exc
    deriv=np.array([u,Uprime,psi,psiprime,pdot])
    if diagnostics:
        mat=np.array([[tL,0.,tPsi],[sL,sU,1.],[cL,cU,cPsi]])
        vec=np.array([-t0,-s0,-c0])
        return deriv, {'expLambda':x,'LambdaPrime':Lprime,
                       'rreq':4*A*x*x+2*B*x+C,
                       'linear_resid':np.max(np.abs(mat@np.array([Lprime,Uprime,psiprime])-vec)),
                       'cond':np.linalg.cond(mat)}
    return deriv


def center_coefficients(eps0,p0,phi0,c:Couplings):
    """Solve the 3 central series equations for Gamma_2,Lambda_2,phi_2.

    Gamma=Gamma0+Gamma2 r^2+..., Lambda=Lambda2 r^2+...,
    phi=phi0+phi2 r^2+...,  p=p0-0.5*(eps0+p0)*Gamma2*r^2+...
    Roots must be tracked by continuation when varying couplings.
    """
    a,b,k=c.alpha,c.beta,c.kappa
    lg=k*eps0/3
    gg=k*(eps0+3*p0)/6
    sg=phi0*(b*(lg-gg)/2+4*a*gg*lg)
    def fun(v):
        g,l,s=v
        return [96*a*k*l*phi0*s+3*b*k*l*phi0**2-12*b*k*phi0*s+2*eps0*k-6*l,
                64*a*g*k*phi0*s+2*b*g*k*phi0**2-b*k*l*phi0**2+8*b*k*phi0*s-4*g+2*k*p0+2*l,
                2*s+b*(g-l)*phi0-8*a*g*l*phi0]
    result=root(fun,[gg,lg,sg],tol=1e-11)
    if not result.success and max(np.abs(fun(result.x)))>1e-9:
        raise ValueError('No reliable GR-connected center root from this initial guess.')
    return result.x


def initial_state(r0,eps0,p0,phi0,c,Gamma0=0):
    g,l,s=center_coefficients(eps0,p0,phi0,c)
    return np.array([Gamma0+g*r0**2,2*g*r0,phi0+s*r0**2,
                     2*s*r0,p0-0.5*(eps0+p0)*g*r0**2]),(g,l,s)



def regular_center(eps_c,p_c,phi_c,c,r0=1e-3):
    # Continuation of the regular central series along phi_c from 0.
    if phi_c==0:
        # At phi=0 the central coefficients are exactly their GR values;
        # there is no need to invoke a nonlinear root finder.
        g=c.kappa*(eps_c+3*p_c)/6.
        l=c.kappa*eps_c/3.
        return (np.array([g*r0*r0,2*g*r0,0.,0.,
                          p_c-(eps_c+p_c)*g*r0*r0]),(g,l,0.))
    a,b,k=c.alpha,c.beta,c.kappa
    gg=k*(eps_c+3*p_c)/6  # GR: Gamma_2=kappa*(epsilon_c+3*p_c)/6
    ll=k*eps_c/3
    # The central rr equation gives g=kappa*(epsilon_c+3*p_c)/6.
    def f(v,ph):
        g,l,s=v
        return np.array([96*a*k*l*ph*s+3*b*k*l*ph*ph-12*b*k*ph*s+2*k*eps_c-6*l,
               64*a*k*g*ph*s+2*b*k*g*ph*ph-b*k*l*ph*ph+8*b*k*ph*s-4*g+2*k*p_c+2*l,
               2*s+b*(g-l)*ph-8*a*g*l*ph])
    v=np.array([gg,ll,0.])
    for ph in np.linspace(0,phi_c,max(7,int(abs(phi_c)/.015)+1))[1:]:
        rr=root(lambda x:f(x,ph),v,tol=1e-11)
        if np.max(np.abs(f(rr.x,ph)))>1e-9 or not np.all(np.isfinite(rr.x)) or rr.x[1]<=0:
            raise ValueError(f'No GR-connected regular center root at phi={ph}')
        v=rr.x
    g,l,s=v
    y0=np.array([g*r0**2,2*g*r0,phi_c+s*r0*r0,2*s*r0,p_c-(eps_c+p_c)*g*r0*r0])
    return y0,v

