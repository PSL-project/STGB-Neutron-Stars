"""GR-background zero-frequency scalar eigenmodes (bifurcation thresholds).

These linear zero modes are not maps of nonlinear solution existence or stability.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from .stellar import gr_curvature

def zero_mode(gr,alpha,beta,rmax=250.,rtol=4e-7):
    rs=gr['radius'];mstar=gr['mass'];r0=float(gr['sol'].t[0]);epsc=gr['eps_c'];pc=gr['p_c']
    k=8*np.pi
    gc=-(k*k)*(4*epsc*epsc/3+4*epsc*pc)
    meff0=beta*k*(epsc-3*pc)/2-alpha*gc
    def radial(r,y,isvac):
        phi,psi=y
        if isvac:
            m=mstar;p=eps=0.
        else:
            m,nu,p=gr['sol'].sol(r)
            eps=gr['eos'].epsilon(max(p,gr['eos'].p_min))
        A=1-2*m/r
        if A<=0:raise ValueError('Trapped surface in zero-mode background')
        mprime=4*np.pi*r*r*eps
        nuprime=(m+4*np.pi*r**3*p)/(r*(r-2*m))
        GammaPrime=2*nuprime
        LambdaPrime=(2*mprime/r-2*m/r**2)/A
        Ric=k*(eps-3*p)
        GB=48*m*m/r**6-16*k*m*eps/r**3-4*k*k*eps*p
        pprime=-(2/r+(GammaPrime-LambdaPrime)/2)*psi+(beta*Ric/2-alpha*GB)*phi/A
        return [psi,pprime]
    y0=[1+meff0*r0*r0/6, meff0*r0/3]
    si=solve_ivp(lambda r,y:radial(r,y,False),(r0,rs),y0,rtol=rtol,atol=1e-9,max_step=.35,dense_output=True)
    if not si.success:raise RuntimeError('zero mode interior failed')
    so=solve_ivp(lambda r,y:radial(r,y,True),(rs,rmax),si.y[:,-1],rtol=rtol,atol=1e-9,max_step=1.,dense_output=True)
    if not so.success:raise RuntimeError('zero mode exterior failed')
    ph,ps=so.y[:,-1]
    # Schwarzschild static massless scalar asymptotes: phi_inf = phi + small corrections.
    # subtract asymptotically decaying part with its exact Schwarzschild radial profile
    # phi(r)=phi_inf + C ln(1-2M/r), psi=C*2M/(r(r-2M)).
    C=ps*rmax*(rmax-2*mstar)/(2*mstar)
    phi_inf=ph-C*np.log1p(-2*mstar/rmax)
    vals=np.r_[si.y[0],so.y[0]]
    nodes=int(np.count_nonzero(vals[:-1]*vals[1:]<0))
    return float(phi_inf),nodes

def find_critical_beta(gr,alpha,beta_min=-40.,beta_max=15.,samples=19):
    """Find fundamental n=0 zero-mode boundary, else NaN; no interpolation of missing roots."""
    grid=np.linspace(beta_min,beta_max,samples)
    vals=[]
    for b in grid:
        try: vals.append(zero_mode(gr,alpha,b)[0])
        except (ValueError,RuntimeError,OverflowError):vals.append(np.nan)
    roots=[]
    for i in range(len(grid)-1):
        if not np.isfinite(vals[i]) or not np.isfinite(vals[i+1]) or vals[i]*vals[i+1]>=0:continue
        try:
            beta=brentq(lambda b:zero_mode(gr,alpha,b)[0],grid[i],grid[i+1],xtol=1e-4)
            mismatch,nodes=zero_mode(gr,alpha,beta)
            # At exact zero mode, 1/r tail decays and no nodes => fundamental threshold.
            if nodes==0:roots.append(beta)
        except (RuntimeError,ValueError):pass
    return float(roots[0]) if roots else np.nan
