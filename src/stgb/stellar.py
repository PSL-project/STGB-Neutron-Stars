"""GR backgrounds, nonlinear scalarized stars, branch shooting and diagnostics."""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq,root
from .gravity import Couplings,rhs,algebraic_x,regular_center
from .eos import MEV_TO_GEO,MSUN_KM

@dataclass
class Settings:
    r0:float=1e-3
    rmax:float=400.
    max_step:float=.8
    rtol:float=3e-8
    shoot_low:float=1e-5
    shoot_high:float=.12
    shoot_points:int=14
    verbose:bool=False

def solve_gr(eos,rho_c,st:Settings):
    ec,pc=eos.central(rho_c)
    def ode(r,y):
        m,nu,p=y;ep=eos.epsilon(p)
        if r<=2*m:raise ValueError('GR metric singular')
        upr=(m+4*np.pi*r**3*p)/(r*(r-2*m))
        return [4*np.pi*r*r*ep,upr,-(ep+p)*upr]
    def surface(r,y):return y[2]-eos.p_min
    surface.terminal=True;surface.direction=-1
    sol=solve_ivp(ode,(st.r0,65.),[4*np.pi*ec*st.r0**3/3,0.,pc],
                  dense_output=True,events=surface,rtol=st.rtol,
                  atol=[1e-11,1e-11,1e-19],max_step=min(st.max_step,.2))
    if not sol.success or len(sol.t_events[0])==0:raise RuntimeError('GR stellar surface not found: '+sol.message)
    radius=float(sol.t_events[0][0]);mass=float(sol.y_events[0][0][0]);
    return dict(eos=eos,rho_c=rho_c,eps_c=ec,p_c=pc,sol=sol,radius=radius,mass=mass,
                M_solar=mass/MSUN_KM)

def solve_scalar(eos,rho_c,phi_c,c:Couplings,st:Settings):
    ec,pc=eos.central(rho_c)
    y0,_=regular_center(ec,pc,float(phi_c),c,st.r0)
    def surface(r,y):return y[4]-eos.p_min
    surface.terminal=True;surface.direction=-1
    sol=solve_ivp(lambda r,y:rhs(r,y,eos.epsilon,c),(st.r0,min(st.rmax,65.)),y0,
                  events=surface,dense_output=True,rtol=st.rtol,
                  atol=[1e-11,1e-11,1e-11,1e-11,1e-19],max_step=min(.16,st.max_step))
    if not sol.success or not len(sol.t_events[0]):raise RuntimeError('Scalar interior failed: '+sol.message)
    rs=float(sol.t_events[0][0]);ys=sol.y_events[0][0].copy();ys[4]=0.
    vac=solve_ivp(lambda r,y:rhs(r,y,lambda p:0.,c),(rs,st.rmax),ys,
                  dense_output=True,rtol=st.rtol,
                  atol=[1e-11,1e-11,1e-11,1e-11,1e-11],max_step=st.max_step)
    if not vac.success:raise RuntimeError('Scalar exterior failed: '+vac.message)
    y=vac.y[:,-1];x=algebraic_x(st.rmax,y[1],y[2],y[3],0.,c)
    m=st.rmax*(1-1/x)/2
    q=-st.rmax**2*y[3]/(1+2*m/st.rmax)
    phi_inf=float(y[2]+st.rmax*y[3]+m*q/st.rmax**2)
    return dict(eos=eos,rho_c=rho_c,phi_c=phi_c,interior=sol,exterior=vac,
                radius=rs,phi_inf=phi_inf,settings=st,params=c)

def observables(sol):
    c=sol['params'];so=sol['exterior'];r=np.linspace(.65*so.t[-1],so.t[-1],100)
    y=so.sol(r)
    H=np.stack([np.ones_like(r),1/r,1/r**2],axis=1)
    phiinf,q,q2=np.linalg.lstsq(H,y[2],rcond=None)[0]
    mm=np.array([rr*(1.-1./algebraic_x(rr,yy[1],yy[2],yy[3],0.,c))/2
                 for rr,yy in zip(r,y.T)])
    M=np.linalg.lstsq(np.stack([np.ones_like(r),1/r],axis=1),mm,rcond=None)[0][0]
    if M<=0 or not np.isfinite(M):raise ValueError('Nonpositive or invalid ADM mass fit')
    return dict(M_km=float(M),M_solar=float(M/MSUN_KM),Q_km=float(q),
        Q_over_M=float(q/M),phi_inf_fit=float(phiinf),phi_inf_est=float(sol['phi_inf']))

def shoot(eos,rho_c,c,st:Settings):
    """Scan positive phi_c, choose verified n=0 root, never turn a failed scan into GR."""
    grid=np.geomspace(st.shoot_low,st.shoot_high,st.shoot_points)
    def mismatch(z):return solve_scalar(eos,rho_c,z,c,st)['phi_inf']
    samples=[]
    for z in grid:
        try: samples.append(mismatch(float(z)))
        except (ValueError,RuntimeError,OverflowError,np.linalg.LinAlgError):samples.append(np.nan)
    candidates=[]
    for j in range(len(grid)-1):
        if not np.isfinite(samples[j]) or not np.isfinite(samples[j+1]) or samples[j]*samples[j+1]>=0:continue
        try:
            phi_c=brentq(mismatch,grid[j],grid[j+1],xtol=2e-10,rtol=1e-8,maxiter=60)
            s=solve_scalar(eos,rho_c,phi_c,c,st)
            ri=np.linspace(st.r0,s['radius'],280)
            re=np.linspace(s['radius'],st.rmax,280)
            vals=np.concatenate([s['interior'].sol(ri)[2],s['exterior'].sol(re)[2]])
            nodes=int(np.count_nonzero(vals[1:]*vals[:-1]<0))
            s['nodes']=nodes;s['fit']=observables(s)
            # Reject non-asymptotically flat branches and spuriously bad shooting roots.
            if abs(s['fit']['phi_inf_fit'])>max(2e-4,0.02*abs(phi_c)):continue
            if nodes==0:candidates.append(s)
        except (ValueError,RuntimeError,np.linalg.LinAlgError):continue
    return (min(candidates,key=lambda s:abs(s['phi_c'])) if candidates else None,
            {'phi_c':grid,'phi_infinity':np.array(samples)})

def gr_curvature(background,r):
    """GR R, G and scalarization terms without differentiating noisy EOS arrays.

    R=kappa*(epsilon-3p)
    GB=48 m^2/r^6 -16 kappa m epsilon/r^3 -4 kappa^2 epsilon p
    """
    r=np.asarray(r,dtype=float)
    rr=np.maximum(r,background['sol'].t[0]);yy=background['sol'].sol(np.minimum(rr,background['radius']))
    m=yy[0];p=np.maximum(yy[2],0)
    eps=np.array([background['eos'].epsilon(z) for z in p])
    outside=rr>background['radius'];m=np.where(outside,background['mass'],m)
    eps=np.where(outside,0,eps);p=np.where(outside,0,p)
    k=8*np.pi
    Ric=k*(eps-3*p)
    GB=48*(m/r**3)**2-16*k*(m/r**3)*eps-4*k*k*eps*p
    return Ric,GB
