"""Compute datasets behind all 7 figure families. Missing physics stays NaN."""
from dataclasses import dataclass
import numpy as np
from .gravity import Couplings
from .stellar import Settings,solve_gr,shoot,gr_curvature
from .onset import find_critical_beta
from .eos import NAMES

@dataclass
class Plan:
    eos:str='SLy4'
    rho_c:float=8.1e14
    beta:float=-10.
    alpha:float=30.
    sequence_alphas:tuple=(-100.,0.,100.,200.)
    coupling_alphas:tuple=(-200.,-100.,0.,100.,200.)
    rho_low:float=4e14
    rho_high:float=1.8e15
    sequence_points:int=8
    threshold_alpha:tuple=(-150.,-75.,0.,75.,150.,225.)
    threshold_beta_min:float=-35.
    threshold_beta_max:float=5.
    threshold_samples:int=16

class Study:
    def __init__(self,all_eos,plan,settings,log=print):
        self.eos=all_eos;self.plan=plan;self.settings=settings;self.log=log
        self.gr_cache={};self.scalar_cache={};self.scan_cache={}
        self.datasets={};self.errors=[]
    def GR(self,name,rho):
        key=(name,float(rho))
        if key not in self.gr_cache:
            self.gr_cache[key]=solve_gr(self.eos[name],rho,self.settings)
        return self.gr_cache[key]
    def ST(self,name,rho,alpha=None,beta=None):
        a=self.plan.alpha if alpha is None else alpha
        b=self.plan.beta if beta is None else beta
        key=(name,float(rho),float(a),float(b))
        if key not in self.scalar_cache:
            try:
                scalar,trace=shoot(self.eos[name],rho,Couplings(a,b),self.settings)
            except (ValueError,RuntimeError,OverflowError,np.linalg.LinAlgError) as exc:
                scalar=None;trace=None;self.errors.append(f'{key}: {exc}')
            self.scalar_cache[key]=scalar;self.scan_cache[key]=trace
            if self.settings.verbose:self.log(f'{name} rho={rho:.3g}, alpha={a:g} beta={b:g}: '+ ('n=0 found' if scalar else 'no verified n=0 branch'))
        return self.scalar_cache[key]
    def sequence(self,name,alpha):
        p=self.plan
        densities=np.geomspace(p.rho_low,p.rho_high,p.sequence_points)
        rows=[]
        for j,rho in enumerate(densities):
            try:
                gr=self.GR(name,float(rho))
                st=self.ST(name,float(rho),alpha)
                o=st['fit'] if st else {}
                rows.append(dict(rho_c=float(rho),M_GR=float(gr['M_solar']),R_GR=float(gr['radius']),
                    M_ST=o.get('M_solar',np.nan),R_ST=float(st['radius']) if st else np.nan,
                    Q=o.get('Q_km',np.nan),Q_over_M=o.get('Q_over_M',np.nan),
                    phi_c=float(st['phi_c']) if st else np.nan, nodes=float(st['nodes']) if st else np.nan))
            except (ValueError,RuntimeError) as exc:
                self.errors.append(f'{name} rho={rho:.4g}: {exc}')
                rows.append(dict(rho_c=float(rho),M_GR=np.nan,R_GR=np.nan,M_ST=np.nan,R_ST=np.nan,Q=np.nan,Q_over_M=np.nan,phi_c=np.nan,nodes=np.nan))
            self.log(f'  {name}: alpha={alpha:g}, density {j+1}/{len(densities)} -> GR='+
                     (f"{rows[-1]['M_GR']:.3f} M_sun" if np.isfinite(rows[-1]['M_GR']) else 'failed')+
                     (f", scalar={rows[-1]['M_ST']:.3f} M_sun" if np.isfinite(rows[-1]['M_ST']) else ', scalar not verified'))
        return rows
    def calculate(self,figures):
        p=self.plan;fs=set(figures)
        if fs & {1,6}:
            self.log('Figure 1/6: central configuration')
            g=self.GR(p.eos,p.rho_c)
            s=self.ST(p.eos,p.rho_c)
            self.datasets['central']=dict(GR=g,ST=s)
        if fs & {2,3}:
            self.log('Figures 2/3: mass-radius and charge sequences')
            self.datasets['sequences']={alpha:self.sequence(p.eos,alpha) for alpha in p.sequence_alphas}
        if 4 in fs:
            self.log('Figure 4: coupling parameter scan')
            g=self.GR(p.eos,p.rho_c);data=[]
            for a in p.coupling_alphas:
                self.log(f'  scanning alpha={a:g} km^2')
                s=self.ST(p.eos,p.rho_c,a)
                o=s['fit'] if s else {}
                data.append(dict(alpha=a,delta_M_pct=100*(o['M_solar']/g['M_solar']-1) if s else np.nan,
                    Q_over_M=o.get('Q_over_M',np.nan),M_ST=o.get('M_solar',np.nan),M_GR=g['M_solar']))
            self.datasets['coupling']=data
        if 5 in fs:
            self.log('Figure 5: six EOS comparisons')
            self.datasets['eos_sequences']={name:self.sequence(name,p.alpha) for name in self.eos}
        if 7 in fs:
            self.log('Figure 7: true linear zero-mode thresholds (not nonlinear existence domains)')
            g=self.GR(p.eos,p.rho_c);data=[]
            for a in p.threshold_alpha:
                try:
                    bc=find_critical_beta(g,a,p.threshold_beta_min,p.threshold_beta_max,p.threshold_samples)
                except (ValueError,RuntimeError) as exc:
                    bc=np.nan;self.errors.append(f'onset alpha={a}: {exc}')
                data.append(dict(alpha=float(a),beta_critical=float(bc)))
                self.log(f'  onset alpha={a:.0f}: beta_critical={bc:.5g}')
            self.datasets['threshold']=data
        return self.datasets
