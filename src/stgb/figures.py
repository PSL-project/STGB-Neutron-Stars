"""Seven figure families. No plotting routine writes files or selects a GUI backend."""
import numpy as np
import matplotlib.pyplot as plt
from .eos import MEV_TO_GEO,MSUN_KM
from .stellar import gr_curvature

GR_COLOR='#21252b'
PAL=['#1768ac','#d1495b','#00978d','#ae6a24','#8155a5','#718096']

def style():
    plt.rcParams.update({'figure.dpi':110,'savefig.dpi':220,'axes.labelsize':11,
        'axes.titlesize':12,'axes.linewidth':.85,'font.size':10,'legend.fontsize':9,
        'grid.alpha':.20,'axes.spines.top':False,'axes.spines.right':False,
        'lines.linewidth':2.15,'figure.constrained_layout.use':True})

def segment(ax,rows,x,y,label,**kw):
    """Do NOT connect unverified/disconnected branches across missing density points."""
    xs=np.array([row[x] for row in rows],float);ys=np.array([row[y] for row in rows],float)
    finite=np.isfinite(xs)&np.isfinite(ys)
    inds=np.flatnonzero(finite)
    if not len(inds):return False
    splits=np.split(inds,np.where(np.diff(inds)>1)[0]+1)
    for i,s in enumerate(splits):
        ax.plot(xs[s],ys[s],label=label if i==0 else '_nolegend_',
                marker='o',markersize=3.2,**kw)
    return True

def annotate(ax,txt):
    ax.text(.5,.5,txt,transform=ax.transAxes,ha='center',va='center',
            fontsize=10,color='#7c7c7c',bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})

def mech(study):
    p=study.plan;gr=study.datasets['central']['GR'];st=study.datasets['central']['ST']
    fig,(a,b)=plt.subplots(1,2,figsize=(12,4.5),num='01 — Scalarization mechanism')
    r=np.linspace(study.settings.r0,max(2*gr['radius'],25),500)
    R,GB=gr_curvature(gr,r)
    ric=p.beta*R/2;gb=-p.alpha*GB;total=ric+gb
    a.plot(r,1e3*ric,color=PAL[0],ls='-',label=r'$\beta_V R/2$ (GR geometry)')
    a.plot(r,1e3*gb,color=PAL[2],ls='-.',label=r'$-\alpha\mathcal{G}$ (GR geometry)')
    a.plot(r,1e3*total,color=PAL[1],label=r'$m_{\mathrm{eff}}^2$ (total)')
    a.axhline(0,color=GR_COLOR,ls='--',lw=1);a.axvline(gr['radius'],color='.6',ls=':',lw=1)
    a.set(xlabel='Areal radius [km]',ylabel=r'$10^3\, m_{\rm eff}^2$ [km$^{-2}$]',title='(a) Linear effective mass on the GR background')
    a.legend(fontsize=8)
    b.axhline(0,color=GR_COLOR,ls='--',label='GR: zero scalar field')
    if st:
        r=np.linspace(study.settings.r0,min(study.settings.rmax,max(3*st['radius'],40)),600)
        rin=r<=st['radius'];ph=np.empty_like(r)
        ph[rin]=st['interior'].sol(r[rin])[2];ph[~rin]=st['exterior'].sol(r[~rin])[2]
        b.plot(r,ph,color=PAL[0],label='Scalarized, node 0')
        b.axvline(st['radius'],color='.6',ls=':',lw=1,label='Scalarized surface')
    else:annotate(b,'No verified n=0 branch in the requested search')
    b.set(xlabel='Areal radius [km]',ylabel=r'$\phi_V(r)$',title='(b) Nonlinear scalar profile')
    b.legend(fontsize=8);fig.suptitle(f'{p.eos} | $\\rho_c$={p.rho_c:.2g} g cm$^{{-3}}$, $\\beta_V$={p.beta:g}, $\\alpha$={p.alpha:g} km$^2$')
    return fig

def mr(study):
    p=study.plan;fig,ax=plt.subplots(figsize=(8,5.5),num='02 — Mass-radius relation')
    seq=study.datasets['sequences'];first=next(iter(seq.values()))
    segment(ax,first,'R_GR','M_GR','GR (all $\\alpha$)',color=GR_COLOR,ls='--',lw=2.5)
    got=False
    for i,(alpha,rows) in enumerate(seq.items()):
        got |= segment(ax,rows,'R_ST','M_ST',rf'Scalarized $\alpha={alpha:g}$ km$^2$',color=PAL[i%len(PAL)])
    if not got:annotate(ax,'No verified scalarized points in the selected scan')
    ax.set(xlabel='Areal stellar radius [km]',ylabel=r'ADM mass [$M_\odot$]',title=f'{p.eos} | $\\beta_V={p.beta:g}$ | density sequence')
    ax.legend(loc='best');ax.grid(True);ax.text(.02,.02,'Equilibrium curves; radial stability not assessed',transform=ax.transAxes,fontsize=8,color='.45');return fig

def charge(study):
    p=study.plan;fig,ax=plt.subplots(figsize=(8,5.3),num='03 — Scalar charge')
    rows=next(iter(study.datasets['sequences'].values()))
    masses=np.array([z['M_GR'] for z in rows],float);masses=masses[np.isfinite(masses)]
    if len(masses):ax.plot([masses.min(),masses.max()],[0,0],ls='--',color=GR_COLOR,label='GR: $Q=0$')
    got=False
    for i,(a,seq) in enumerate(study.datasets['sequences'].items()):
        tmp=[dict(M_ST=row['M_ST'],absQ=abs(row['Q_over_M'])) for row in seq]
        got|=segment(ax,tmp,'M_ST','absQ',rf'$\alpha={a:g}$ km$^2$',color=PAL[i%len(PAL)])
    if not got:annotate(ax,'No verified scalarized points in the selected scan')
    ax.set(xlabel=r'ADM mass [$M_\odot$]',ylabel=r'$|Q|/M$',title=f'{p.eos} | $\\beta_V={p.beta:g}$')
    ax.set_ylim(bottom=-.003);ax.legend();ax.grid(True);return fig

def coupling(study):
    p=study.plan;rows=study.datasets['coupling']
    a=np.array([z['alpha'] for z in rows]);d=np.array([z['delta_M_pct'] for z in rows]);q=np.abs([z['Q_over_M'] for z in rows])
    fig,(ax,bx)=plt.subplots(1,2,figsize=(11.4,4.5),num='04 — Coupling scan')
    for x,y,axis,lab in [(a,d,ax,r'$100(M_{\rm STGB}-M_{\rm GR})/M_{\rm GR}$ [%]'),
                         (a,q,bx,r'$|Q|/M_{\rm STGB}$')]:
        good=np.isfinite(y)
        axis.axhline(0,color=GR_COLOR,ls='--',label='GR baseline')
        segment(axis,[dict(x=float(xx),y=float(yy)) for xx,yy in zip(x,y)],'x','y','Verified n=0 solution',color=PAL[1])
        if not np.any(good):annotate(axis,'No verified scalarized points')
        axis.set(xlabel=r'Gauss–Bonnet coupling $\alpha$ [km$^2$]',ylabel=lab)
        axis.grid(True);axis.legend(fontsize=8)
    fig.suptitle(f'{p.eos}: fixed $\\beta_V={p.beta:g}$ and $\\rho_c={p.rho_c:.2g}$ g cm$^{{-3}}$')
    return fig

def eos_compare(study):
    p=study.plan;fig,axes=plt.subplots(2,3,figsize=(13,8),sharex=True,sharey=True,num='05 — EOS comparison')
    for ax,(name,rows) in zip(axes.flat,study.datasets['eos_sequences'].items()):
        segment(ax,rows,'R_GR','M_GR','GR',ls='--',color=GR_COLOR)
        ok=segment(ax,rows,'R_ST','M_ST','Scalarized, n=0',color=PAL[0])
        if not ok:annotate(ax,'No n=0 branch found')
        ax.set_title(name);ax.grid(True)
    for ax in axes[-1]:ax.set_xlabel('Radius [km]')
    for ax in axes[:,0]:ax.set_ylabel(r'Mass [$M_\odot$]')
    axes[0,0].legend(loc='best')
    fig.text(.5,.01,'Equilibrium sequences only; stability not assessed.',ha='center',fontsize=9,color='.45')
    fig.suptitle(rf'Identical couplings for all EOS: $\beta_V={p.beta:g}$, $\alpha={p.alpha:g}$ km$^2$')
    return fig

def internal(study):
    p=study.plan;gr=study.datasets['central']['GR'];st=study.datasets['central']['ST']
    fig,(a,b)=plt.subplots(1,2,figsize=(11.8,4.6),num='06 — Internal structure')
    rgr=np.linspace(study.settings.r0,gr['radius'],350);yg=gr['sol'].sol(rgr)
    epsgr=np.array([gr['eos'].epsilon(z) for z in yg[2]])/MEV_TO_GEO
    a.plot(rgr,yg[2]/MEV_TO_GEO,'--',color=GR_COLOR,label='GR')
    b.plot(rgr,epsgr,'--',color=GR_COLOR,label='GR')
    if st:
        rr=np.linspace(study.settings.r0,st['radius'],350);ys=st['interior'].sol(rr)
        eps=np.array([gr['eos'].epsilon(z) for z in ys[4]])/MEV_TO_GEO
        a.plot(rr,ys[4]/MEV_TO_GEO,color=PAL[0],label='Scalarized')
        b.plot(rr,eps,color=PAL[0],label='Scalarized')
        for ax in (a,b):ax.axvline(st['radius'],ls=':',color=PAL[0],alpha=.6)
    else:
        for ax in (a,b):annotate(ax,'No verified n=0 branch')
    # Overlay dashed GR again after the scalar curve: near-identical solutions remain visually distinguishable.
    a.plot(rgr,yg[2]/MEV_TO_GEO,'--',color=GR_COLOR,lw=1.8,alpha=.95,zorder=7,label='_nolegend_')
    b.plot(rgr,epsgr,'--',color=GR_COLOR,lw=1.8,alpha=.95,zorder=7,label='_nolegend_')
    for ax in (a,b):
        ax.axvline(gr['radius'],color=GR_COLOR,ls=':',alpha=.5);ax.set_xlabel('Areal radius [km]');ax.legend();ax.grid(True)
    a.set(ylabel=r'$p$ [MeV fm$^{-3}$]',title='(a) Pressure profile')
    b.set(ylabel=r'$\varepsilon$ [MeV fm$^{-3}$]',title='(b) Total energy-density profile')
    fig.suptitle(f'{p.eos} | same central energy density | $\\alpha={p.alpha:g}$ km$^2$, $\\beta_V={p.beta:g}$')
    return fig

def onset(study):
    p=study.plan;rows=study.datasets['threshold']
    fig,ax=plt.subplots(figsize=(8,5.5),num='07 — Linear scalarization threshold')
    got=segment(ax,[dict(alpha=d['alpha'],beta=d['beta_critical']) for d in rows],
                'alpha','beta','Fundamental (node-0) zero mode',color=PAL[2])
    if not got:annotate(ax,'No node-0 threshold in searched coupling window')
    ax.axhline(0,color=GR_COLOR,ls='--',lw=1,alpha=.6,label=r'$\beta_V=0$ reference')
    ax.set(xlabel=r'$\alpha$ [km$^2$]',ylabel=r'Critical $\beta_V$',
        title=f'{p.eos} | $\\rho_c={p.rho_c:.2g}$ g cm$^{{-3}}$')
    ax.text(.02,.03,'Linear GR-background threshold; not a nonlinear existence or stability map',
        transform=ax.transAxes,fontsize=8,color='.4')
    ax.grid(True);ax.legend();return fig

BUILD={1:mech,2:mr,3:charge,4:coupling,5:eos_compare,6:internal,7:onset}
NAMES={1:'01_effective_mass_and_scalar',2:'02_mass_radius',3:'03_scalar_charge',4:'04_alpha_dependence',
       5:'05_six_eos_mass_radius',6:'06_interior_structure',7:'07_linear_zero_mode_threshold'}

def build(study,figures):
    style();return {n:BUILD[n](study) for n in figures}
