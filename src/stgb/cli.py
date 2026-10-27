"""Interactive research CLI: preview first, save only after preview/approval."""
import argparse
import csv
import json
import os
import sys
import tempfile
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from .eos import load_eos,NAMES,MEV_TO_GEO
from .analysis import Study,Plan
from .stellar import Settings,gr_curvature
from .figures import build,NAMES as FIG_NAMES

REPO=Path(__file__).resolve().parents[2]

def parse(argv=None):
    pa=argparse.ArgumentParser(description='Ricci–Gauss–Bonnet neutron stars: seven research figures')
    pa.add_argument('--eos',choices=NAMES,default='SLy4')
    pa.add_argument('--rho-c',type=float,default=8.1e14,help='epsilon_c/c^2 in g/cm^3 (NOT baryon rest-mass density)')
    pa.add_argument('--alpha',type=float,default=30.,help='Gauss–Bonnet coupling [km^2]')
    pa.add_argument('--beta',type=float,default=-10.,help='Ricci coupling in canonical scalar convention')
    pa.add_argument('--quality',choices=['quick','thesis'],default='quick')
    pa.add_argument('--figures',default='all',help='all or comma-separated integers, e.g. 1,4,6')
    pa.add_argument('--data-dir',type=Path,default=REPO/'data')
    pa.add_argument('--output',type=Path,default=Path.home()/'Documents'/'STGB_Results')
    pa.add_argument('--no-show',action='store_true',help='Headless mode (for servers/CI)')
    pa.add_argument('--save',action='store_true',help='Save plots and data after preview; with --no-show, save immediately')
    pa.add_argument('--preview-only',action='store_true',help='Display and never save')
    pa.add_argument('--verbose',action='store_true')
    pa.add_argument('--rmax',type=float,default=400.,help='Outer scalar integration radius [km]')
    pa.add_argument('--rtol',type=float,default=3e-8)
    pa.add_argument('--shoot-high',type=float,default=.12)
    pa.add_argument('--sequence-alphas',default=None,help='Comma-separated alpha [km^2] for plots 2 and 3')
    pa.add_argument('--coupling-alphas',default=None,help='Comma-separated alpha scan for plot 4')
    pa.add_argument('--threshold-alphas',default=None,help='Comma-separated alpha scan for plot 7')
    pa.add_argument('--density-low',type=float,default=None,help='Minimum central epsilon/c^2 [g/cm3]')
    pa.add_argument('--density-high',type=float,default=None,help='Maximum central epsilon/c^2 [g/cm3]')
    pa.add_argument('--sequence-points',type=int,default=None,help='Number of central densities in M-R sequences')
    pa.add_argument('--beta-min',type=float,default=None,help='Min beta for zero-mode root search')
    pa.add_argument('--beta-max',type=float,default=None,help='Max beta for zero-mode root search')

    a=pa.parse_args(argv)
    a.figures=list(range(1,8)) if a.figures.strip().lower()=='all' else sorted(set(int(x) for x in a.figures.split(',')))
    if not a.figures or any(x not in range(1,8) for x in a.figures):pa.error('--figures must be all or 1..7')
    if a.preview_only and a.save:pa.error('Choose either --preview-only or --save')
    if a.no_show and not a.save:pa.error('Headless mode requires --save (or remove --no-show)')
    return a

def csv_rows(path,rows):
    if not rows:return
    with path.open('w',newline='',encoding='utf8') as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0].keys()));w.writeheader()
        for row in rows:w.writerow(row)

def rows_for_save(study):
    p=study.plan;d=study.datasets;tables={}
    if 'sequences' in d:
        for alpha,rows in d['sequences'].items():tables[f'sequence_{p.eos}_alpha_{alpha:+g}.csv']=rows
    if 'eos_sequences' in d:
        for name,rows in d['eos_sequences'].items():tables[f'eos_compare_{name}_alpha_{p.alpha:+g}.csv']=rows
    if 'coupling' in d:tables['coupling_scan.csv']=d['coupling']
    if 'threshold' in d:tables['linear_zero_mode_threshold.csv']=d['threshold']
    if 'central' in d:
        g=d['central']['GR'];s=d['central']['ST']
        r=np.linspace(study.settings.r0,max(g['radius']*2,25),450)
        Ric,GB=gr_curvature(g,r)
        be=p.beta*Ric/2;ga=-p.alpha*GB
        tables['effective_mass_gr_background.csv']=[dict(r_km=float(x),R_km_minus2=float(v),GB_km_minus4=float(w),
            ricci_term_km_minus2=float(z),gb_term_km_minus2=float(t),
            meff_sq_km_minus2=float(z+t)) for x,v,w,z,t in zip(r,Ric,GB,be,ga)]
        rr=np.linspace(study.settings.r0,g['radius'],350);yy=g['sol'].sol(rr)
        tables['selected_gr_profile.csv']=[dict(r_km=float(x),m_Msun=float(m)/1.4766250385,
             p_MeV_fm3=float(pp)/MEV_TO_GEO,epsilon_MeV_fm3=float(g['eos'].epsilon(pp))/MEV_TO_GEO)
             for x,m,pp in zip(rr,yy[0],yy[2])]
        if s is not None:
            rr=np.linspace(study.settings.r0,min(study.settings.rmax,max(s['radius']*3,40)),400)
            data=[]
            for r in rr:
                v=s['interior'].sol(r) if r<=s['radius'] else s['exterior'].sol(r)
                data.append(dict(r_km=float(r),phi=float(v[2]),phi_prime=float(v[3]),Gamma=float(v[0]),
                        p_MeV_fm3=float(max(0,v[4])/MEV_TO_GEO)
                        if r<=s['radius'] else 0.))
            tables['selected_scalar_profile.csv']=data
    return tables

def save_all(dest,plots,study,params):
    """Try explicit destination, then a safe user-home fallback; never write next to Python script."""
    fallback=Path(tempfile.gettempdir())/'STGB_Results'
    targets=[Path(dest).expanduser(),Path.home()/'STGB_Results',fallback]
    failures=[]
    out=None
    for candidate in targets:
        try:
            candidate.mkdir(parents=True,exist_ok=True)
            marker=candidate/'.write_test'
            marker.write_text('ok',encoding='utf-8');marker.unlink()
            out=candidate;break
        except (OSError,PermissionError) as ex:failures.append(f'{candidate}: {ex}')
    if out is None:raise OSError('No writable output folder. '+ '; '.join(failures))
    if out.resolve()!=Path(dest).expanduser().resolve():
        print('Requested results folder is not writable; using:',out)
        print('Original errors:',failures)
    run_folder=out/('STGB_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
    run_folder.mkdir(parents=True,exist_ok=True)
    figures_dir=run_folder/'figures';tables_dir=run_folder/'data'
    figures_dir.mkdir();tables_dir.mkdir()
    for n,fig in plots.items():
        for ext in ('png','pdf'):
            fig.savefig(figures_dir/(FIG_NAMES[n]+'.'+ext),dpi=240,bbox_inches='tight')
    for name,rows in rows_for_save(study).items():csv_rows(tables_dir/name,rows)
    meta=dict(created_utc=datetime.now(timezone.utc).isoformat(),settings=vars(params).copy(),
            plan=vars(study.plan),solver_settings=vars(study.settings),
            eos_quality={n:dict(removed_rows=e.filtered_rows,**e.causality_summary()) for n,e in study.eos.items()},
            failures=study.errors,
            limitations=['Scalar branch absence means not detected by limited shooting scan, NOT proven nonexistence',
                'Surface is minimum EOS table pressure, not mathematically p=0',
                'Onset plot is linear zero-mode threshold, not nonlinear solution domain',
                'No full pointwise publication benchmark reproduction yet',
                'Confirm data provenance and permissions before distributing EOS tables publicly'])
    meta['settings']['data_dir']=str(params.data_dir);meta['settings']['output']=str(params.output)
    (run_folder/'run_metadata.json').write_text(json.dumps(meta,indent=2,default=float),encoding='utf-8')
    print('\nSaved',len(plots),'figures (PNG + PDF) and CSVs to:',run_folder.resolve())
    return run_folder

def main(argv=None):
    args=parse(argv)
    if args.no_show:
        import matplotlib
        matplotlib.use('Agg')
    else:
        import matplotlib
        if 'agg' in matplotlib.get_backend().lower():
            print('WARNING: non-interactive Matplotlib backend detected. Use Windows Python/TkAgg for figure windows.')
    settings=Settings(rmax=args.rmax,rtol=args.rtol,shoot_high=args.shoot_high,verbose=args.verbose)
    plan=Plan(eos=args.eos,rho_c=args.rho_c,alpha=args.alpha,beta=args.beta)
    if args.quality=='quick':
        plan.sequence_points=5
        plan.sequence_alphas=(-100.,30.,200.)
        plan.coupling_alphas=(-200.,0.,30.,200.)
        plan.threshold_alpha=(-120.,0.,120.,225.)
        plan.threshold_samples=14
        settings.shoot_points=11
    else:
        plan.sequence_points=17
        plan.sequence_alphas=(-200.,0.,30.,100.,200.)
        plan.coupling_alphas=tuple(np.linspace(-250,250,21))
        plan.threshold_alpha=tuple(np.linspace(-200,250,19))
        plan.threshold_samples=24
        settings.shoot_points=24
    def alpha_values(raw):
        return tuple(float(x) for x in raw.split(',') if x.strip())
    if args.sequence_alphas:plan.sequence_alphas=alpha_values(args.sequence_alphas)
    else:plan.sequence_alphas=tuple(sorted(set((*plan.sequence_alphas,float(args.alpha)))))
    if args.coupling_alphas:plan.coupling_alphas=alpha_values(args.coupling_alphas)
    if args.threshold_alphas:plan.threshold_alpha=alpha_values(args.threshold_alphas)
    if args.density_low is not None:plan.rho_low=args.density_low
    if args.density_high is not None:plan.rho_high=args.density_high
    if args.sequence_points is not None:plan.sequence_points=args.sequence_points
    if args.beta_min is not None:plan.threshold_beta_min=args.beta_min
    if args.beta_max is not None:plan.threshold_beta_max=args.beta_max
    if plan.sequence_points<2 or not (0<plan.rho_low<plan.rho_high):raise ValueError('Invalid central-density sequence range')
    if plan.threshold_beta_min>=plan.threshold_beta_max:raise ValueError('Invalid beta threshold range')

    names=list(NAMES) if 5 in args.figures else [args.eos]
    eos={n:load_eos(n,args.data_dir) for n in names}
    print('\nSTGB research package | model: Ricci + Gauss–Bonnet, canonical scalar')
    print('EOS:',names,'| figure IDs:',args.figures,'|',args.quality)
    print('GR baseline dashed; scalarized n=0 solid. No missing-point interpolation.')
    study=Study(eos,plan,settings)
    study.calculate(args.figures)
    plots=build(study,args.figures)
    print('\nFigure calculations complete. Figures are prepared in memory.')
    if study.errors:
        print(f'Note: {len(study.errors)} failed/nonregular configurations. See metadata if saved.')
    if not args.no_show:
        print('Previewing figures now; CLOSE the plot windows when finished examining them.')
        import matplotlib.pyplot as plt
        plt.show(block=True)
    should_save=args.save
    if not args.save and not args.preview_only and not args.no_show:
        if sys.stdin.isatty():
            try:should_save=input('\nSave the reviewed figures and numerical CSVs? [y/N]: ').strip().lower() in ('y','yes')
            except EOFError:should_save=False
        else:print('No terminal input; nothing saved. To save next time, use --save.')
    if should_save:save_all(args.output,plots,study,args)
    else:print('No files saved. Use --save to request saving after preview.')
    if study.errors:
        print('First calculation issues:',*study.errors[:8],sep='\n - ')
    return 0
