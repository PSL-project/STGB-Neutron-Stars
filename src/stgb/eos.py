"""EOS ingestion, dimensional checks, and shape-preserving interpolation."""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

MEV_TO_GEO = 1.602176634e33*6.67430e-8/(2.99792458e10**4)*1e10  # MeV/fm³ -> km^-2
MSUN_KM = 1.4766250385
RHO_CONVERT = 1.602176634e33/(2.99792458e10**2)  # MeV/fm³ -> g/cm³
NAMES = ('APR','SLy4','DD2','FSUGarnet','G3','IOPB')

class EOS:
    """Tabulated epsilon(P), preserving monotone log-log segments via PCHIP.

    rho_g_cm-3 in supplied tables is epsilon/c², NOT baryonic density.
    Surface is a finite-pressure cutoff at the lowest tabulated point.
    We raise on nonmonotonic pairs instead of silently removing any EOS point.
    """
    def __init__(self,path,name=None):
        source=pd.read_csv(path)
        if 'EOS' in source.columns:
            if name is None: raise ValueError('Specify an EOS name for combined table')
            source=source.loc[source['EOS'].astype(str).str.strip()==name]
        if len(source)<15: raise ValueError('Too few EOS rows')
        self.name=name or Path(path).stem.split('_')[0]
        fields=['rho_g_cm-3','epsilon_MeV_fm-3','P_MeV_fm-3']
        data=source[fields].to_numpy(dtype=float)
        if not np.all(np.isfinite(data)) or np.any(data<=0):raise ValueError(f'{self.name}: nonpositive/nonfinite EOS values')
        data=data[np.argsort(data[:,1])]
        if np.any(np.diff(data[:,1])<=0):raise ValueError(f'{self.name}: epsilon not strictly increasing')
        # Small supplied crust-table reversals are not invertible as epsilon(P).
        # Filter the offending points, with a recorded WARNING rather than hiding changes.
        bad=int(np.count_nonzero(np.diff(data[:,2])<=0))
        if bad:
            retained=[]; last=-float('inf')
            for row in data:
                if row[2]>last:
                    retained.append(row);last=row[2]
            removed=len(data)-len(retained)
            warnings.warn(f'{self.name}: {bad} local pressure reversals; excluded {removed} points to obtain invertible epsilon(P). Inspect original crust EOS.',RuntimeWarning)
            data=np.asarray(retained)
        self.filtered_rows=int(len(source)-len(data))
        rho,eps,prs=data.T
        ratio=rho/(eps*RHO_CONVERT)
        if np.max(abs(ratio-1))>0.01:
            raise ValueError(f'{self.name}: supplied rho is not epsilon/c² to 1%; check metadata')
        self.rho=rho;self.energy_mev=eps;self.pressure_mev=prs
        self.energy=eps*MEV_TO_GEO;self.pressure=prs*MEV_TO_GEO
        self.p_min,self.p_max=float(self.pressure[0]),float(self.pressure[-1])
        self.eps_min,self.eps_max=float(self.energy[0]),float(self.energy[-1])
        self.ep=PchipInterpolator(np.log(self.pressure),np.log(self.energy),extrapolate=False)
        self.pe=PchipInterpolator(np.log(self.energy),np.log(self.pressure),extrapolate=False)
        self.pr=PchipInterpolator(np.log(rho),np.log(self.pressure),extrapolate=False)
    def epsilon(self,p):
        if p<=self.p_min: return self.eps_min
        if p>self.p_max*(1+1e-11):raise ValueError(f'{self.name}: EOS pressure out of bounds: {p}')
        return float(np.exp(self.ep(np.log(min(p,self.p_max)))))
    def central(self,rho_c):
        if not (self.rho[0]<rho_c<self.rho[-1]):raise ValueError(f'{self.name}: central density outside EOS table')
        p=float(np.exp(self.pr(np.log(rho_c))))
        return self.epsilon(p),p
    def causality_summary(self):
        # Approximate sound-speed squared via finite secants; diagnostic not a physical EOS audit.
        cs2=np.diff(self.pressure)/np.diff(self.energy)
        return dict(max_cs2=float(np.max(cs2)),fraction_superluminal=float(np.mean(cs2>1)))


def load_eos(name,data_dir):
    path=Path(data_dir)/f'{name}_standardized.csv'
    if not path.exists():raise FileNotFoundError(f'EOS missing: {path}')
    return EOS(path,name)
