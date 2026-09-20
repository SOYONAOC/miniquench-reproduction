"""Audit Figure 2 rebounds using saved forces and finite prescribed checks."""
import json
from pathlib import Path
import time
import numpy as np
from miniquench.physics import Config,ROOT
from miniquench.dynamics import run

out={'saved_trajectories':[], 'checks':[]}
for mass in [1e8,1e9]:
 d=np.genfromtxt(ROOT/f'data/curves/fig2_{mass:.0e}_full_sn_only.csv',delimiter=',',names=True)
 v=d['velocity_km_s']; indices=np.flatnonzero((v[:-1]<0)&(v[1:]>=0))+1
 q=d['pressure_erg_cm3']*d['radius_cm']**5*np.exp(d['time_myr']/(120*(7/10)**-4))
 sel=d['time_myr']>80
 out['saved_trajectories'].append({'mass':mass,'bounce_samples':[{k:float(d[k][i]) for k in ['time_myr','r_over_R0','pressure_dyne','gravity_dyne','accretion_dyne','sn_per_myr','sfr_msun_yr']} for i in indices], 'post_sn_Pr5_exp_invariant_relative_span':float(np.ptp(q[sel])/np.mean(q[sel]))})
for mass,branch in [(1e8,'outward_only'),(1e9,'outward_only'),(1e8,'fine_until_return')]:
 cfg=Config(mass=mass,model='sn_only',loading='outward_only' if branch=='outward_only' else 'until_return',history_dt=.025 if branch=='fine_until_return' else .05,max_step=.25 if branch=='fine_until_return' else .5,rtol=2e-9 if branch=='fine_until_return' else 2e-7,atol=1e-11 if branch=='fine_until_return' else 1e-9)
 start=time.perf_counter()
 try:
  result,curve=run(cfg)
  v=curve['velocity_km_s'];result['bounce_count']=int(np.sum((v[:-1]<0)&(v[1:]>=0)))
 except ArithmeticError as exc:
  result={'status':'numerical_failure','error':str(exc),'mass':mass,'wall_seconds':time.perf_counter()-start}
 result['branch']=branch
 out['checks'].append(result)
 print(json.dumps(result),flush=True)
(ROOT/'data/figure2_rebound_audit.json').write_text(json.dumps(out,indent=2)+'\n')
