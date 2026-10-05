"""After building class and the in-place Python extension, run from any directory.
Requires NumPy. Checks data discovery, missing-path errors and CLI/Python agreement.
"""
import sys,json,subprocess,tempfile,os
from pathlib import Path
import numpy as np
repo=Path(__file__).resolve().parents[1];sys.path.insert(0,str(repo/'python'))
from classy import Class,CosmoComputationError
p={'output':'tCl,pCl,lCl,mPk','lensing':'yes','l_max_scalars':800,'P_k_max_1/Mpc':1,'h':.6736,'omega_b':.02237,'omega_cdm':.12,'A_s':2.1e-9,'n_s':.965,'tau_reio':.054}
results={}
# Check ordinary LCDM and accurate full-table lensing against the CLI, retaining default algorithms.
for accurate in [0,1]:
 q={**p,'accurate_lensing':accurate};c=Class();c.set(q);c.compute();cl=c.lensed_cl(800)
 assert all(np.isfinite(cl[k]).all() for k in ['tt','ee','te','bb','pp'])
 with tempfile.TemporaryDirectory() as tmp:
  d=Path(tmp);args={**q,'root':str(d/'out'),'overwrite_root':'yes'};f=d/'input.ini';f.write_text('\n'.join(f'{k} = {v}' for k,v in args.items()))
  r=subprocess.run([str(repo/'class'),str(f)],capture_output=True,text=True,timeout=90,env={**os.environ,'OMP_NUM_THREADS':'1'});assert r.returncode==0,r.stdout+r.stderr
  a=np.loadtxt(d/'out_cl_lensed.dat');ell=a[:,0].astype(int);errs=[]
  for i,k in enumerate(['tt','ee','te','bb','pp'],1):
   value=cl[k][ell]*ell*(ell+1)/(2*np.pi);norm=np.sqrt(a[:,1]*a[:,2]) if k=='te' else abs(a[:,i]);errs.append(float(max(abs(value-a[:,i])/np.maximum(norm,1e-30))))
  assert max(errs)<1e-7;results[f'LCDM_accurate_{accurate}_max_CLI_difference']=max(errs)
 c.struct_cleanup();c.empty()
# Test a small conformal model through the original shooting route.
c=Class();q={**p,'output':'mPk','lensing':'no','omega_cdm':1e-10,'omega_qcdm':.12,'Omega_scf':-1,'Omega_Lambda':0,'Omega_fld':0,'scf_potential':'exp','scf_coupling_type':'conformal','scf_C0':1,'scf_beta':.03,'scf_lambda':.1,'scf_shooting_target':'scf_V0','scf_V0':1,'attractor_ic_scf':'no','scf_phi_ini':1,'scf_phi_prime_ini':-1e-8}
q.pop('l_max_scalars');c.set(q);c.compute();assert np.isfinite(c.pk_lin(.1,0)) and c.pk_lin(.1,0)>0;c.struct_cleanup();c.empty();results['legacy_conformal_shooting']=True
c=Class();c.set({**p,'YHe':0.2454,'base_path':'/tmp/'+'missing_hyrec_data_'*6})
try:c.compute();raise AssertionError('Missing path unexpectedly succeeded')
except CosmoComputationError as e:assert 'could not open file' in str(e);results['long_missing_path_exception']=True
finally:c.struct_cleanup();c.empty()
results['module']=str(sys.modules['classy'].__file__);print(json.dumps(results,indent=2))
