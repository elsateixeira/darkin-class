"""Run with: python test/test_ide_numerics.py --executable ./class (NumPy required)."""
import argparse,json,os,re,subprocess,tempfile
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser();parser.add_argument('--executable',default='./class');args=parser.parse_args();exe=Path(args.executable).resolve()
fixtures=json.loads((Path(__file__).parent/'ide_numerics/inputs.json').read_text())
with tempfile.TemporaryDirectory(prefix='ide-numerics-') as tmp:
 def run(name,p,expect_success=True):
  d=Path(tmp)/name;d.mkdir();p={**p,'root':str(d/'out'),'overwrite_root':'yes','headers':'yes'}
  ini=d/'input.ini';ini.write_text('\n'.join(f'{k} = {v}' for k,v in p.items())+'\n')
  r=subprocess.run([str(exe),str(ini)],capture_output=True,text=True,timeout=180,env={**os.environ,'OMP_NUM_THREADS':'1'})
  if expect_success:assert r.returncode==0,r.stdout+r.stderr
  return d,r
 for model in ['momentum','conformal','continuation']:
  final_fields=[]
  for seed in ([1e-6,1.,1e6] if model!='continuation' else [1.]):
   p={**fixtures[model],'scf_V0':seed};d,r=run(f'{model}_{seed}',p)
   f=d/'out_background.dat';a=np.loadtxt(f);assert np.isfinite(a).all()
   header=[s for s in f.read_text().splitlines() if s.startswith('#')][-1]
   c={v:int(k)-1 for k,v in re.findall(r'(\d+):(\S+)',header)};h=float(p['h']);h0=h/2997.92458
   assert abs((a[-1,c['H']]/h0)**2-1)<5e-8
   assert abs(a[-1,c['(.)rho_qcdm']]/h0**2-float(p['omega_qcdm'])/h**2)<3e-8
   final_fields.append(a[-1,c['phi_scf']])
  assert max(final_fields)-min(final_fields)<1e-5
  print(model,'shooting/closure passed',flush=True)
 spectra=[]
 for chunk in [0,64,256]:
  p={**fixtures['cmb'],'lensing_mu_chunk_size':chunk,'l_max_scalars':1000,'delta_l_max':500,'accurate_lensing':1}
  d,r=run(f'chunk_{chunk}',p);a=np.loadtxt(d/'out_cl_lensed.dat');assert np.isfinite(a).all();spectra.append(a)
 for a in spectra[1:]:
  b=spectra[0];assert np.max(abs(a[:,1:]-b[:,1:])/np.maximum(abs(b[:,1:]),1e-30))<1e-7
 print('lensing chunk equivalence passed',flush=True)
 p={**fixtures['cmb'],'scf_gamma0':.3,'scf_lambda':1.}
 d,r=run('pole',p,False);assert r.returncode!=0 and 'Momentum Euler denominator' in r.stdout+r.stderr
 print('momentum pole rejection passed',flush=True)
