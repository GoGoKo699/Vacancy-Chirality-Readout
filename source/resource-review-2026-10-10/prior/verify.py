#!/usr/bin/env python3
"""Run new and preserved scientific checks into a new external directory.

Canonical byte differences and scientific assertions are separate receipt facts.
The archived sources and tolerances are never changed to obtain a passing run.
"""
from __future__ import annotations
import argparse,hashlib,json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def inventory():return {p.relative_to(ROOT).as_posix():sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file()}
def diff(a,b,path=''):
    out=[]
    if isinstance(a,dict) and isinstance(b,dict):
        if a.keys()!=b.keys():out.append({'path':path,'type':'keys'})
        for k in sorted(a.keys()&b.keys()):out+=diff(a[k],b[k],path+'/'+k)
    elif isinstance(a,list) and isinstance(b,list):
        if len(a)!=len(b):out.append({'path':path,'type':'length'})
        for i,(x,y) in enumerate(zip(a,b)):out+=diff(x,y,path+'/'+str(i))
    elif a!=b or type(a)!=type(b):
        row={'path':path,'reference':a,'current':b}
        if type(a) in (int,float) and type(b) in (int,float):row['absolute_difference']=abs(a-b)
        out.append(row)
    return out

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
    out=args.output_dir.expanduser().resolve()
    if out.exists() or out==ROOT or ROOT in out.parents:ap.error('Use a new directory outside this package.')
    before=inventory();incoming=json.loads((ROOT/'INCOMING_SHA256.json').read_text())
    for n,h in incoming.items():
        if sha(ROOT/'prior'/n)!=h:raise RuntimeError('Changed incoming source: '+n)
    for m in (ROOT/'prior').rglob('MANIFEST.json'):
        for n,item in json.loads(m.read_text()).items():
            if sha(m.parent/n)!=item['sha256']:raise RuntimeError('Archived manifest mismatch: '+n)
    out.mkdir(parents=True)
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    rows=[]
    for name,script,reference in [('followup','check_followup.py','evidence/first.json'),('exact_series','exact_short_time.py','evidence/exact_first.json'),('pilot','prior/check_pilot.py','prior/evidence/first.json')]:
        output=out/(name+'.json')
        with (out/(name+'.log')).open('xb') as log:
            p=subprocess.run([sys.executable,str(ROOT/script),'--output',str(output)],env=env,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=120)
        data=json.loads(output.read_text()) if output.exists() else {}
        changes=diff(json.loads((ROOT/reference).read_text()),data)
        rows.append({'suite':name,'returncode':p.returncode,'status':data.get('status'),'groups':data.get('test_groups',0),'canonical_byte_identity':output.exists() and output.read_bytes()==(ROOT/reference).read_bytes(),'differences':changes})
    source_ok=before==inventory();assertions=all(r['returncode']==0 and r['status']=='PASS' for r in rows)
    identity=all(r['canonical_byte_identity'] for r in rows)
    receipt={'assertions_pass':assertions,'source_unchanged':source_ok,'all_reports_byte_identical':identity,'total_scientific_groups':sum(r['groups'] for r in rows),'incoming_files':len(incoming),'source_sha256':before,'runs':rows,'scope':'Author-side checks; not external proof, numerical interval, or device certification.'}
    with (out/'RECEIPT.json').open('x') as f:json.dump(receipt,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('source_sha256','runs')},indent=2))
    raise SystemExit(0 if assertions and source_ok and identity else 1)
if __name__=='__main__':main()
