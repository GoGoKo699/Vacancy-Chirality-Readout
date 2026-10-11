#!/usr/bin/env python3
"""Verify sources, five new resource groups and all eleven original groups.

Output must be a new directory outside the package. Canonical byte identity,
scientific assertion results and source integrity are recorded separately.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parent


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def inventory():
    out={}
    for p in sorted(ROOT.rglob('*')):
        if p.is_symlink():raise ValueError('Unexpected source symlink: '+str(p))
        if p.is_file():out[p.relative_to(ROOT).as_posix()]=digest(p)
    return out


def differences(a,b,path=''):
    rows=[]
    if isinstance(a,dict) and isinstance(b,dict):
        if a.keys()!=b.keys():rows.append({'path':path,'kind':'keys'})
        for k in sorted(a.keys()&b.keys()):rows+=differences(a[k],b[k],path+'/'+str(k))
    elif isinstance(a,list) and isinstance(b,list):
        if len(a)!=len(b):rows.append({'path':path,'kind':'length'})
        for i,(x,y) in enumerate(zip(a,b)):rows+=differences(x,y,path+'/'+str(i))
    elif a!=b or type(a)!=type(b):
        item={'path':path,'reference':a,'rerun':b}
        if type(a) in (float,int) and type(b) in (float,int):item['absolute_difference']=abs(a-b)
        rows.append(item)
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
    out=args.output_dir.expanduser().resolve()
    if out.exists() or out==ROOT or ROOT in out.parents:ap.error('Use a new output directory outside this package.')
    before=inventory()
    incoming=json.loads((ROOT/'INCOMING_HASHES.json').read_text())
    actual={p.relative_to(ROOT/'prior').as_posix():digest(p) for p in (ROOT/'prior').rglob('*') if p.is_file()}
    if incoming!=actual:raise RuntimeError('Incoming archive membership or bytes changed.')
    manifests=[]
    for m in sorted((ROOT/'prior').rglob('*MANIFEST*.json')):
        record=json.loads(m.read_text())
        for n,v in record.items():
            p=(m.parent/n).resolve()
            if not p.is_relative_to(ROOT/'prior'):raise RuntimeError('Archive manifest traversal.')
            if digest(p)!=v['sha256'] or p.stat().st_size!=v['bytes']:raise RuntimeError('Archive manifest mismatch: '+n)
        manifests.append({'path':m.relative_to(ROOT).as_posix(),'members':len(record)})
    out.mkdir(parents=True)
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    with (out/'resources.log').open('xb') as f:
        new=subprocess.run([sys.executable,str(ROOT/'check_resources.py'),'--output',str(out/'resources.json')],env=env,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=120)
    with (out/'preserved_wrapper.log').open('xb') as f:
        old=subprocess.run([sys.executable,str(ROOT/'prior/verify.py'),'--output-dir',str(out/'preserved')],env=env,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=400)
    data=json.loads((out/'resources.json').read_text()) if (out/'resources.json').exists() else {}
    oldpath=out/'preserved/RECEIPT.json';previous=json.loads(oldpath.read_text()) if oldpath.exists() else {}
    ref=ROOT/'evidence/new_first.json'
    identity=(out/'resources.json').exists() and (out/'resources.json').read_bytes()==ref.read_bytes()
    diff=differences(json.loads(ref.read_text()),data)
    newpass=new.returncode==0 and data.get('status')=='PASS' and data.get('test_groups')==5
    oldpass=old.returncode==0 and previous.get('assertions_pass') is True and previous.get('source_unchanged') is True and previous.get('total_scientific_groups')==11
    sourceok=before==inventory()
    allidentity=identity and previous.get('all_reports_byte_identical') is True
    receipt={'new_assertions_pass':newpass,'preserved_assertions_pass':oldpass,
             'all_reports_byte_identical':allidentity,'new_report_differences':diff,
             'total_scientific_groups':data.get('test_groups',0)+previous.get('total_scientific_groups',0),
             'source_unchanged':sourceok,'incoming_files':len(incoming),'verified_manifests':manifests,
             'new_returncode':new.returncode,'preserved_wrapper_returncode':old.returncode,'source_sha256':before,
             'scope':'Author-side verification only; no apparatus validation, interval certification, external theorem audit, or priority clearance.'}
    passed=newpass and oldpass and sourceok and allidentity
    receipt['status']='PASS' if passed else 'FAIL'
    with (out/'RECEIPT.json').open('x') as f:json.dump(receipt,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('source_sha256','new_report_differences')},indent=2))
    raise SystemExit(0 if passed else 1)
if __name__=='__main__':main()
