#!/usr/bin/env python3
"""Docking-free cumulative scoring for REINVENT4 generated.csv libraries.

Dependencies: pandas, numpy, rdkit. Input is never modified.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors

RDLogger.DisableLog('rdApp.warning')
DEFAULT_WEIGHTS = {'QED': .25, 'SA': .20, 'SCScore': .15,
                   'ParentSimilarity': .15, 'Diversity': .15, 'Properties': .10}
PARENTS = [
'NCC(CN)NCC(O)COc1cccc2ccccc12',
'Cc1nn(-c2ccc(F)cc2)c2sc(C(=O)N3C[C@@H](N)[C@H](N)C3)cc12',
'NCCCNCCOc1ccc(-c2ccccc2)cc1',
'NCCCC[C@H](N)C(=O)Nc1ccc(-c2nc3ccccc3o2)cc1',
'NCCCC[C@H](N)C(=O)N1CCc2onc(-c3cccc4ccccc34)c2C1',
'Cc1cc(N2CCC(CN3CCN(C)CC3)CC2)nc2cc3c(cc12)CCC3',
'COc1cc2c(cc1OC)-c1cc(-c3ccccc3)c(C#N)c(=N)n1CC2',
'CCn1cc(C(=O)NCCN(C)C)c(=O)c2cc(F)c(N3CCN(C)CC3)cc21',
]

def clip01(x):
    return min(1., max(0., float(x)))

def trapezoid(x, low, ideal_low, ideal_high, high):
    if not math.isfinite(x): return float('nan')
    if x <= low or x >= high: return 0.
    if ideal_low <= x <= ideal_high: return 1.
    if x < ideal_low: return (x-low)/(ideal_low-low)
    return (high-x)/(high-ideal_high)

def property_score(mol):
    # Configurable starting preferences, not universal RNA-binding requirements.
    metrics = [
        (Descriptors.MolWt(mol), 150, 250, 500, 650),
        (Descriptors.MolLogP(mol), -2, 0, 4, 6),
        (rdMolDescriptors.CalcTPSA(mol), 0, 30, 120, 180),
        (Lipinski.NumHDonors(mol), -1, 0, 3, 7),
        (Lipinski.NumHAcceptors(mol), -1, 0, 8, 13),
    ]
    return float(np.mean([trapezoid(*x) for x in metrics]))

def fingerprint(mol):
    return AllChem.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(mol)

def num(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else float('nan')
    except (TypeError, ValueError):
        return float('nan')

def score_column(name, value):
    x = num(value)
    if not math.isfinite(x): return float('nan')
    if name == 'QED': return clip01(x)
    if name == 'SA': return clip01((10-x)/9)
    if name == 'SCScore': return clip01((5-x)/4)
    if name == 'ParentSimilarity':
        return math.exp(-0.5*((x-.55)/.20)**2)
    raise KeyError(name)

def hit_from_path(path):
    for parent in path.parents:
        if parent.name.startswith('Hit_'):
            return parent.name
    return 'Unknown'

def parent_for_hit(hit):
    try:
        idx = int(hit.split('_')[-1])-1
        return Chem.MolFromSmiles(PARENTS[idx]) if 0 <= idx < 8 else None
    except (ValueError, IndexError):
        return None

def load_weights(path):
    if path is None: return dict(DEFAULT_WEIGHTS)
    with open(path) as f: data = json.load(f)
    weights = data.get('weights', data)
    if set(weights) != set(DEFAULT_WEIGHTS):
        raise ValueError(f'Weights must contain exactly: {list(DEFAULT_WEIGHTS)}')
    weights = {k:float(v) for k,v in weights.items()}
    if any(not math.isfinite(v) or v < 0 for v in weights.values()) or sum(weights.values()) <= 0:
        raise ValueError('Weights must be finite, nonnegative, and sum to >0')
    return {k:v/sum(weights.values()) for k,v in weights.items()}

def run(args):
    weights = load_weights(args.weights)
    inp = Path(args.input)
    files = [inp] if inp.is_file() else sorted(inp.rglob('generated.csv'))
    if not files: raise FileNotFoundError(f'No generated.csv files found under {inp}')
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    all_frames=[]
    for file in files:
        df = pd.read_csv(file, low_memory=False)
        if 'SMILES' not in df: raise ValueError(f'Missing SMILES: {file}')
        df['Hit'] = hit_from_path(file)
        df['SourceFile'] = str(file)
        all_frames.append(df)
    df = pd.concat(all_frames, ignore_index=True)
    # Fingerprints retained only in memory, not written into CSV.
    mols = [Chem.MolFromSmiles(str(s)) if pd.notna(s) else None for s in df['SMILES']]
    fps = [fingerprint(m) if m is not None else None for m in mols]
    parent_fps = {}
    for hit in df['Hit'].unique():
        p = parent_for_hit(hit)
        parent_fps[hit] = fingerprint(p) if p else None
    normalized = {k:[] for k in DEFAULT_WEIGHTS if k != 'Diversity'}
    for i, (m,fp) in enumerate(zip(mols,fps)):
        row = df.iloc[i]
        for k in ('QED','SA','SCScore'):
            normalized[k].append(score_column(k, row.get(k, float('nan'))))
        parent_fp = parent_fps[row['Hit']]
        similarity = (DataStructs.TanimotoSimilarity(fp,parent_fp)
                      if fp is not None and parent_fp is not None else num(row.get('ParentSimilarity')))
        normalized['ParentSimilarity'].append(score_column('ParentSimilarity',similarity))
        normalized['Properties'].append(property_score(m) if m is not None else float('nan'))
    for k, vals in normalized.items(): df['Score_'+k] = vals
    df['ValidSMILES'] = [m is not None for m in mols]
    df['MissingComponents'] = [','.join(k for k in normalized if not math.isfinite(num(df.at[i,'Score_'+k]))) for i in range(len(df))]
    # Static chemical-quality score; diversity is intentionally excluded here.
    intrinsic_weights = {k:v for k,v in weights.items() if k != 'Diversity'}
    intrinsic_sum = sum(intrinsic_weights.values())
    if intrinsic_sum <= 0: raise ValueError('At least one non-diversity weight must be positive')
    def intrinsic(i):
        if not df.at[i,'ValidSMILES']: return float('nan')
        values = [(w,num(df.at[i,'Score_'+k])) for k,w in intrinsic_weights.items() if w > 0]
        if args.missing == 'reject' and any(not math.isfinite(x) for w,x in values): return float('nan')
        return sum(w*(x if math.isfinite(x) else 0) for w,x in values)/intrinsic_sum
    df['IntrinsicScore'] = [intrinsic(i) for i in range(len(df))]
    df['CumulativeScore'] = np.nan
    df['Score_Diversity'] = np.nan
    df['SelectionRank'] = pd.Series([pd.NA]*len(df), dtype='Int64')
    df['Selected'] = False
    # Greedy diversity-aware ranking, independently for each hit. At each step,
    # score diversity against the molecules already selected in that hit.
    for hit, group in df.groupby('Hit', sort=True):
        candidates = [int(i) for i in group.index if math.isfinite(num(df.at[i,'IntrinsicScore']))]
        selected=[]
        while candidates and len(selected) < args.select_per_hit:
            best_i=None; best_s=-float('inf'); best_div=0.
            for i in candidates:
                div = 1.-max(DataStructs.BulkTanimotoSimilarity(fps[i], [fps[j] for j in selected])) if selected else 1.
                intrinsic_value = df.at[i,'IntrinsicScore']
                s = intrinsic_sum*intrinsic_value + weights['Diversity']*div
                if s > best_s + 1e-12 or (abs(s-best_s)<=1e-12 and (best_i is None or i<best_i)):
                    best_i,best_s,best_div=i,s,div
            df.at[best_i,'CumulativeScore'] = 100*best_s
            df.at[best_i,'Score_Diversity'] = best_div
            df.at[best_i,'SelectionRank'] = len(selected)+1
            df.at[best_i,'Selected'] = True
            selected.append(best_i)
            candidates.remove(best_i)
        print(f'{hit}: total={len(group)}, valid_scored={group.IntrinsicScore.notna().sum()}, selected={len(selected)}')
    # Every valid molecule gets a stable, comparable static score; only selected
    # molecules receive a dynamic diversity-aware cumulative score.
    df['IntrinsicScore'] *= 100
    df.to_csv(out/'scored_molecules.csv', index=False)
    selected = df[df['Selected']].sort_values(['Hit','SelectionRank'])
    selected.to_csv(out/'prioritized_compounds.csv', index=False)
    summary = df.groupby('Hit',as_index=False).agg(Molecules=('SMILES','size'),
             ValidScored=('IntrinsicScore','count'),Selected=('Selected','sum'),
             MeanIntrinsicScore=('IntrinsicScore','mean'))
    summary.to_csv(out/'summary.csv',index=False)
    with open(out/'scoring_config.json','w') as f:
        json.dump({'weights':weights,'missing':args.missing,'select_per_hit':args.select_per_hit,
                   'score_scale':'0-100','notes':{'IntrinsicScore':'static, non-diversity components renormalized',
                   'CumulativeScore':'dynamic, assigned only to greedily selected molecules',
                   'Diversity':'1 minus max Tanimoto to previously selected molecules; first selection = 1',
                   'ParentSimilarity':'Gaussian desirability centered at 0.55, sigma 0.20'}},f,indent=2)
    print(f'COMPLETE: {len(df)} molecules; {len(selected)} selected; output={out}')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',required=True,help='Root containing Hit_*/generated.csv or one CSV')
    ap.add_argument('--output',required=True)
    ap.add_argument('--weights',help='Optional JSON weights file')
    ap.add_argument('--select-per-hit',type=int,default=30)
    ap.add_argument('--missing',choices=['reject','zero'],default='reject')
    a=ap.parse_args()
    if a.select_per_hit < 1: ap.error('--select-per-hit must be >= 1')
    run(a)
