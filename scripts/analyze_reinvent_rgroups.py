#!/usr/bin/env python3
"""Post-generation RDKit R-group analysis and structure galleries for REINVENT4."""
import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw, rdRGroupDecomposition, rdDepictor
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFingerprintGenerator, DataStructs

RDLogger.DisableLog('rdApp.warning')


def mol_from_smiles(value):
    if pd.isna(value):
        return None
    return Chem.MolFromSmiles(str(value))


def gallery(df, path, title, max_mols=24):
    mols, legends = [], []
    for _, row in df.head(max_mols).iterrows():
        mol = mol_from_smiles(row['SMILES'])
        if mol is None:
            continue
        mols.append(mol)
        label = str(row.get('Hit', ''))
        for key in ('QED', 'SA', 'SCScore'):
            val = pd.to_numeric(row.get(key), errors='coerce')
            if pd.notna(val):
                label += f'  {key}:{val:.2f}'
        legends.append(label.strip())
    if mols:
        img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(360, 260),
                                   legends=legends, useSVG=False)
        img.save(str(path))


def run(root, output, top_n, max_cores, min_group):
    output.mkdir(parents=True, exist_ok=True)
    tables = output / 'tables'
    galleries = output / 'galleries'
    tables.mkdir(exist_ok=True)
    galleries.mkdir(exist_ok=True)
    summary, all_assignments, all_freqs = [], [], []
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    for csv_file in sorted(root.glob('Hit_*/generated.csv')):
        hit = csv_file.parent.name
        df = pd.read_csv(csv_file, keep_default_na=True)
        if 'SMILES' not in df.columns:
            print(f'{hit}: missing SMILES, skipped')
            continue
        df['Hit'] = hit
        df['Mol'] = df['SMILES'].map(mol_from_smiles)
        invalid = int(df['Mol'].isna().sum())
        df = df[df['Mol'].notna()].copy().reset_index(drop=True)
        if df.empty:
            continue
        # Group by the actual mapped generation core, not by the Murcko scaffold.
        core_col = 'mapped_scaffold' if 'mapped_scaffold' in df else 'scaffold'
        if core_col not in df:
            raise ValueError(f'{hit}: no scaffold or mapped_scaffold column')
        df['_core'] = df[core_col].fillna('').astype(str)
        counts = df['_core'].value_counts()
        cores = [(core, int(n)) for core, n in counts.items() if core and core != 'nan'][:max_cores]
        print(f'{hit}: {len(df)} valid, {invalid} invalid, {len(cores)} cores selected')
        assigned_indices = set()
        for ci, (core_smi, ncore) in enumerate(cores, start=1):
            group = df[df['_core'] == core_smi].copy()
            core = Chem.MolFromSmiles(core_smi)
            core_id = f'{hit}_core{ci:02d}'
            if core is None or len(group) < min_group:
                summary.append(dict(Hit=hit, CoreID=core_id, CoreSMILES=core_smi,
                                    GroupSize=len(group), Matched=0, Unmatched=len(group),
                                    Status='invalid_core_or_small_group'))
                continue
            # RDKit performs an explicit core match; failed molecules are recorded.
            params = rdRGroupDecomposition.RGroupDecompositionParameters()
            params.removeHydrogensPostMatch = True
            params.onlyMatchAtRGroups = False
            rgd = rdRGroupDecomposition.RGroupDecomposition([core], params)
            accepted = []
            for idx, row in group.iterrows():
                try:
                    result = rgd.Add(row['Mol'])
                except Exception:
                    result = -1
                if result >= 0:
                    accepted.append(idx)
            matched = len(accepted)
            if matched:
                try:
                    rgd.Process()
                    rows = rgd.GetRGroupsAsRows(asSmiles=True)
                except Exception as exc:
                    print(f'  {core_id}: decomposition failed: {exc}')
                    rows = []
                if len(rows) != len(accepted):
                    print(f'  {core_id}: unexpected output length, skipping assignments')
                    rows = []
                for idx, rgroup in zip(accepted, rows):
                    rec = {'Hit': hit, 'CoreID': core_id, 'SourceRow': idx + 2,
                           'SMILES': df.loc[idx, 'SMILES']}
                    for key, val in rgroup.items():
                        rec[key] = val
                    all_assignments.append(rec)
                    assigned_indices.add(idx)
                if rows:
                    rcols = sorted({key for row in rows for key in row if key != 'Core'})
                    for rcol in rcols:
                        vals = [str(row.get(rcol, '')) for row in rows]
                        for val, freq in Counter(vals).most_common():
                            all_freqs.append({'Hit': hit, 'CoreID': core_id, 'RGroup': rcol,
                                              'SubstituentSMILES': val, 'Count': freq,
                                              'PercentOfMatched': round(100 * freq / len(rows), 3)})
            summary.append(dict(Hit=hit, CoreID=core_id, CoreSMILES=core_smi,
                                GroupSize=len(group), Matched=matched,
                                Unmatched=len(group)-matched,
                                Status='ok' if matched else 'no_matches'))
        # Representative compounds: high QED, low SA and SCScore; deterministic ranking.
        for col in ('QED', 'SA', 'SCScore'):
            df[col] = pd.to_numeric(df[col], errors='coerce') if col in df else float('nan')
        df['_rank'] = (df['QED'].rank(pct=True, na_option='bottom')
                       + (1 - df['SA'].rank(pct=True, na_option='bottom'))
                       + (1 - df['SCScore'].rank(pct=True, na_option='bottom')))
        best = df.sort_values(['_rank', 'SMILES'], ascending=[False, True]).head(top_n)
        best.drop(columns=['Mol', '_core', '_rank']).to_csv(tables / f'{hit}_top_candidates.csv', index=False)
        gallery(best, galleries / f'{hit}_top_candidates.png', hit, top_n)
        # Diverse representatives selected greedily from a high-QED candidate pool.
        pool = df.sort_values(['QED', 'SMILES'], ascending=[False, True]).head(min(250, len(df)))
        fps = [generator.GetFingerprint(m) for m in pool['Mol']]
        chosen = [0]
        while len(chosen) < min(top_n, len(pool)):
            remaining = [i for i in range(len(pool)) if i not in chosen]
            best_idx = max(remaining, key=lambda i: (min(1 - DataStructs.TanimotoSimilarity(fps[i], fps[j]) for j in chosen), -i))
            chosen.append(best_idx)
        diverse = pool.iloc[chosen]
        diverse.drop(columns=['Mol', '_core', '_rank']).to_csv(tables / f'{hit}_diverse_candidates.csv', index=False)
        gallery(diverse, galleries / f'{hit}_diverse_candidates.png', hit, top_n)
    pd.DataFrame(summary).to_csv(tables / 'core_decomposition_summary.csv', index=False)
    pd.DataFrame(all_assignments).to_csv(tables / 'rgroup_assignments.csv', index=False)
    pd.DataFrame(all_freqs).to_csv(tables / 'rgroup_frequencies.csv', index=False)
    (output / 'settings.json').write_text(json.dumps({'input': str(root), 'top_n': top_n,
         'max_cores_per_hit': max_cores, 'min_group': min_group}, indent=2))
    print(f'Done. Outputs: {output}')
    print('Note: R-group frequencies are per mapped core, not globally comparable positions.')

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--input', default='results/v5/REINVENT_V5_All_8_Hits_Relaxed')
    p.add_argument('--output', default='results/v5/REINVENT_V5_RGroup_Analysis')
    p.add_argument('--top-n', type=int, default=24)
    p.add_argument('--max-cores', type=int, default=3)
    p.add_argument('--min-group', type=int, default=5)
    args = p.parse_args()
    run(Path(args.input), Path(args.output), args.top_n, args.max_cores, args.min_group)
