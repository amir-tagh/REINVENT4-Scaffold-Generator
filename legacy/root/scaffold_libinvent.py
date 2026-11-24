#!/usr/bin/env python3
"""Generate LibInvent analogs from SMILES using RDKit-selected scaffold sites.

Tested design targets REINVENT4 sampling TOML and its LibInvent prior.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import math
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, QED, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--smiles', help='Single input SMILES')
    group.add_argument('--input', type=Path, help='Input .smi/.txt or .csv file')
    p.add_argument('--smiles-column', default='SMILES', help='CSV column containing SMILES')
    p.add_argument('--num-molecules', '-n', type=int, default=100, help='Target number of unique valid molecules')
    p.add_argument('--model', default='priors/libinvent.prior', help='Path to LibInvent prior')
    p.add_argument('--out', type=Path, default=Path('libinvent_results'))
    p.add_argument('--device', default='cuda:0', help='Torch device (cuda:0 or cpu)')
    p.add_argument('--reinvent-executable', default='reinvent')
    p.add_argument('--max-designs-per-scaffold', type=int, default=3)
    p.add_argument('--max-scaffolds', type=int, default=30)
    p.add_argument('--max-rounds', type=int, default=5)
    p.add_argument('--batch-size', type=int, default=128, help='Sampling count per design per round')
    p.add_argument('--seed', type=int, default=None, help='Optional RNG seed for REINVENT')
    p.add_argument('--include-only-existing-sites', action='store_true', help='Use only core atoms bearing substituents in at least one input')
    return p.parse_args()


def load_inputs(args):
    if args.smiles:
        return [('input_1', args.smiles)]
    path = args.input
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.suffix.lower() == '.csv':
        with path.open(newline='') as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames or args.smiles_column not in reader.fieldnames:
                raise ValueError(f'CSV must contain column {args.smiles_column!r}')
            return [(str(row.get('ID') or row.get('id') or f'input_{i}'), row[args.smiles_column].strip())
                    for i, row in enumerate(reader, 1) if row.get(args.smiles_column, '').strip()]
    rows = []
    for i, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        parts = line.split()
        rows.append((parts[1] if len(parts) > 1 else f'input_{i}', parts[0]))
    return rows


def write_csv(path, rows, columns):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def make_scaffold(mol):
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(scaffold, canonical=True)


def get_existing_substituted_sites(mol, core):
    """Return core atom indices with external neighbors for any substructure match."""
    positions = set()
    for match in mol.GetSubstructMatches(core, uniquify=True, maxMatches=100):
        in_core = set(match)
        for i, mol_i in enumerate(match):
            if any(n.GetIdx() not in in_core for n in mol.GetAtomWithIdx(mol_i).GetNeighbors()):
                positions.add(i)
    return positions


def eligible_sites(core):
    """Sites with a replaceable implicit/explicit hydrogen or external valence."""
    sites = []
    for atom in core.GetAtoms():
        if atom.GetAtomicNum() not in (6, 7, 8, 16):
            continue
        if atom.GetTotalNumHs() > 0:
            sites.append(atom.GetIdx())
    return sites


def make_design(core_smiles, sites):
    """Append two mapped dummy atoms to selected core atoms, preserving the core."""
    core = Chem.MolFromSmiles(core_smiles)
    if core is None or len(sites) != 2 or sites[0] == sites[1]:
        return None
    rw = Chem.RWMol(core)
    for map_no, site in enumerate(sites):
        atom = Chem.Atom(0)
        atom.SetAtomMapNum(map_no)
        idx = rw.AddAtom(atom)
        rw.AddBond(site, idx, Chem.BondType.SINGLE)
    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        result = Chem.MolToSmiles(mol, canonical=True)
        check = Chem.MolFromSmiles(result)
        if check is None or sum(a.GetAtomicNum() == 0 for a in check.GetAtoms()) != 2:
            return None
        return result
    except Exception:
        return None


def choose_designs(inputs, args):
    by_scaffold = defaultdict(list)
    valid_inputs = []
    for input_id, smi in inputs:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f'WARNING: invalid SMILES for {input_id}: {smi}', file=sys.stderr)
            continue
        canonical = Chem.MolToSmiles(mol)
        scaffold = make_scaffold(mol)
        if not scaffold:
            print(f'WARNING: no Bemis-Murcko scaffold for {input_id}; skipped', file=sys.stderr)
            continue
        by_scaffold[scaffold].append((input_id, mol))
        valid_inputs.append({'input_id': input_id, 'input_smiles': canonical, 'scaffold': scaffold})

    # Larger input scaffold families first; preserve deterministic ordering.
    families = sorted(by_scaffold.items(), key=lambda kv: (-len(kv[1]), kv[0]))[:args.max_scaffolds]
    candidates = []
    designs = []
    for scaffold, members in families:
        core = Chem.MolFromSmiles(scaffold)
        existing = Counter()
        for _, mol in members:
            existing.update(get_existing_substituted_sites(mol, core))
        allowed = eligible_sites(core)
        if args.include_only_existing_sites:
            allowed = [i for i in allowed if existing[i] > 0]
        scored = []
        for site in allowed:
            atom = core.GetAtomWithIdx(site)
            # Explainable heuristic, not a binding or synthesizability predictor.
            score = 10 * existing[site] + (2 if atom.GetAtomicNum() == 6 else 0)
            scored.append((site, score))
            candidates.append({
                'scaffold': scaffold, 'site_atom_index_0_based': site,
                'atom_symbol': atom.GetSymbol(), 'existing_substituent_count': existing[site],
                'score': score,
                'rationale': 'Observed substituent in input series' if existing[site] else 'Core atom with available H',
                'input_ids': ';'.join(m[0] for m in members),
            })
        ranked = sorted(scored, key=lambda item: (-item[1], item[0]))
        pairs = sorted(itertools.combinations(ranked, 2),
                       key=lambda pair: (-(pair[0][1] + pair[1][1]),
                                         min(pair[0][0], pair[1][0]), max(pair[0][0], pair[1][0])))
        selected = 0
        for pair in pairs:
            a, b = sorted([pair[0][0], pair[1][0]])
            mapped = make_design(scaffold, (a, b))
            if mapped is None:
                continue
            designs.append({
                'design_id': f'design_{len(designs)+1:04d}',
                'scaffold': scaffold, 'mapped_scaffold': mapped,
                'site_0': a, 'site_1': b,
                'input_ids': ';'.join(m[0] for m in members),
                'input_count': len(members),
            })
            selected += 1
            if selected >= args.max_designs_per_scaffold:
                break
        if selected == 0:
            print(f'WARNING: no valid two-point design for scaffold {scaffold}', file=sys.stderr)
    return valid_inputs, candidates, designs


def quote_toml(value):
    return '"' + str(value).replace('\\', '\\\\').replace('"', '\\"') + '"'


def run_design(design, count, run_dir, args):
    run_dir.mkdir(parents=True, exist_ok=True)
    scaffold_file = (run_dir / 'scaffold.smi').resolve()
    output_file = (run_dir / 'sampled.csv').resolve()
    config_file = (run_dir / 'sampling.toml').resolve()
    log_file = (run_dir / 'sampling.log').resolve()
    scaffold_file.write_text(design['mapped_scaffold'] + '\n')
    config = [
        'run_type = "sampling"',
        f'device = {quote_toml(args.device)}',
        '', '[parameters]',
        f'model_file = {quote_toml(Path(args.model).resolve())}',
        f'smiles_file = {quote_toml(scaffold_file)}',
        f'output_file = {quote_toml(output_file)}',
        f'num_smiles = {count}',
        'unique_molecules = true',
        'randomize_smiles = true',
    ]
    if args.seed is not None:
        # Sampling seed is deliberately not set here because the supported
        # REINVENT4 seed parameter varies across releases.
        print('WARNING: --seed provided, but not applied; REINVENT4 version-specific seed setting not verified', file=sys.stderr)
    config_file.write_text('\n'.join(config) + '\n')
    result = subprocess.run([args.reinvent_executable, '-l', str(log_file), str(config_file)],
                            text=True, capture_output=True)
    (run_dir / 'stdout.txt').write_text(result.stdout)
    (run_dir / 'stderr.txt').write_text(result.stderr)
    if result.returncode != 0 or not output_file.is_file():
        print(f'WARNING: failed {design["design_id"]}: {result.stderr[-500:]} (see {run_dir})', file=sys.stderr)
        return []
    with output_file.open(newline='') as f:
        return list(csv.DictReader(f))


def properties(mol):
    return {
        'MW': round(Descriptors.MolWt(mol), 2),
        'cLogP': round(Descriptors.MolLogP(mol), 2),
        'TPSA': round(rdMolDescriptors.CalcTPSA(mol), 2),
        'HBD': rdMolDescriptors.CalcNumHBD(mol),
        'HBA': rdMolDescriptors.CalcNumHBA(mol),
        'QED': round(QED.qed(mol), 3),
    }


def main():
    args = parse_args()
    if args.num_molecules < 1 or args.max_rounds < 1 or args.batch_size < 1:
        raise ValueError('num-molecules, max-rounds, and batch-size must be positive')
    if args.max_designs_per_scaffold < 1 or args.max_scaffolds < 1:
        raise ValueError('max-designs-per-scaffold and max-scaffolds must be positive')
    model = Path(args.model)
    if not model.is_file():
        raise FileNotFoundError(f'Model not found: {model}. Run from the REINVENT4 directory or specify --model.')
    args.out.mkdir(parents=True, exist_ok=True)
    valid, candidates, designs = choose_designs(load_inputs(args), args)
    write_csv(args.out / 'inputs.csv', valid, ['input_id', 'input_smiles', 'scaffold'])
    write_csv(args.out / 'candidates.csv', candidates,
              ['scaffold', 'site_atom_index_0_based', 'atom_symbol', 'existing_substituent_count',
               'score', 'rationale', 'input_ids'])
    write_csv(args.out / 'designs.csv', designs,
              ['design_id', 'scaffold', 'mapped_scaffold', 'site_0', 'site_1', 'input_ids', 'input_count'])
    if not designs:
        raise SystemExit('No valid two-point LibInvent scaffold designs. See candidates.csv and warnings.')
    print(f'{len(valid)} valid inputs; {len(designs)} scaffold designs')
    generated = {}
    # A round visits each design once; stop early at requested count.
    for round_i in range(1, args.max_rounds + 1):
        for design in designs:
            remaining = args.num_molecules - len(generated)
            if remaining <= 0:
                break
            # Oversample moderately to compensate for duplicate/invalid structures.
            count = min(args.batch_size, max(16, math.ceil(remaining / len(designs) * 1.5)))
            run_dir = args.out / 'runs' / f'{design["design_id"]}_round_{round_i:02d}'
            raw = run_design(design, count, run_dir, args)
            for row in raw:
                smi = row.get('SMILES', '').strip()
                mol = Chem.MolFromSmiles(smi) if smi else None
                if mol is None or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
                    continue
                canonical = Chem.MolToSmiles(mol, canonical=True)
                if canonical in generated:
                    continue
                generated[canonical] = {
                    'SMILES': canonical,
                    'design_id': design['design_id'],
                    'input_ids': design['input_ids'],
                    'scaffold': design['scaffold'],
                    'mapped_scaffold': design['mapped_scaffold'],
                    'R-groups': row.get('R-groups', ''),
                    'NLL': row.get('NLL', ''),
                    **properties(mol),
                }
                if len(generated) >= args.num_molecules:
                    break
            print(f'Round {round_i}, {design["design_id"]}: {len(generated)}/{args.num_molecules} unique molecules', flush=True)
        if len(generated) >= args.num_molecules:
            break
    rows = list(generated.values())
    columns = ['SMILES', 'design_id', 'input_ids', 'scaffold', 'mapped_scaffold',
               'R-groups', 'NLL', 'MW', 'cLogP', 'TPSA', 'HBD', 'HBA', 'QED']
    write_csv(args.out / 'generated.csv', rows, columns)
    if rows:
        try:
            mols = [Chem.MolFromSmiles(r['SMILES']) for r in rows[:24]]
            legends = [f'{r["design_id"]} | QED {r["QED"]}' for r in rows[:24]]
            img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(340, 260),
                                       legends=legends, useSVG=False)
            img.save(str(args.out / 'molecules.png'))
        except Exception as exc:
            print(f'WARNING: could not create structure image: {exc}', file=sys.stderr)
    print(f'Done: {len(rows)} unique molecules saved to {args.out / "generated.csv"}')
    if len(rows) < args.num_molecules:
        print('WARNING: target not reached; try more rounds or additional designs.', file=sys.stderr)


if __name__ == '__main__':
    main()
