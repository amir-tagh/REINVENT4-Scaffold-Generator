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
from rdkit.Chem import Descriptors, Draw, QED, rdMolDescriptors, FilterCatalog
from rdkit.Contrib.SA_Score import sascorer
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator
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
    for name, default in [('qed-min',0.7), ('sa-max',4.0), ('sc-max',3.5),
                          ('mw-min',150), ('mw-max',500), ('logp-min',-1), ('logp-max',5),
                          ('tpsa-max',140), ('hbd-max',5), ('hba-max',10),
                          ('rotb-max',10), ('charge-max',2), ('lipinski-max',1),
                          ('similarity-min',0), ('similarity-max',1)]:
        p.add_argument('--'+name, type=float, default=default,
                       help=f'Filter threshold (default {default})')
    p.add_argument('--scscore-model', type=Path, help='Official SCScore JSON.gz model weights')
    p.add_argument('--scscore-repo', type=Path, help='Path to cloned SCScore repository')
    p.add_argument('--allow-pains', action='store_true')
    p.add_argument('--exclude-brenk', action='store_true')
    p.add_argument('--vendor-catalog', type=Path, help='CSV with SMILES and optional ID; exact structure matching')
    p.add_argument('--vendor-only', action='store_true', help='Require exact match in supplied catalog')
    p.add_argument('--no-filters', action='store_true', help='Disable property and alert thresholds; retain scaffold and novelty checks')
    p.add_argument('--core-smiles', help='Explicit connected core to preserve in every parent')
    p.add_argument('--attachment-sites', help='Two zero-based atom indices in canonical core SMILES')
    p.add_argument('--require-smarts', action='append', default=[], help='Required substructure; repeat for multiple features')
    p.add_argument('--exclude-smarts', action='append', default=[], help='Forbidden substructure; repeat as needed')
    p.add_argument('--retro-csv', type=Path, help='External retrosynthesis results: SMILES, RouteFound (true/false), optional RAscore')
    p.add_argument('--require-route', action='store_true')
    p.add_argument('--rascore-min', type=float)
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
    explicit = Chem.MolFromSmiles(args.core_smiles) if args.core_smiles else None
    if args.core_smiles and (explicit is None or len(Chem.GetMolFrags(explicit)) != 1 or any(a.GetAtomicNum() == 0 for a in explicit.GetAtoms())):
        raise ValueError('Provide a valid connected core without dummy atoms')
    by_scaffold = defaultdict(list)
    valid_inputs = []
    for input_id, smi in inputs:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f'WARNING: invalid SMILES for {input_id}: {smi}', file=sys.stderr)
            continue
        canonical = Chem.MolToSmiles(mol)
        if explicit is not None and not mol.HasSubstructMatch(explicit, useChirality=True):
            raise ValueError(f'Explicit core not present in parent {input_id}')
        scaffold = Chem.MolToSmiles(explicit) if explicit is not None else make_scaffold(mol)
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
        if args.attachment_sites:
            sites = tuple(int(x) for x in args.attachment_sites.split(','))
            if len(sites) != 2 or len(set(sites)) != 2 or not set(sites).issubset(allowed):
                raise ValueError('Specify two distinct eligible canonical-core atom indices; see candidates or inspect canonical core')
            scores = dict(scored)
            pairs = [tuple((i, scores[i]) for i in sites)]
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
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    return {'MW': mw, 'cLogP': logp, 'TPSA': rdMolDescriptors.CalcTPSA(mol),
            'HBD': hbd, 'HBA': hba, 'QED': QED.qed(mol),
            'SA': sascorer.calculateScore(mol),
            'RotB': rdMolDescriptors.CalcNumRotatableBonds(mol),
            'FormalCharge': Chem.GetFormalCharge(mol),
            'FractionCSP3': rdMolDescriptors.CalcFractionCSP3(mol),
            'LipinskiViolations': sum([mw > 500, logp > 5, hbd > 5, hba > 10])}


class LibraryFilter:
    def __init__(self, args, inputs):
        self.args = args
        self.parents = {r['input_smiles'] for r in inputs}
        self.fp = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        self.parent_fps = defaultdict(list)
        for r in inputs:
            self.parent_fps[r['scaffold']].append(self.fp.GetFingerprint(Chem.MolFromSmiles(r['input_smiles'])))
        self.cores = {r['scaffold']: Chem.MolFromSmiles(r['scaffold']) for r in inputs}
        self.catalogs = {}
        for name in ('PAINS', 'BRENK'):
            params = FilterCatalog.FilterCatalogParams()
            params.AddCatalog(getattr(FilterCatalog.FilterCatalogParams.FilterCatalogs, name))
            self.catalogs[name] = FilterCatalog.FilterCatalog(params)
        self.sc = None
        if args.scscore_model:
            if args.scscore_repo:
                sys.path.insert(0, str(args.scscore_repo.resolve()))
            from scscore.standalone_model_numpy import SCScorer
            self.sc = SCScorer()
            self.sc.restore(str(args.scscore_model.resolve()))
        self.required = [Chem.MolFromSmarts(x) for x in args.require_smarts]
        self.forbidden = [Chem.MolFromSmarts(x) for x in args.exclude_smarts]
        if any(q is None for q in self.required + self.forbidden):
            raise ValueError('Invalid SMARTS query')
        self.retro = {}
        if (args.require_route or args.rascore_min is not None) and not args.retro_csv:
            raise ValueError('Retrosynthesis filters require --retro-csv')
        if args.retro_csv:
            with args.retro_csv.open(newline='') as f:
                reader = csv.DictReader(f)
                if not {'SMILES','RouteFound'}.issubset(reader.fieldnames or []):
                    raise ValueError('Retrosynthesis CSV requires SMILES and RouteFound')
                for row in reader:
                    m = Chem.MolFromSmiles(row['SMILES'])
                    if m is None: raise ValueError('Invalid retrosynthesis SMILES')
                    key = Chem.MolToSmiles(m)
                    status = row['RouteFound'].strip().lower()
                    if status not in ('true','false'): raise ValueError('RouteFound must be true or false')
                    score = float(row['RAscore']) if row.get('RAscore') else None
                    if score is not None and not 0 <= score <= 1: raise ValueError('RAscore outside [0,1]')
                    value = (status, score)
                    if key in self.retro and self.retro[key] != value: raise ValueError('Conflicting retrosynthesis entries')
                    self.retro[key] = value
        self.vendor = defaultdict(set)
        if args.vendor_only and not args.vendor_catalog:
            raise ValueError('--vendor-only requires --vendor-catalog')
        if args.vendor_catalog:
            with args.vendor_catalog.open(newline='') as f:
                reader = csv.DictReader(f)
                if 'SMILES' not in (reader.fieldnames or []):
                    raise ValueError('Vendor CSV requires SMILES column')
                for row in reader:
                    m = Chem.MolFromSmiles(row['SMILES'])
                    if m is not None:
                        self.vendor[Chem.MolToSmiles(m)].add(row.get('ID') or 'catalog_match')

    def evaluate(self, mol, canonical, scaffold):
        a = self.args
        p = properties(mol)
        p['SCScore'] = float(self.sc.get_score_from_smi(canonical)[1]) if self.sc else ''
        p['ParentSimilarity'] = max(DataStructs.BulkTanimotoSimilarity(self.fp.GetFingerprint(mol), self.parent_fps[scaffold]))
        p['ScaffoldRetained'] = mol.HasSubstructMatch(self.cores[scaffold], useChirality=True)
        for name, cat in self.catalogs.items():
            p[name] = ';'.join(sorted({x.GetDescription() for x in cat.GetMatches(mol)}))
        p['VendorStatus'] = ('catalog_match' if canonical in self.vendor else 'not_matched') if a.vendor_catalog else 'not_checked'
        p['VendorIDs'] = ';'.join(sorted(self.vendor.get(canonical, [])))
        route, ra = self.retro.get(canonical, ('not_checked', None))
        p['RouteFound'] = route
        p['RAscore'] = ra if ra is not None else ''
        p['RequiredFeaturesRetained'] = all(mol.HasSubstructMatch(q, useChirality=True) for q in self.required)
        p['ForbiddenFeaturesPresent'] = any(mol.HasSubstructMatch(q, useChirality=True) for q in self.forbidden)
        reasons = []
        if not p['RequiredFeaturesRetained']: reasons.append('required_feature_loss')
        if p['ForbiddenFeaturesPresent']: reasons.append('forbidden_feature')
        if a.require_route and route != 'true': reasons.append('route_not_confirmed')
        if a.rascore_min is not None and (ra is None or ra < a.rascore_min): reasons.append('RAscore')
        if canonical in self.parents: reasons.append('parent_duplicate')
        if not p['ScaffoldRetained']: reasons.append('scaffold_loss')
        if a.vendor_only and canonical not in self.vendor: reasons.append('vendor_not_matched')
        if not a.no_filters:
            checks = [('QED', p['QED'] >= a.qed_min), ('SA', p['SA'] <= a.sa_max),
                      ('MW', a.mw_min <= p['MW'] <= a.mw_max),
                      ('cLogP', a.logp_min <= p['cLogP'] <= a.logp_max),
                      ('TPSA', p['TPSA'] <= a.tpsa_max), ('HBD', p['HBD'] <= a.hbd_max),
                      ('HBA', p['HBA'] <= a.hba_max), ('RotB', p['RotB'] <= a.rotb_max),
                      ('FormalCharge', abs(p['FormalCharge']) <= a.charge_max),
                      ('Lipinski', p['LipinskiViolations'] <= a.lipinski_max),
                      ('Similarity', a.similarity_min <= p['ParentSimilarity'] <= a.similarity_max),
                      ('PAINS', a.allow_pains or not p['PAINS']),
                      ('BRENK', not a.exclude_brenk or not p['BRENK'])]
            if self.sc: checks.append(('SCScore', math.isfinite(p['SCScore']) and 0 < p['SCScore'] <= a.sc_max))
            reasons.extend(k for k, ok in checks if not ok)
        return p, reasons


def main():
    args = parse_args()
    if args.num_molecules < 1 or args.max_rounds < 1 or args.batch_size < 1:
        raise ValueError('num-molecules, max-rounds, and batch-size must be positive')
    if args.max_designs_per_scaffold < 1 or args.max_scaffolds < 1:
        raise ValueError('max-designs-per-scaffold and max-scaffolds must be positive')
    model = Path(args.model)
    if not model.is_file():
        raise FileNotFoundError(f'Model not found: {model}. Run from the REINVENT4 directory or specify --model.')
    if args.out.exists() and any(args.out.iterdir()):
        raise ValueError("Output directory must be empty or new")
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
    filters = LibraryFilter(args, valid)
    if filters.sc is None:
        print('WARNING: SCScore not enabled; supply --scscore-model to apply it.', file=sys.stderr)
    generated = {}
    rejected = []
    seen = set()
    stats = Counter()
    # A round visits each design once; stop early at requested count.
    for round_i in range(1, args.max_rounds + 1):
        for design in designs:
            remaining = args.num_molecules - len(generated)
            if remaining <= 0:
                break
            # Oversample moderately to compensate for duplicate/invalid structures.
            count = min(args.batch_size, max(16, math.ceil(remaining / len(designs) * 5)))
            run_dir = args.out / 'runs' / f'{design["design_id"]}_round_{round_i:02d}'
            raw = run_design(design, count, run_dir, args)
            for row in raw:
                stats['sampled'] += 1
                smi = row.get('SMILES', '').strip()
                mol = Chem.MolFromSmiles(smi) if smi else None
                if mol is None or len(Chem.GetMolFrags(mol)) != 1 or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
                    stats["invalid"] += 1
                    continue
                canonical = Chem.MolToSmiles(mol, canonical=True)
                key = (canonical, design['scaffold'])
                if key in seen or canonical in generated:
                    stats['duplicate'] += 1
                    continue
                seen.add(key)
                props, reasons = filters.evaluate(mol, canonical, design['scaffold'])
                record = {
                    'SMILES': canonical,
                    'design_id': design['design_id'],
                    'input_ids': design['input_ids'],
                    'scaffold': design['scaffold'],
                    'mapped_scaffold': design['mapped_scaffold'],
                    'R-groups': row.get('R-groups', ''),
                    'NLL': row.get('NLL', ''),
                    **props,
                }
                if reasons:
                    rejected.append({**record, 'RejectionReasons': ';'.join(reasons)})
                    stats['rejected'] += 1
                    stats.update('fail_' + reason for reason in reasons)
                    continue
                generated[canonical] = record
                stats['accepted'] += 1
                if len(generated) >= args.num_molecules:
                    break
            print(f'Round {round_i}, {design["design_id"]}: {len(generated)}/{args.num_molecules} unique molecules', flush=True)
        if len(generated) >= args.num_molecules:
            break
    rows = list(generated.values())
    columns = ['SMILES', 'design_id', 'input_ids', 'scaffold', 'mapped_scaffold',
               'R-groups', 'NLL', 'MW', 'cLogP', 'TPSA', 'HBD', 'HBA', 'QED',
               'SA', 'SCScore', 'RotB', 'FormalCharge', 'FractionCSP3', 'LipinskiViolations',
               'ParentSimilarity', 'ScaffoldRetained', 'PAINS', 'BRENK', 'VendorStatus', 'VendorIDs', 'RouteFound', 'RAscore', 'RequiredFeaturesRetained', 'ForbiddenFeaturesPresent']
    write_csv(args.out / 'generated.csv', rows, columns)
    write_csv(args.out / 'rejected.csv', rejected, columns + ['RejectionReasons'])
    write_csv(args.out / 'filter_summary.csv', [{'metric': k, 'count': v} for k,v in sorted(stats.items())], ['metric','count'])
    import json
    (args.out / 'filter_settings.json').write_text(json.dumps(vars(args), default=str, indent=2))
    if rows:
        try:
            mols = [Chem.MolFromSmiles(r['SMILES']) for r in rows[:24]]
            legends = [f'{r["design_id"]} | QED {r["QED"]:.2f}' for r in rows[:24]]
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
