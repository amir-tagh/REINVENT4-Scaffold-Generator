#!/usr/bin/env python3
"""Generate LibInvent analogs from SMILES using RDKit-selected scaffold sites.

Tested design targets REINVENT4 sampling TOML and its LibInvent prior.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import math
import time
import json
import copy
import re
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, QED, rdMolDescriptors, FilterCatalog
from rdkit.Contrib.SA_Score import sascorer
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--smiles', help='Single input SMILES')
    group.add_argument('--input', type=Path, help='Input .smi/.txt or .csv file')
    p.add_argument('--input-charge-policy', choices=['neutralize', 'preserve'], default='neutralize',
                   help='Neutralize removable formal charges on input molecules; preserve original SMILES in audit file')
    p.add_argument('--input-largest-fragment', action='store_true',
                   help='Optionally retain largest organic fragment (e.g., remove salt counterions)')
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
    p.add_argument("--min-samples", type=int, default=128)
    p.add_argument("--resume", action="store_true")
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



def standardize_parent_smiles(smiles, args):
    """Sanitize and optionally neutralize a parent without changing the recorded original."""
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        raise ValueError('RDKit could not parse or sanitize input SMILES')
    original_charge = Chem.GetFormalCharge(mol)
    original_fragments = len(Chem.GetMolFrags(mol))
    if args.input_largest_fragment:
        mol = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(mol)
    if args.input_charge_policy == 'neutralize':
        mol = rdMolStandardize.Uncharger().uncharge(mol)
    Chem.SanitizeMol(mol)
    normalized = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    check = Chem.MolFromSmiles(normalized)
    if check is None:
        raise ValueError('Standardized SMILES failed round-trip validation')
    return normalized, {
        'OriginalSMILES': smiles,
        'StandardizedSMILES': normalized,
        'OriginalFormalCharge': original_charge,
        'StandardizedFormalCharge': Chem.GetFormalCharge(check),
        'OriginalFragments': original_fragments,
        'StandardizedFragments': len(Chem.GetMolFrags(check)),
        'ChargePolicy': args.input_charge_policy,
        'LargestFragment': args.input_largest_fragment,
        'Changed': Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True, isomericSmiles=True) != normalized,
        'Status': 'ok',
        'Error': '',
    }


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


class PersistentSampler:
    def __init__(self, args):
        import torch
        from reinvent.runmodes.create_adapter import create_adapter
        from reinvent.runmodes.setup_sampler import setup_sampler
        self.torch = torch
        self.setup = setup_sampler
        self.device = torch.device(args.device)
        if args.seed is not None:
            random.seed(args.seed)
            import numpy as np
            np.random.seed(args.seed)
            torch.manual_seed(args.seed)
        start = time.perf_counter()
        self.adapter, _, self.model_type = create_adapter(str(Path(args.model).resolve()), 'inference', self.device)
        # REINVENT 4.8.24 LibinventSampler obtains scaffold token IDs from a
        # CPU DataLoader. Its LibinventAdapter.sample forwards them unchanged,
        # while the network's embedding weights are on CUDA. Move the token
        # tensor to the model device at the adapter boundary. Keep sequence
        # lengths on CPU for RNN packing compatibility.
        if self.model_type == 'Libinvent':
            original_sample = self.adapter.sample

            def sample_on_model_device(scaffold_seqs, scaffold_seq_lengths):
                scaffold_seqs = scaffold_seqs.to(self.device)
                return original_sample(scaffold_seqs, scaffold_seq_lengths)

            self.adapter.sample = sample_on_model_device

            # REINVENT4 4.8.24 LibInvent RNN can construct initial LSTM
            # hidden/cell states on CPU even when its encoder is on CUDA.
            # A forward pre-hook reconciles the hidden state device without
            # modifying the installed REINVENT source or the prior.
            def align_lstm_state(module, inputs):
                if len(inputs) < 2 or inputs[1] is None:
                    return None
                device = next(module.parameters()).device
                hidden = inputs[1]
                if not isinstance(hidden, tuple) or len(hidden) != 2:
                    return None
                moved = tuple(t.to(device) for t in hidden)
                return (inputs[0], moved, *inputs[2:])

            self._lstm_device_hooks = [
                module.register_forward_pre_hook(align_lstm_state)
                for module in self.adapter.model.network.modules()
                if isinstance(module, torch.nn.LSTM)
            ]
            # Decoder autoregressive token IDs are allocated on CPU by the
            # older LibInvent sampling implementation. Move embedding indices
            # to the embedding weight device at every encoder/decoder call.
            # This leaves CPU sequence lengths untouched for pack_padded_sequence.
            def align_embedding_indices(module, inputs):
                if not inputs:
                    return None
                indices = inputs[0]
                device = module.weight.device
                if isinstance(indices, torch.Tensor) and indices.device != device:
                    return (indices.to(device), *inputs[1:])
                return None

            self._embedding_device_hooks = [
                module.register_forward_pre_hook(align_embedding_indices)
                for module in self.adapter.model.network.modules()
                if isinstance(module, torch.nn.Embedding)
            ]
            print(f'Installed device-alignment hooks on {len(self._lstm_device_hooks)} LSTM and '
                  f'{len(self._embedding_device_hooks)} embedding modules', flush=True)
        if self.model_type not in ('Libinvent', 'LibinventTransformer'):
            raise ValueError(f'Unexpected model type: {self.model_type}')
        self.load_seconds = time.perf_counter() - start
        print(f'Model loaded once: {self.model_type}, {self.load_seconds:.2f} s', flush=True)

    def sample(self, design, count, run_dir):
        from reinvent.models.model_factory.sample_batch import SmilesState
        from reinvent_plugins.normalizers.rdkit_smiles import normalize
        run_dir.mkdir(parents=True, exist_ok=True)
        scaffold = design['mapped_scaffold']
        (run_dir / 'scaffold.smi').write_text(scaffold + '\n')
        sampler, _ = self.setup(self.model_type, dict(batch_size=count, randomize_smiles=True), self.adapter)
        with self.torch.no_grad():
            batch = sampler.sample([scaffold])
        smiles = normalize(batch.smilies, keep_all=True) if self.model_type == 'Libinvent' else batch.smilies
        rows = [dict(SMILES=s, **{'SMILES_state': state.value, 'Scaffold': core, 'R-groups': rg, 'NLL': round(nll,2)})
                for s,state,core,rg,nll in zip(smiles,batch.states,batch.items1,batch.items2,batch.nlls.cpu().tolist())
                if state == SmilesState.VALID]
        write_csv(run_dir / 'sampled.csv', rows, ['SMILES','SMILES_state','Scaffold','R-groups','NLL'])
        return rows


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
    def __init__(self, args, inputs, shared=None):
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
        self.sc = shared.sc if shared is not None else None
        if args.scscore_model and shared is None:
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


def process_hit(args, engine, shared):
    started = time.perf_counter()
    sampling_seconds = filtering_seconds = 0.0
    calls = requested = 0
    history = defaultdict(lambda: [0,0])
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
    filters = LibraryFilter(args, valid, shared)
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
            attempted, accepted = history[design['design_id']]
            rate = (accepted + 2) / (attempted + 10)
            count = min(args.batch_size, max(args.min_samples, math.ceil(remaining / len(designs) / rate)))
            run_dir = args.out / 'runs' / f'{design["design_id"]}_round_{round_i:02d}'
            tick = time.perf_counter()
            raw = engine.sample(design, count, run_dir)
            sampling_seconds += time.perf_counter() - tick
            calls += 1
            requested += count
            before = len(generated)
            tick = time.perf_counter()
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
            filtering_seconds += time.perf_counter() - tick
            history[design['design_id']][0] += count
            history[design['design_id']][1] += len(generated) - before
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
    (args.out / 'timing.json').write_text(json.dumps(dict(sampling_seconds=sampling_seconds, filtering_seconds=filtering_seconds, total_seconds=time.perf_counter()-started, sampling_calls=calls, requested_samples=requested), indent=2))
    if len(rows) < args.num_molecules:
        print('WARNING: target not reached; try more rounds or additional designs.', file=sys.stderr)


def add_cumulative_scores(hit_dir):
    """Add static CMOS to every generated row; dynamic CMOS to selected rows.

    Does not affect sampling or acceptance. Preserves existing columns/order.
    """
    import argparse
    import importlib.util
    import pandas as pd

    csv_path = hit_dir / 'generated.csv'
    if not csv_path.is_file():
        return
    original = pd.read_csv(csv_path, low_memory=False)
    if 'IntrinsicScore' in original.columns and 'CumulativeScore' in original.columns:
        return
    if original.empty:
        original['IntrinsicScore'] = pd.Series(dtype=float)
        original['CumulativeScore'] = pd.Series(dtype=float)
        original.to_csv(csv_path, index=False)
        return

    scorer_path = Path(__file__).resolve().parent / 'cumulative_score.py'
    if not scorer_path.is_file():
        raise FileNotFoundError(f'CMOS scoring module missing: {scorer_path}')
    spec = importlib.util.spec_from_file_location('cumulative_score', scorer_path)
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)

    score_dir = hit_dir / 'cumulative_scoring'
    scorer.run(argparse.Namespace(
        input=str(csv_path), output=str(score_dir),
        weights=None, select_per_hit=30, missing='reject',
    ))
    scored = pd.read_csv(score_dir / 'scored_molecules.csv', low_memory=False)
    if len(scored) != len(original) or not scored['SMILES'].equals(original['SMILES']):
        raise RuntimeError(f'CMOS output does not match original rows: {csv_path}')
    score_columns = [c for c in scored.columns if c.startswith('Score_')]
    score_columns += ['ValidSMILES', 'MissingComponents', 'IntrinsicScore',
                      'CumulativeScore', 'SelectionRank', 'Selected']
    for column in score_columns:
        original[column] = scored[column].to_numpy()
    tmp = csv_path.with_suffix('.scored.tmp')
    original.to_csv(tmp, index=False)
    tmp.replace(csv_path)
    print(f'CMOS: added scores to {csv_path}', flush=True)


def main():
    args = parse_args()
    if min(args.num_molecules,args.max_rounds,args.batch_size,args.min_samples) < 1:
        raise ValueError('Counts must be positive')
    root = args.out
    if args.smiles or args.input.suffix.lower() == '.csv':
        parents = load_inputs(args)
    else:
        parents = []
        for line in args.input.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith('#'): continue
            parts = line.split()
            if len(parts) == 1: parents.append((f'Hit_{len(parents)+1:04d}',parts[0]))
            elif len(parts) == 2: parents.append(tuple(parts))
            else: raise ValueError('Expected SMILES or ID SMILES')
    if not parents: raise ValueError('No inputs')
    raw_parents = list(parents)
    audit = []
    standardized_parents = []
    for name, smi in raw_parents:
        try:
            standardized, record = standardize_parent_smiles(smi, args)
            standardized_parents.append((name, standardized))
            print(f'Input {name}: charge {record["OriginalFormalCharge"]:+d} -> {record["StandardizedFormalCharge"]:+d}; '
                  f'{smi} -> {standardized}', flush=True)
        except Exception as exc:
            record = dict(OriginalSMILES=smi, StandardizedSMILES='', OriginalFormalCharge='',
                          StandardizedFormalCharge='', OriginalFragments='', StandardizedFragments='',
                          ChargePolicy=args.input_charge_policy, LargestFragment=args.input_largest_fragment,
                          Changed='', Status='error', Error=str(exc))
            print(f'ERROR: invalid input {name}: {exc}', file=sys.stderr)
        audit.append(dict(ParentID=name, **record))
    if any(r['Status'] != 'ok' for r in audit):
        raise ValueError('Input standardization failed; correct invalid SMILES before generation')
    parents = standardized_parents
    settings = dict(vars(args)); settings.pop('resume'); settings = json.loads(json.dumps(settings,default=str))
    settings['parents'] = [list(x) for x in raw_parents]
    if root.exists() and any(root.iterdir()):
        if not args.resume or json.loads((root/'batch_settings.json').read_text()) != settings:
            raise ValueError('Existing output: resume requires identical settings')
    root.mkdir(parents=True,exist_ok=True)
    (root/'batch_settings.json').write_text(json.dumps(settings,indent=2))
    write_csv(root/'input_standardization.csv', audit,
              ['ParentID','OriginalSMILES','StandardizedSMILES','OriginalFormalCharge',
               'StandardizedFormalCharge','OriginalFragments','StandardizedFragments',
               'ChargePolicy','LargestFragment','Changed','Status','Error'])
    engine = PersistentSampler(args)
    shared = LibraryFilter(args, [])
    all_rows=[]; summary=[]; unique={}; fields=[]
    failed=False
    for i,(name,smi) in enumerate(parents,1):
        original_smi = raw_parents[i-1][1]
        folder=f'Hit_{i:04d}'
        dest=root/folder
        local=copy.copy(args); local.smiles=smi; local.input=None; local.out=dest
        rc=0; timing={}
        try:
            if not (args.resume and (dest/'timing.json').is_file() and (dest/'generated.csv').is_file()):
                if dest.exists():
                    # Preserve interrupted output rather than deleting it.
                    dest.rename(root/f'{folder}_interrupted_{time.time_ns()}')
                print(f'Processing {name}: {smi}',flush=True)
                process_hit(local,engine,shared)
            # CMOS: score the completed hit before it is merged into batch CSVs.
            # On --resume, score only if the columns are not already present.
            add_cumulative_scores(dest)
            timing=json.loads((dest/'timing.json').read_text())
        except (Exception,SystemExit) as exc:
            import traceback
            traceback.print_exc(); rc=1
        count=0
        if rc==0:
            with (dest/'generated.csv').open(newline='') as f:
                reader=csv.DictReader(f)
                for col in reader.fieldnames or []:
                    if col not in fields: fields.append(col)
                for row in reader:
                    key=Chem.MolToSmiles(Chem.MolFromSmiles(row['SMILES']),isomericSmiles=True)
                    count+=1
                    row.update(ParentID=name,ParentSMILES=original_smi,StandardizedParentSMILES=smi,ParentFolder=folder,GeneratedID=f'{folder}_{count:06d}')
                    all_rows.append(row); unique.setdefault(key,[]).append(row)
        status='failed' if rc else ('complete' if count>=args.num_molecules else 'shortfall')
        failed |= bool(rc)
        summary.append(dict(ParentID=name,ParentSMILES=original_smi,StandardizedParentSMILES=smi,ParentFolder=folder,Requested=args.num_molecules,Accepted=count,Status=status,ReturnCode=rc,**timing))
        meta=['ParentID','ParentSMILES','StandardizedParentSMILES','ParentFolder','GeneratedID']
        write_csv(root/'combined_generated.csv',all_rows,meta+fields)
        rows=[]
        for j,(key,members) in enumerate(unique.items(),1):
            row=dict(members[0],SMILES=key,LibraryID=f'LIB_{j:07d}',AllParentIDs=json.dumps(list(dict.fromkeys(x['ParentID'] for x in members))),AllParentFolders=json.dumps(list(dict.fromkeys(x['ParentFolder'] for x in members))),OccurrenceCount=len(members))
            rows.append(row)
        write_csv(root/'combined_unique_generated.csv',rows,['LibraryID','AllParentIDs','AllParentFolders','OccurrenceCount']+meta+fields)
        (root/'combined_unique_generated.smi').write_text(''.join(f"{r['SMILES']}\t{r['LibraryID']}\n" for r in rows))
        write_csv(root/'overall_summary.csv',summary,['ParentID','ParentSMILES','StandardizedParentSMILES','ParentFolder','Requested','Accepted','Status','ReturnCode','sampling_seconds','filtering_seconds','total_seconds','sampling_calls','requested_samples'])
    (root/'batch_timing.json').write_text(json.dumps(dict(model_load_seconds=engine.load_seconds,unique_molecules=len(unique),accepted_records=len(all_rows)),indent=2))
    if failed: sys.exit(1)

if __name__ == '__main__':
    main()
