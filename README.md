# Automated scaffold-to-LibInvent workflow

Requires the existing REINVENT4 installation and RDKit environment. Run from the REINVENT4 repository directory so relative model paths resolve.

## Examples

```bash
python scaffold_libinvent.py --smiles 'Cc1ccccc1' --num-molecules 100 --model priors/libinvent.prior --out results_toluene
python scaffold_libinvent.py --input hits.smi --num-molecules 500 --out results_hits
python scaffold_libinvent.py --input hits.csv --smiles-column SMILES --num-molecules 1000 --out results_hits
```

The number requested is the **target number of unique valid generated molecules** across all selected scaffold designs. The pipeline may return fewer if the model cannot produce enough valid, unique structures within `--max-rounds`.

Input files: `.smi`/`.txt` (SMILES first column, optional identifier second column), `.csv` (column default `SMILES`).

Outputs:
- `inputs.csv`: validated input molecules and their Bemis–Murcko scaffolds
- `candidates.csv`: ranked candidate scaffold positions and their rationale
- `designs.csv`: selected two-point LibInvent designs and mapped scaffold SMILES
- `generated.csv`: unique, valid output structures with input/design provenance and RDKit properties
- `molecules.png`: 2D grid of example generated molecules (if Pillow is installed)
- `runs/`: TOML configs, scaffold inputs, raw outputs and logs for each sampling round

## How attachment points are selected

RDKit extracts a Bemis–Murcko scaffold. Candidate sites are atoms on that core with an available hydrogen (or, optionally, atoms bearing substituents in the input compound). Substituents attached to the core are **not preserved** in the LibInvent input; they are considered variable R-groups. Two attachment sites are selected per scaffold design. Candidate sites are ranked heuristically by frequency of attachment in the input set, whether a site is an existing substituent position, and a modest preference for carbon sites. This is **not** a prediction of binding affinity, synthetic feasibility, or pocket accessibility.

If a core has fewer than two eligible positions, the input is reported as skipped. When no shared scaffold exists across inputs, each distinct scaffold is processed separately. A single shared scaffold is not forced across unrelated compounds.

## Caveats

- The LibInvent prior in this workflow expects **two attachment points** per scaffold.
- Bemis–Murcko scaffolds may be too minimal for a medicinal-chemistry series; inspect `designs.csv` before trusting generated structures.
- No RNA/protein structure or docking is used for attachment-site selection or scoring.
- A generation count is a *target*, not a guarantee.
- The pipeline assumes the installed REINVENT4 CLI accepts `reinvent -l LOG CONFIG`, as tested in your environment.
