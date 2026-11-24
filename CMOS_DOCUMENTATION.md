# Cumulative Molecular Optimization Score (CMOS)

## Overview

CMOS is a docking-independent, multi-objective molecular scoring and diversity-aware selection framework integrated into the REINVENT4 LibInvent scaffold-based molecular generation pipeline.

The scoring module runs automatically after generation of each parent hit.

## Scoring Components

| Component | Default weight |
|---|---:|
| Quantitative Estimate of Drug-likeness (QED) | 0.25 |
| Synthetic Accessibility (SA) | 0.20 |
| Synthetic Complexity Score (SCScore) | 0.15 |
| Parent-molecule similarity | 0.15 |
| Chemical diversity | 0.15 |
| Molecular properties | 0.10 |

**IntrinsicScore (0–100):** Combines the molecular quality components, excluding diversity, with the remaining weights renormalized.

**CumulativeScore (0–100):** Includes chemical diversity relative to previously selected molecules using a greedy selection strategy.

By default, CMOS selects up to 30 compounds per parent hit.

## Pipeline Integration

The main generator, `scaffold_libinvent_v5_devicefix.py`, automatically invokes `cumulative_score.py` after each parent hit is processed.

Scoring results are merged into each hit's `generated.csv` and subsequently into the combined generation outputs.

Additional CSV columns:

- `Score_QED`
- `Score_SA`
- `Score_SCScore`
- `Score_ParentSimilarity`
- `Score_Properties`
- `Score_Diversity`
- `ValidSMILES`
- `MissingComponents`
- `IntrinsicScore`
- `CumulativeScore`
- `SelectionRank`
- `Selected`

All eligible molecules receive an intrinsic score. Only selected molecules receive a cumulative score and selection rank.

## HyperGator Execution

The pipeline is executed in the configured REINVENT4 Conda environment using GPU-enabled Slurm nodes.

From the HyperGator REINVENT4 project root:

```bash
sbatch scripts/run_multiple_hits_v4_standardized.sh \
    inputs/Tau-hits.txt \
    1000 \
    results/v5/REINVENT_V5_All_8_Hits_CMOS
```

The HyperGator project stores the custom scripts in `scripts/`. This GitHub repository stores them in its root directory.

## Validation

The integrated pipeline was validated using the first MAPT RNA parent hit.

| Metric | Result |
|---|---:|
| Generated molecules | 1,000 |
| Molecules with IntrinsicScore | 1,000 |
| Molecules with CumulativeScore | 30 |
| CMOS-selected compounds | 30 |

The scoring module successfully updated `Hit_0001/generated.csv` with molecular scores and selection ranks.

## Current Limitations

- CMOS currently performs post-generation scoring and ranking, not reinforcement-learning-guided generation.
- Docking scores are not included in the current CMOS formulation.
- The validated HyperGator Conda environment includes a CUDA device-alignment fix applied to the installed LibInvent package. This environment-level modification is not included in the repository.
- Molecular scoring and selection do not establish RNA binding affinity, biological activity, or experimental selectivity.

---

## CMOS V2 — Multiparameter Molecular Optimization

CMOS V2 extends the original scoring framework with
solubility, CNS, BBB, ADME and toxicity predictions.

### Default CNS Optimization Weights

| Component | Weight |
|---|---:|
| QED | 0.10 |
| Synthetic accessibility | 0.10 |
| Synthetic complexity | 0.10 |
| Parent similarity | 0.10 |
| Structural diversity | 0.15 |
| Consensus solubility | 0.10 |
| CNS desirability | 0.10 |
| BBB desirability | 0.05 |
| ADME desirability | 0.05 |
| Toxicity desirability | 0.15 |
| Total | 1.00 |

### Methodology

- Normalized component desirability scores between 0 and 1.
- Weighted intrinsic molecular quality scoring.
- Consensus ESOL and ADMET-AI solubility prediction.
- Five-component CNS desirability scoring.
- Multiparameter ADME and toxicity assessment.
- Morgan fingerprint-based structural diversity.
- Greedy diversity-aware compound selection per parent.
- Configurable optimization profiles.
- Missing-component eligibility checks.

Predicted pKa is reported but excluded from active
optimization pending further model validation.

### Implementation

- scripts/cumulative_score_v2.py
- scripts/molecular_admet.py
- scripts/molecular_pka.py
- scripts/molecular_cns.py
- scripts/molecular_adme.py
- scripts/molecular_toxicity.py
- config/optimization_profiles.json
- config/cns_sensitivity_profiles.json

### Pipeline Integration

The V6 pipeline uses CMOS V2 for final compound
prioritization through scripts/run_reinvent_v6.sh.

CMOS V1 remains available for backward compatibility
and methodological benchmarking.

The maintained implementation is located in scripts/.

CMOS V2 scores range from 0 to 1 and are not directly
comparable with CMOS V1 scores.

The workflow has been computationally validated on
HyperGator. Experimental validation and portable
environment configuration remain necessary.
