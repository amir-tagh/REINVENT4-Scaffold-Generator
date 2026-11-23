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
