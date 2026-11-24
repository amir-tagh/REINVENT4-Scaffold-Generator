# Installation Guide

## 1. System Requirements

- Linux x86-64
- Python 3.11
- Conda or Miniconda
- Git
- GNU Bash
- GNU Patch
- NVIDIA GPU with a compatible CUDA driver
- Recommended starting memory allocation: 32 GB
- Slurm (optional, for HPC execution)

The original V6 workflow was validated on Linux using
an NVIDIA GPU. CPU-only execution and other operating
systems have not been validated.

## 2. Clone the Repository

```bash
git clone https://github.com/amir-tagh/REINVENT4-Scaffold-Generator.git
cd REINVENT4-Scaffold-Generator
```

## 3. Create Computational Environments

The platform uses three independent Conda environments
to avoid dependency conflicts.

```bash
conda env create -f environments/reinvent4.yml
conda env create -f environments/admet_ai.yml
conda env create -f environments/cns_mpo.yml
```

| Environment | Purpose |
|---|---|
| reinvent4 | Molecular generation, cheminformatics and CMOS scoring |
| admet_ai | ADMET-AI and solubility prediction |
| cns_mpo | MolGpKa-based pKa prediction |

The supplied YAML files are dependency snapshots from
the working development environment.

Fresh environment creation has not yet been independently
validated. CUDA dependencies may require adjustment
for different GPU drivers and hardware.

## 4. Reference Software Versions

| Package | Version |
|---|---|
| Python | 3.11.18 |
| REINVENT4 | 4.8.24 |
| ADMET-AI | 2.0.1 |
| Chemprop | 2.3.1 |
| prop-profiler | 0.2.1 |
| PyTorch (generation) | 2.12.0 |
| PyTorch (prediction) | 2.15.0 |
| RDKit | 2026.9.1 / 2026.9.2 |
| PyTorch Geometric | 2.8.1 |

## 5. Install External Models

### 5.1. LibInvent Pretrained Model

Obtain a compatible LibInvent prior from the official
REINVENT4 project:

https://github.com/MolecularAI/REINVENT4

Create the model directory:

```bash
mkdir -p priors
```

Place the compatible pretrained model at:

```text
priors/libinvent.prior
```

Model weights are not distributed with this repository.

### 5.2. SCScore

Clone the official SCScore repository:

```bash
git clone https://github.com/connorcoley/scscore.git
```

The default expected model location is:

```text
scscore/models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz
```

Verify that this model exists. Obtain any missing model
resources from the upstream project.

### 5.3. ADMET-AI and MolGpKa

ADMET-AI and prop-profiler are included in their
respective environment specifications.

Additional model resources may be required by the
upstream packages during initialization.

Consult their official documentation if model
initialization fails.

## 6. LibInvent CUDA Compatibility Patch

The validated REINVENT4 installation contains a local
device-alignment modification to sample_decorations().

The patch initializes sampling tensors on the same
device as the neural network.

Activate the generation environment:

```bash
conda activate reinvent4
```

Locate the installed package directory:

```bash
SITE_PACKAGES=$(python -c 'import reinvent; from pathlib import Path; print(Path(reinvent.__file__).resolve().parent.parent)')
```

Check patch compatibility:

```bash
patch --dry-run -p1 -d "$SITE_PACKAGES" \
    < patches/libinvent_cuda_device.patch
```

Apply only if the dry run succeeds:

```bash
patch -p1 -d "$SITE_PACKAGES" \
    < patches/libinvent_cuda_device.patch
```

This patch was derived from REINVENT4 version 4.8.24.

Do not force the patch onto incompatible versions or
apply it repeatedly to an already modified installation.

## 7. Runtime Configuration

```bash
cp config/runtime.example.conf config/runtime.conf
```

Edit the configuration if necessary:

```bash
nano config/runtime.conf
```

Available settings include:

| Parameter | Description |
|---|---|
| REINVENT_ENV | Generation environment name |
| ADMET_ENV | ADMET prediction environment name |
| PKA_ENV | pKa prediction environment name |
| REINVENT_PY | Optional explicit Python executable |
| ADMET_PY | Optional explicit Python executable |
| PKA_PY | Optional explicit Python executable |
| LIBINVENT_PRIOR | LibInvent pretrained model |
| SCSCORE_REPO | SCScore repository location |
| SCSCORE_MODEL | SCScore pretrained model |
| REINVENT_DEVICE | Molecular generation device |
| PKA_DEVICE | pKa prediction device |
| CPU_THREADS | Number of CPU threads |

The default configuration uses CUDA-enabled generation
and pKa prediction.

## 8. Verify the Installation

```bash
conda run -n reinvent4 python -c \
    "import torch, rdkit, reinvent; print(torch.__version__, rdkit.__version__); print('CUDA:', torch.cuda.is_available())"

conda run -n admet_ai python -c \
    "from admet_ai import ADMETModel; print('ADMET-AI import successful')"

conda run -n cns_mpo python -m pip show prop-profiler
```

Check external model files:

```bash
ls -lh priors/libinvent.prior

ls -lh scscore/models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz
```

Check launcher syntax:

```bash
bash -n scripts/run_reinvent_v6.sh
```

These checks verify basic installation requirements.
They do not replace an end-to-end computational test.

## 9. Prepare Input Molecules

The input is a whitespace-delimited text file containing
a SMILES string followed by an optional compound identifier.

Example:

```text
CC(=O)Oc1ccccc1C(=O)O Example_001
CC(C)Cc1ccc(C(C)C(=O)O)cc1 Example_002
```

Use unique identifiers and valid molecular SMILES.

## 10. Local GPU Execution

Run from the repository root:

```bash
bash scripts/run_reinvent_v6.sh \
    examples/example_hits.txt \
    100 \
    results/example \
    cns \
    10
```

Positional arguments:

| Argument | Description | Default |
|---|---|---|
| 1 | Input molecular library | Required |
| 2 | Target molecules per parent | 1000 |
| 3 | Output directory | results/v6/REINVENT_V6_Results |
| 4 | Optimization profile | cns |
| 5 | Selected molecules per parent | 30 |

## 11. HPC Execution with Slurm

Submit from the repository root:

```bash
sbatch \
    --account=YOUR_ACCOUNT \
    --partition=YOUR_GPU_PARTITION \
    scripts/run_reinvent_v6.sh \
    examples/example_hits.txt \
    1000 \
    results/example_hpc \
    cns \
    30
```

Adjust account, partition, QoS, memory, GPU and CPU
resources according to your computing facility.

## 12. Optimization Profiles

Available profiles include:

- cns
- cns_safety
- general
- oral
- experimental

Profile weights are defined in:

```text
config/optimization_profiles.json
```

See [CMOS Documentation](CMOS_DOCUMENTATION.md)
for the scoring methodology.

## 13. Output Files

```text
results/example/
├── generation/
│   ├── combined_generated.csv
│   ├── combined_unique_generated.csv
│   └── combined_unique_generated.smi
├── admet/
│   ├── ADMET_Predictions.csv
│   ├── ADMET_pKa.csv
│   ├── ADMET_CNS.csv
│   ├── ADMET_CNS_ADME.csv
│   └── ADMET_Complete.csv
├── cmos/
│   ├── CMOS_V2.csv
│   ├── CMOS_V2_selected.csv
│   ├── CMOS_V2_selected.smi
│   └── CMOS_V2.summary.json
├── optimization_profiles.json
└── Final_Property_Comparison.csv
```

The selected SMILES file contains the prioritized
molecular library for downstream computational and
experimental evaluation.

## 14. Troubleshooting

### CUDA device mismatch

Verify GPU availability and check whether the
LibInvent compatibility patch is required.

### Missing Python dependencies

Check that the correct Conda environment is being used.

### Missing pretrained models

Verify the model locations in config/runtime.conf.

### Insufficient generated molecules

Some scaffolds may not produce the requested number
of unique molecules under the specified filters.

### Environment installation conflicts

Review the exported dependency versions and the
installation requirements of the upstream packages.

## 15. Limitations

- The portable launcher has passed syntax checks but
  has not yet undergone independent end-to-end testing.
- Fresh installation from the YAML files remains untested.
- Predicted ADMET properties require experimental validation.
- Predicted pKa is excluded from active CMOS optimization.
- RNA binding affinity is not currently an optimization objective.
- Commercial availability is not verified.
- Third-party software and pretrained models remain subject
  to their respective licenses.
