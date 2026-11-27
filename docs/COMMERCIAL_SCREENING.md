# V7 Commercial Compound Screening

## Overview

V7 extends the molecular generation and multiparameter optimization
platform with commercial catalog screening through SwissSimilarity.

The module identifies exact molecular structures and structural analogs
reported in commercial compound libraries, connecting computational
molecular design with potential experimental acquisition.

The original V6 computational workflow and CMOS V2 scoring remain unchanged.

## Supported Commercial Libraries

- Enamine
- ChemBridge
- ChemDiv
- Life Chemicals
- Maybridge
- SPECS
- Asinex
- Otava

The default screening profile uses Enamine and ChemBridge.
The vendors profile includes all eight enabled libraries.

## Features

- Exact molecular identity verification using canonical SMILES and InChIKey.
- Stereochemistry-aware molecular matching.
- ECFP-based structural similarity screening.
- Configurable similarity thresholds and analog limits.
- Molecular query deduplication.
- Persistent SQLite screening cache.
- Interrupted-session recovery.
- Explicit remote submission authorization and submission limits.
- Offline cache-only screening.
- Commercial annotation of CMOS V2 selected molecules.
- Per-molecule and per-parent screening summaries.

## Match Classification

| Classification | Description |
|---|---|
| EXACT_MATCH | Identical molecular identity without unresolved stereochemistry |
| EXACT_REPORTED_STRUCTURE | Identical reported structure with unresolved stereochemistry |
| STEREOCHEMISTRY_MISMATCH | Matching connectivity but different stereochemical representation |
| FINGERPRINT_ONLY_MATCH | Identical fingerprint similarity without identical molecular identity |
| STRUCTURAL_ANALOG | Nonidentical structure returned by similarity screening |

## Full Pipeline

Run in the configured computational environment:

    bash scripts/run_reinvent_v7.sh full \
        examples/example_hits.txt \
        100 results/v7/example cns 10

## Commercial-Only Screening

Screen an existing CMOS V2 selected molecular library:

    bash scripts/run_reinvent_v7.sh commercial-only \
        results/example/cmos/CMOS_V2_selected.csv \
        results/v7/commercial_only \
        --libraries Enamine ChemBridge

## Remote Screening

Remote submissions are disabled by default.

Explicit authorization is required:

    bash scripts/run_reinvent_v7.sh commercial-only \
        tests/commercial_screening/mock_cmos.csv \
        results/v7/remote_example \
        --libraries Enamine \
        --allow-remote \
        --maximum-new-jobs 1

The submission limit applies to new jobs, not resumed sessions.

Only submit molecular structures authorized for transmission
to external services.

## Screening Options

| Option | Description |
|---|---|
| --libraries | Select specific commercial libraries |
| --profile | Select a predefined commercial library profile |
| --minimum-similarity | Minimum analog similarity, default 0.5 |
| --maximum-analogs | Maximum analogs per library, default 20 |
| --allow-remote | Authorize remote screening |
| --maximum-new-jobs | Limit new remote submissions |
| --skip-commercial | Run the original computational pipeline without screening |
| --mode both | Report exact matches and structural analogs |

The integrated V7 launcher requires screening mode both.
The standalone screening engine also supports exact and similarity modes.

## Configuration

Commercial settings are defined in:

    config/commercial_sources.json

The existing computational runtime configuration is reused:

    config/runtime.conf

Each screening run creates a separate commercial runtime JSON.
The base configuration remains unchanged.

## Outputs

Commercial screening:

    commercial/
        Commercial_Availability.csv
        exact_matches.csv
        commercial_analogs.csv
        special_matches.csv
        screening_jobs.csv
        invalid_queries.csv
        screening_summary.json

Integrated molecular library:

    final/
        CMOS_V2_Commercial_Annotated.csv
        CMOS_V2_Commercial_Exact.csv
        CMOS_V2_Commercial_Strict_Exact.csv
        CMOS_V2_Commercial_Analog_Queries.csv
        CMOS_V2_Commercial_Unscreened.csv
        Commercial_Per_Parent_Summary.csv
        Commercial_Screening_Summary.json

The annotated library preserves every original CMOS V2 column
and adds 21 commercial screening annotations.

Original optimization scores and selection rankings are unchanged.

## Persistent Cache

The SQLite database is stored locally:

    databases/commercial/screening_cache.sqlite

Completed screening results can be reused without additional
remote submissions. The database is excluded from version control.

## Validation

The original V7 implementation was validated on the University
of Florida HyperGator HPC infrastructure through:

- Molecular identity classification tests.
- Molecular query deduplication tests.
- Checkpoint and recovery tests.
- Offline commercial screening benchmarks.
- Controlled live commercial library screening.
- Full molecular generation and optimization benchmarks.
- Commercial integration and data integrity tests.

The repository includes synthetic fixtures for offline recovery testing.

Run within the configured computational environment:

    python tests/commercial_screening/test_checkpoint_recovery.py

The portable launcher has passed Bash syntax validation.
Independent end-to-end deployment remains to be validated.

## Limitations

- Catalog matches do not establish current stock availability.
- Commercial stock status remains UNVERIFIED.
- Similarity screening is not exhaustive catalog searching.
- Absence of an exact returned match does not establish catalog absence.
- Structural analogs do not necessarily retain biological activity.
- Generated molecular scores are not transferable to commercial analogs.
- Concurrent execution and large-scale remote throughput remain unvalidated.
- External service usage is subject to provider policies.

## External Service

SwissSimilarity:
https://www.swisssimilarity.ch/

Command-line interface:
https://www.swisssimilarity.ch/command-line.php
